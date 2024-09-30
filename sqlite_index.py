import xml.etree.ElementTree as ET
import os
import re
import sqlite3
from nltk.tokenize import word_tokenize,sent_tokenize

def parse_pubmed_xml(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    pmid = root.find(".//PMID").text if root.find(".//PMID") is not None else "Unknown PMID"
    title = root.find(".//ArticleTitle").text if root.find(".//ArticleTitle") is not None else "No Title"
    pubdate = extract_pubdate(root)
    authors = extract_authors(root)


    abstract_list = []
    for elem in root.iter('AbstractText'):
        if elem is not None:
            elem_text = ''.join(elem.itertext())
        abstract_list.append(elem_text)
    abstract = "\n".join(abstract_list)
    
    
    return {"pmid": pmid, "title": title, "pubdate":pubdate, "authors": authors,"abstract": abstract}

def check_pubmed(xml_file,db_name):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    title = root.find(".//ArticleTitle").text if root.find(".//ArticleTitle") is not None else "No Title"

    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    c.execute("SELECT 1 FROM pubmed_articles WHERE title = ?", (title,))
    result = c.fetchone()
    conn.close()
    # Return True if the title exists, False otherwise
    return result is not None

def extract_authors(root):
    authors = []
    for author in root.findall(".//Author"):
        last_name = author.find("LastName").text
        fore_name = author.find("ForeName").text 
        full_name = f"{fore_name} {last_name} "
        authors.append(full_name)
        authors_string = ",".join(authors)
    return authors_string

def extract_pubdate(root):
    pubdate_element = root.find(".//PubDate")
    # Check if PubDate exists
    if pubdate_element is not None:
        year = pubdate_element.find("Year").text if pubdate_element.find("Year") is not None else "Unknown Year"
        month = pubdate_element.find("Month").text if pubdate_element.find("Month") is not None else "Unknown Month"
        pubdate = f"{year}-{month}"
    else:
        pubdate = "Unknown Published Date"
    return pubdate

def parse_all_xml_files(xml_dir,db_name):
    articles = []
    for file_name in os.listdir(xml_dir):
        if file_name.endswith(".xml"):
            file_path = os.path.join(xml_dir, file_name)
            if not check_pubmed(file_path,db_name): 
                articles.append(parse_pubmed_xml(file_path))
                print(f'{file_name} is added to data')
    return articles

# Function to create a SQLite database and table
def create_database(db_name):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    c.execute('''
        CREATE TABLE IF NOT EXISTS pubmed_articles (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            pmid TEXT UNIQUE,
            title TEXT,
            pubdate TEXT,
            authors TEXT,
            abstract TEXT
        )
    ''')
    conn.commit()
    conn.close()


# Function to insert parsed articles into the SQLite database
def insert_articles_to_db(db_name, articles):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    for article in articles:
        c.execute('''
            INSERT OR IGNORE INTO pubmed_articles (pmid, title, pubdate,authors, abstract)
            VALUES (?, ?, ?, ?, ?)
        ''', (article['pmid'], article['title'], article['pubdate'],article['authors'], article['abstract']))
    conn.commit()
    conn.close()

# Function to search the SQLite database
def search_db_by_query(db_name, query,selected_year):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    
    keywords = query.split()
    keyword_conditions = ' OR '.join(['title LIKE ? OR abstract LIKE ?' for _ in keywords])
    sql_query = f'''
            SELECT id,pmid, title, pubdate,abstract
            FROM pubmed_articles
            WHERE ({keyword_conditions})
            AND CAST(SUBSTR(pubdate, 1, 4) AS INTEGER) >= ?
        '''
    # Build the parameter list, each keyword applies to both title and abstract
    parameters = []
    for keyword in keywords:
        parameters.extend(['%' + keyword + '%', '%' + keyword + '%'])
    parameters.append(selected_year)
    c.execute(sql_query, parameters)   

    results = c.fetchall()
    conn.close()

    # Debugging step: check if `pubdate` was formatted correctly
    for result in results:
        print(f"Pubdate: {result[3]}")  # Assuming `pubdate` is the 4th column

    return results

import streamlit as st

def highlight_text(text, query):

    # Highlight search query in the text.
    keywords = query.split()
    pattern = '|'.join([re.escape(query) for query in keywords])
    # Function to wrap matched text in <mark> tags
    def highlight(match):
        return f'<mark>{match.group(0)}</mark>'
    # Substitute matches with highlighted text
    highlighted_text = re.sub(pattern, highlight, text, flags=re.IGNORECASE)

    
    # highlighted_text = re.sub(
    #     re.escape(query),
    #     lambda match: f'<mark>{match.group(0)}</mark>',
    #     text,
    #     flags=re.IGNORECASE
    # )

    return highlighted_text
    

# Define CSS to set background color
page_bg_css = '''
<style>
    [data-testid="stAppViewContainer"] {
        background-color: #E3F2FD;  /* Light blue background color */
    }
    mark {
        background-color: yellow;  /* Highlight color */
    }
</style>
'''

def main_page(db_name):

    
    # initialize home page
    if 'page' not in st.session_state:
        st.session_state.page = 'home'

    print(f"session state: {st.session_state.page}")

    # initialize search_clicked
    if 'search_clicked' not in st.session_state:
        st.session_state.search_clicked = False
    # sava current query
    if 'query' not in st.session_state:
        st.session_state.query = ""


    # If page is "main" show the main page
    if st.session_state.page == 'home': 
        st.markdown(page_bg_css, unsafe_allow_html=True)

        st.title("PubMed Search Engine")
        col1, col2 = st.columns([3, 1]) 

        with col1:
            st.session_state.query = st.text_input("Enter search query", value=st.session_state.query,placeholder="Search...")

        # Search button
        with col2:
            if st.button("Search"):
                st.session_state.search_clicked = True

        # Layout 
        layout_type = st.selectbox("Choose layout type", ["Vertical", "Horizontal"])
                                   
        # Uploaded function
        uploaded_file = st.file_uploader("Choose a file", type=["xml"])
        if uploaded_file is not None:
          
            file_name = uploaded_file.name
            file_path = os.path.join("./data/", file_name)

            # Save the file to the local "uploads" folder
            with open(file_path, "wb") as f:
                f.write(uploaded_file.getbuffer())

            st.success(f"File '{file_name}' uploaded and saved to local folder successfully!")
        
        # Filter function
        selected_year = st.selectbox("Show results including and after:", get_available_years())

        # Search function
        if st.session_state.search_clicked:
            if st.session_state.query:
                results = search_db_by_query(db_name,st.session_state.query,selected_year)  

                # if horizontal
                if layout_type == "Horizontal":
                    display_horizontal(results)

                else:                       
                    count = 0
                    for id, pmid, title,pubdate, abstract in results:
                        # Filter
                        count += 1
                        # save selected result id
                        display(id, pmid, title,pubdate, abstract)
                        
                        st.write("---")
                    st.write(f"{count} results found ")   
                
                
    # show detail page
    else:
        detail_page(str(st.session_state.selected_id))
        
# Function to fetch available years from the database
def get_available_years():
    conn = sqlite3.connect('pubmed.db')
    cursor = conn.cursor()

    # Query to get distinct years
    cursor.execute('''
    SELECT DISTINCT CAST(substr(pubdate, 1, 4) AS INTEGER) AS year
    FROM pubmed_articles
    ORDER BY year;
    ''')
    years = [row[0] for row in cursor.fetchall()]

    conn.close()
    return years

# Function to filter results based on the selected year
def filter_by_year(selected_year):
    conn = sqlite3.connect('pubmed.db')
    cursor = conn.cursor()

    # Query to filter by the selected year
    cursor.execute('''
    SELECT *
    FROM pubmed_articles
    WHERE CAST(substr(pubdate, 1, 4) AS INTEGER) = ?;
    ''', (selected_year,))
    
    results = cursor.fetchall()

    conn.close()
    return results

# Function to handle page navigation
def go_page(page_name):
        st.session_state.page = page_name
        st.rerun()  # Force re-run to reflect the state change immediately

def display_horizontal(results):
    col1, col2 = st.columns(2)  # Create two columns
    count = 0
    
    for id, pmid, title,pubdate, abstract in results:
        if count % 2 == 0:
            with col1:
                display(id, pmid, title,pubdate, abstract)
                st.write("---")
        else:
            with col2:
                display(id, pmid, title,pubdate, abstract)
                st.write("---")
        count += 1


def display(id, pmid, title,pubdate, abstract):
        st.session_state.selected_id = id
        st.write(f"**PMID:** {pmid}")
        st.write(f"**Pubdate:** {pubdate}")
        highlighted_title = highlight_text(title, st.session_state.query)
        highlighted_abstract = highlight_text(abstract, st.session_state.query)
        st.markdown(f"**Title:** {highlighted_title}", unsafe_allow_html=True)
        st.markdown(f"**Abstract:** {limit_words(highlighted_abstract,30)}", unsafe_allow_html=True)
        if st.button(f"More",key= id):
            go_page(str(id))


def detail_page(id):
    
    db_id, pmid, title, pubdate, authors, abstract = search_db_by_id(id)
    highlighted_title = highlight_text(title, st.session_state.query)
    highlighted_abstract = highlight_text(abstract, st.session_state.query)
    st.markdown(page_bg_css, unsafe_allow_html=True)
    st.title(f"**More about**")
    st.write(f"**Database ID:** {db_id}")
    st.write(f"**PMID:** {pmid}")
    st.markdown(f"**Title: {highlighted_title}**", unsafe_allow_html=True)
    st.write(f"**Publication Date:** {pubdate}")
    st.write(f"**Authors:** {authors}")
    st.markdown(f"**Abstract:** {highlighted_abstract,50}", unsafe_allow_html=True)


    with st.expander(f"Statistic", expanded=False):
        st.write(f"number of sentences: {count_sentences(abstract)}")
        st.write(f"number of words: {count_words(abstract)}")
        st.write(f"number of characters: {count_characters(abstract)}")
        st.write(f"number of ascii: {count_ascii_and_non_ascii(abstract)[0]}")
        st.write(f"number of non-ascii: {count_ascii_and_non_ascii(abstract)[1]}")

    if st.button("Back to Results"):
        go_page("home")
    
def search_db_by_id(id):
    conn = sqlite3.connect('pubmed.db')
    c = conn.cursor()
    
    # Query to fetch the article with the given id
    c.execute("SELECT id, pmid, title, pubdate, authors, abstract FROM pubmed_articles WHERE id = ?", (id))
    result = c.fetchone()  # Fetch one result
    conn.close()
    
    return result

def limit_words(text, max_words):
    words = text.split()
    if len(words) > max_words:
        return ' '.join(words[:max_words]) + '...'
    return text

def count_sentences(text):
    # Tokenize the text into sentences
    sentences = sent_tokenize(text)
    return len(sentences)
    
def count_characters(text):    
    return len(text)

def count_words(text):
    # Split with space
    words_list = text.split()
    return len(words_list)

def count_ascii_and_non_ascii(text):
    ascii_count = 0
    non_ascii_count = 0
    
    for char in text:
        if ord(char) < 128:
            ascii_count += 1
        else:
            non_ascii_count += 1
            
    return ascii_count, non_ascii_count



if __name__ == "__main__":

    xml_directory = "./data/"  
    db_name = "pubmed.db"
    create_database(db_name)
    articles = parse_all_xml_files(xml_directory,db_name)
    insert_articles_to_db(db_name, articles)
    main_page(db_name)
    



import xml.etree.ElementTree as ET
import os
import re
import sqlite3
import nltk
from nltk.tokenize import word_tokenize,sent_tokenize



def parse_pubmed_xml(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    pmid = root.find(".//PMID").text if root.find(".//PMID") is not None else "Unknown PMID"
    title = root.find(".//ArticleTitle").text if root.find(".//ArticleTitle") is not None else "No Title"
    authors = extract_authors(root)
    abstract = root.find(".//AbstractText").text if root.find(".//AbstractText") is not None else "No Abstract"
    return {"pmid": pmid, "title": title, "authors": authors,"abstract": abstract}


def extract_authors(root):
    authors = []
    for author in root.findall(".//Author"):
        last_name = author.find("LastName").text
        fore_name = author.find("ForeName").text
        # You can choose how to format the author name 
        full_name = f"{fore_name} {last_name} "
        authors.append(full_name)
        authors_string = ",".join(authors)
    
    return authors_string


def parse_all_xml_files(xml_dir):
    articles = []
    for file_name in os.listdir(xml_dir):
        if file_name.endswith(".xml"):
            file_path = os.path.join(xml_dir, file_name)
            articles.append(parse_pubmed_xml(file_path))
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
            authors TEXT,
            abstract TEXT
        )
    ''')
    conn.commit()
    conn.close()
    print("database created")

# Function to insert parsed articles into the SQLite database
def insert_articles_to_db(db_name, articles):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    for article in articles:
        c.execute('''
            INSERT OR IGNORE INTO pubmed_articles (pmid, title, authors, abstract)
            VALUES (?, ?, ?, ?)
        ''', (article['pmid'], article['title'], article['authors'], article['abstract']))
    conn.commit()
    conn.close()

# Function to search the SQLite database
def search_db_by_query(db_name, query):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    
    keywords = query.split()
    keyword_conditions = ' OR '.join(['title LIKE ? OR abstract LIKE ?' for _ in keywords])
    sql_query = f'''
            SELECT id,pmid, title, abstract FROM pubmed_articles
            WHERE {keyword_conditions}
        '''
    # Build the parameter list, each keyword applies to both title and abstract
    parameters = []
    for keyword in keywords:
        parameters.extend(['%' + keyword + '%', '%' + keyword + '%'])
    
    c.execute(sql_query, parameters)   

    
    results = c.fetchall()
    conn.close()
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
            
            
        if st.session_state.search_clicked:
            if st.session_state.query:

                
                results = search_db_by_query(db_name,st.session_state.query)   
        
                st.write(f"{len(results)} results ")   
                for id, pmid, title, abstract in results:

                    # save selected result id
                    st.session_state.selected_id = id
                    st.write(f"**PMID:** {pmid}")
                    
                    highlighted_title = highlight_text(title, st.session_state.query)
                    highlighted_abstract = highlight_text(abstract, st.session_state.query)
                    st.markdown(f"**Title:** {highlighted_title}", unsafe_allow_html=True)
                    st.markdown(f"**Abstract:** {limit_words(highlighted_abstract,50)}", unsafe_allow_html=True)

                    if st.button(f"More",key= id):
                        go_page(str(id))
                    
                    st.write("---")
    # show detail page
    else:
        detail_page(str(st.session_state.selected_id))
        

# Function to handle page navigation
def go_page(page_name):
        st.session_state.page = page_name
        st.rerun()  # Force re-run to reflect the state change immediately


def detail_page(id):
    
    db_id, pmid, title, authors, abstract = search_db_by_id(id)

    st.markdown(page_bg_css, unsafe_allow_html=True)
    st.title(f"**More about**")
    st.write(f"**Database ID:** {db_id}")
    st.write(f"**PMID:** {pmid}")
    st.write(f"**Title:** {title}")
    st.write(f"**Authors:** {authors}")
    st.write(f"**Abstract:** {abstract}")

    with st.expander(f"Statistic", expanded=False):
        st.write(f"number of sentences: {count_sentences(abstract)}")
        st.write(f"number of words: {count_words(abstract)}")
        st.write(f"number of characters: {count_characters(abstract)}")

    if st.button("Back to Results"):
        go_page("home")
    
def search_db_by_id(id):
    conn = sqlite3.connect('pubmed.db')
    c = conn.cursor()
    
    # Query to fetch the article with the given id
    c.execute("SELECT id, pmid, title, authors, abstract FROM pubmed_articles WHERE id = ?", (id))
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

# nltk.download('punkt_tab')
def count_words(text):
    # Split with nltk
    # words_with_punctuation = word_tokenize(text)
    # words_list = [word for word in words_with_punctuation if re.match(r'^[a-zA-Z\-]+$', word)]

    # Split with space
    words_list = text.split()
    # print(words_list)
    return len(words_list)

if __name__ == "__main__":

    # xml_directory = "./data/obesity/xml"  
    # articles = parse_all_xml_files(xml_directory)
    db_name = "pubmed.db"
    # create_database(db_name)
    # insert_articles_to_db(db_name, articles)
    main_page(db_name)
    



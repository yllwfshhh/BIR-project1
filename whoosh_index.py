from whoosh.index import create_in,open_dir,exists_in
from whoosh.fields import Schema, TEXT, ID
from whoosh.qparser import QueryParser
import xml.etree.ElementTree as ET
import os
import shutil
import re




def parse_pubmed_xml(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()

    article_info = {}
    article_info['pmid'] = root.find(".//PMID").text
    article_info['title'] = root.find(".//ArticleTitle").text
    article_info['authors'] = extract_authors(root)
    abstract = root.find(".//AbstractText")
    article_info['abstract'] = abstract.text if abstract is not None else "No Abstract"
    return article_info


def extract_authors(root):
    authors = []
    for author in root.findall(".//Author"):
        last_name = author.find("LastName").text
        fore_name = author.find("ForeName").text
        # You can choose how to format the author name 
        full_name = f"{fore_name} {last_name} "
        authors.append(full_name)
    
    return authors


def parse_all_xml_files(xml_dir):
    articles = []
    for file_name in os.listdir(xml_dir):
        
        if file_name.endswith(".xml"):
            file_path = os.path.join(xml_dir, file_name)
            articles.append(parse_pubmed_xml(file_path))
    return articles


def index_articles(articles,index):
    
    with index.writer(clear=True) as writer:
        for article in articles:

            # print(f"pmif:{article['pmid']}")
            # print(f"title:{article['title']}")
            # print(f"author:{article['authors']}")
            

            writer.add_document(
                pmid=article['pmid'],
                title=article['title'],
                abstract=article['abstract'],
                authors=article['authors']
                
            )
            
        print("Indexing complete.")
        

def search(query_str):
    with index.searcher() as searcher:
        query = QueryParser("abstract", index.schema).parse(query_str)
        results = searcher.search(query, limit=10)  # Return top 10 results
        
        # Collect results before the searcher is closed
        search_results = []
        for result in results:
            search_results.append({
                "pmid": result["pmid"],
                "title": result["title"],
                "authors": result["authors"],
                "abstract": result["abstract"]
            })
    return search_results

import streamlit as st

def highlight_text(text, query):
    """Highlight search query in the text."""
    highlighted_text = re.sub(
        re.escape(query),
        lambda match: f'<mark>{match.group(0)}</mark>',
        text,
        flags=re.IGNORECASE
    )
    return highlighted_text

def streamlit_interface():

    st.title("PubMed Search Engine")
    query = st.text_input("Enter search query")
    print(f"query:{query}") 

    if st.button("Search"):
        print("search button clicked")
        if query:
            print(f"query:{query}") 
            results = search(query)

            st.write(f"{len(results)} results ")   
            for result in results:
                st.write(f"**PMID:** {result['pmid']}")
                
                # Highlight search terms in both title and abstract
                highlighted_title = highlight_text(result['title'], query)
                highlighted_abstract = highlight_text(result['abstract'], query)
                st.markdown(f"**Title:** {highlighted_title}", unsafe_allow_html=True)
                st.write(f"**Authors:** {result['authors']}")
                st.markdown(f"**Abstract:** {highlighted_abstract}", unsafe_allow_html=True)
                st.write("---")

if __name__ == "__main__":

    if not os.path.exists("indexdir"):

        print("index didn't exist, create an index")
        # Define the schema for the index
        schema = Schema(pmid=ID(stored=True), title=TEXT(stored=True), authors=TEXT(stored=True),abstract=TEXT(stored=True))

        # Parsing
        print("parsing")
        articles = parse_all_xml_files("./data/obesity/xml")

        # Recreate the index directory and schema
        print("mkdir indexdir")
        os.mkdir("indexdir")
        index = create_in("indexdir", schema)

        # Index
        print("indexing")
        index_articles(articles,index)

    else: 
        print("Index already exists, skipping indexing.")
        index = open_dir("indexdir")

    streamlit_interface()

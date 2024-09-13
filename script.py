import streamlit as st
from bs4 import BeautifulSoup
import re


# Function to parse the XML and extract text content
def parse_xml(file):
    soup = BeautifulSoup(file, "xml")
    articles = soup.find_all("PubmedArticle")
    print(articles)
    
    data = []
    
    for article in articles:
        title = article.find("ArticleTitle").get_text() if article.find("ArticleTitle") else ""
        abstract = article.find("AbstractText").get_text() if article.find("AbstractText") else ""
        authors = [author.get_text() for author in article.find_all("Author")]
        keywords = [kw.get_text() for kw in article.find_all("Keyword")]
        
        # Save document information
        data.append({
            "title": title,
            "abstract": abstract,
            "authors": authors,
            "keywords": keywords

        })
    
    return data

# Function to count sentences in the text
def count_sentences(text):
    sentences = re.split(r'(?<!\w\.\w.)(?<![A-Z][a-z]\.)(?<=\.|\?)\s', text)
    return len(sentences)

# Function to search documents by keywords
def search_documents(data, keyword):
    results = []
    for doc in data:
        if keyword.lower() in doc["title"].lower() or keyword.lower() in doc["abstract"].lower():
            results.append(doc)
    return results

def main():

    # Upload the XML file
    uploaded_file = st.file_uploader("Upload an XML file", type="xml")


    if uploaded_file is not None:
        st.write("File uploaded successfully!")
        
        # Parse the XML file
        docs = parse_xml(uploaded_file)
        
        # Show stats for each document
        for i, doc in enumerate(docs):
            st.write(f"### Document {i+1}")
            st.write(f"**Title:** {doc['title']}")
            st.write(f"**Keywords:** {', '.join(doc['keywords'])}")
            st.write(f"**Authors:** {', '.join(doc['authors'])}")
            st.write(f"**Sentence Count (Abstract):** {count_sentences(doc['abstract'])}")

        
        # Search functionality
        keyword = st.text_input("Enter a keyword to search:")
        
        if keyword:
            st.write(f"Search Results for '{keyword}':")
            results = search_documents(docs, keyword)
            
            if results:
                for res in results:
                    st.write(f"**Title:** {res['title']}")
                    st.write(f"**Abstract:** {res['abstract']}")
                    st.write(f"**Keywords:** {', '.join(res['keywords'])}")
                    
            else:
                st.write("No results found.")

if __name__ == "__main__":
    main()
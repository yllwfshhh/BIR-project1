from Bio import Entrez
import os
import streamlit as st
import xml.etree.ElementTree as ET

def extract_pubmed_info(xml_file_name, query):
    # Parse the XML file
    tree = ET.parse(xml_file_name)
    root = tree.getroot()

    # Initialize a dictionary to store the extracted information
    info = {}

    # Extract the article title
    article_title = root.find(".//ArticleTitle")
    print(article_title)

    if article_title is not None:
        info['ArticleTitle'] = article_title.text
    else:
        info['ArticleTitle'] = "No Title Found"

    # Add more fields to extract as needed (e.g., abstract, authors)
    # Example: Extracting the abstract
    abstract = root.find(".//AbstractText")
    if abstract is not None:
        info['Abstract'] = abstract.text
    else:
        info['Abstract'] = "No Abstract Found"

    return info

# 设置PubMed数据库的邮箱地址
def Download(query,num) :
    Entrez.email = "yllwfshhh@gmail.com"  # 请替换为您的邮箱地址
    xml_dir = "data"  # 存储XML文件的目录
    # 搜索PubMed数据库
    handle = Entrez.esearch(db="pubmed", term=query, retmax=num)  # 最多检索100篇文章
    record = Entrez.read(handle)
    print(record)
    handle.close()

    if record['IdList'] == None :
        st.markdown("## **No Result **")
    # 下载并保存XML文件
    base_path = os.path.join(xml_dir,query)

    if os.path.isdir(base_path):
        print("folder exists")
    else:
        print("forder didn't exists")
        print(base_path)
        os.makedirs(os.path.join(base_path),exist_ok=True)
        os.makedirs(os.path.join(base_path,"xml"),exist_ok=True)
        # os.makedirs(os.path.join(base_path,"pkl"),exist_ok=True)

    for pubmed_id in record["IdList"]:
        fetch_handle = Entrez.efetch(db="pubmed", id=pubmed_id, rettype="xml", retmode="xml")
        xml_data = fetch_handle.read()
        fetch_handle.close()

        # 构造XML文件的文件名，通常使用PubMed ID
        xml_file_name = os.path.join(base_path,"xml", f"X{pubmed_id}.xml")

        # 将XML数据保存到文件，使用二进制模式 'wb' 来写入字节数据
        with open(xml_file_name, "wb") as xml_file:
            xml_file.write(xml_data)

        temp = extract_pubmed_info(xml_file_name,query)
        st.markdown("Successful Download !! **PMID:** {}".format(pubmed_id))
        st.markdown("**Title :** {}".format(temp['ArticleTitle']))

    return

Download("obesity",5)
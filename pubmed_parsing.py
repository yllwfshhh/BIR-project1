from Bio import Entrez
import os

# 设置PubMed数据库的邮箱地址
def parsing(query,num) :

    Entrez.email = "yllwfshhh@gmail.com"  # 请替换为您的邮箱地址
    xml_dir = "data"  # 存储XML文件的目录
    # 搜索PubMed数据库
    handle = Entrez.esearch(db="pubmed", term=query, retmax=num)  # 最多检索100篇文章
    record = Entrez.read(handle)
    print(record)
    handle.close()

    if record['IdList'] == None :
        print("## **No Result **")
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
        
        print(f"pubmed id:{pubmed_id} downloaded sucessful")

    return

parsing("cancer",5)
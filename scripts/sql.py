import sqlite3

# Connect to your SQLite database
conn = sqlite3.connect('pubmed.db')
c = conn.cursor()

# Define the SQL query
sql_query = '''
SELECT id, pmid, title, pubdate, abstract 
FROM pubmed_articles 
WHERE (title LIKE ? OR abstract LIKE ?)
AND SUBSTR(pubdate, 1, 4) >= ?
'''
# 
# Define the parameters
parameters = ('%ente%', '%ente%','2018')

# Execute the query
c.execute(sql_query, parameters)

# Fetch all results
results = c.fetchall()
print('hello')


# Print the results
for row in results:
    print(f"ID: {row[0]}, PMID: {row[1]}, Title: {row[2]}, Pubdate: {row[3]}")

# Close the connection
conn.close()
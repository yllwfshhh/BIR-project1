import streamlit as st

# Define CSS to set background color for the main content area
page_bg_css = '''
<style>
    [data-testid="stAppViewContainer"] {
        background-color: #ADD8E6;  /* Light blue background color */
    }
</style>
'''

# Apply the custom CSS
st.markdown(page_bg_css, unsafe_allow_html=True)

# Content of your Streamlit app
st.title("Streamlit App with Custom Background Color")
st.write("The main content area has a light blue background.")

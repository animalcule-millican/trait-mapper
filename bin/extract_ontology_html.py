#!/usr/bin/env python3
from bs4 import BeautifulSoup
import pandas as pd
# Open the HTML file and parse it with BeautifulSoup
with open('/home/glbrc.org/millican/repos/trait-mapper/etc/plabse_ontology_map.html', 'r') as f:
    soup = BeautifulSoup(f, 'html.parser')

# Find all elements that have the "value", "data-section", and "data-key" attributes
elements = soup.find_all(attrs={'value': True, 'data-section': True, 'data-key': True})

# Extract the values of these attributes
html_dict = {}
for element in elements:
    value = element['value']
    data_section = element['data-section']
    data_key = element['data-key'].split(">")[-1]
    print(f'value: {value}, data-section: {data_section}, data-key: {data_key}')
    html_dict[value] = {'pgptid': value, 'onotology': data_section, 'name_function': data_key}

# Create a DataFrame from the dictionary
df = pd.DataFrame.from_dict(html_dict, orient='index')
df.to_csv('/home/glbrc.org/millican/repos/trait-mapper/etc/plabse_ontology_map.csv', index=False)
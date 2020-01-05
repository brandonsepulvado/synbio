import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

# set directory
files_to_read  = list(Path('/Users/brandonsepulvado/Documents/synbio/data/web_of_science').glob('*.xlsx'))

# create list of data frames
data_list = [pd.read_excel(f) for f in files_to_read]

# inspect
for i,df in enumerate(data_list):
    df_shape = df.shape
    if df_shape != (500, 67):
        print(df_shape, files_to_read[i])

# there are two that do not have the standard dimensions : first and last

# combine the data frames together
data = pd.concat(data_list)

# clean data ===================================================================

# which variables to keep
vars_to_keep = ['PT', 'AU', 'OA', 'PM', 'PY', 'TI', 'AB', 'DE']

# subset data based upon selected variables
data_subset = data[vars_to_keep]

# how many with missing publication year?
data_subset['PY'].isnull().sum()

# which ones
data_subset[data_subset['PY'].isnull()]

# keep only rows not missing a year
data_subset = data_subset[data_subset['PY'].isnull() == False]

# make year an integer
data_subset.loc[:,'PY'] = data_subset.loc[:,'PY'].astype(int)

# create dictionary with new names
rename_dict = {
    'PT':'pub_type', 
    'AU' : 'authors', 
    'OA' : 'open_access', 
    'PM' : 'pubmed_id', 
    'PY' : 'pub_year', 
    'TI' : 'title', 
    'AB' : 'abstract', 
    'DE' : 'author_keywords'
}

# apply dictionary to give columns substantive names
data_subset = data_subset.rename(columns = rename_dict)

# make abstract and keywords lowercase
data_subset.loc[:,'abstract'] = data_subset.loc[:,'abstract'].str.lower()
data_subset.loc[:,'author_keywords'] = data_subset.loc[:,'author_keywords'].str.lower()

# basic descriptives ===========================================================

# number of publications per year
pubs_per_year = data_subset['pub_year'].value_counts().sort_index()

# plot the value counts
pubs_per_year.plot(kind='line')
plt.xticks(rotation=60, fontsize=8)
plt.yticks(fontsize=8)
plt.show()

# how many of each publication type
data_subset['pub_type'].value_counts()

# how many abstracts contain ethic*
keywords_around_ethic = data_subset[data_subset['abstract'].str.contains('ethic', na = False)] \
    .author_keywords.dropna()

# get single keyword or keyphrase per row
keywords_around_ethic = keywords_around_ethic.str.split(pat=';').explode()

# remove whitespace
keywords_around_ethic = keywords_around_ethic.str.strip()

# get keyword counts
keywords_ethic_counts = keywords_around_ethic.value_counts()

# get pubmed identifiers for ethics-related keywords

# TODO : make a decision about which other keywords are ethics related

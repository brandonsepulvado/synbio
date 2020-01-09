# ==============================================================================
# implementing ner in spacy
# ==============================================================================

# preliminaries ================================================================

# import modules
import pandas as pd
import spacy
import textacy
from pathlib import Path

# import data 

# set data director path
data_folder = Path('./data/web_of_science')

# set file name
file_name = Path('data_subset.csv')

# import data
df = pd.read_csv(data_folder / file_name)

# create a list of abstracts 
abstract_list = [abstract for abstract in df[df['abstract'].str.len() > 0.0]['abstract']]

# out of the box spacy =========================================================

# create corpus
abstract_corpus = textacy.Corpus(lang=textacy.load_spacy_lang('eng_core_web_lg'), data=abstract_list)
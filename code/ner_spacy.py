# ==============================================================================
# implementing ner in spacy
# ==============================================================================

# preliminaries ================================================================

# import modules
import pandas as pd
import textacy
from pathlib import Path

import spacy
from spacy.lang.char_classes import ALPHA, ALPHA_LOWER, ALPHA_UPPER
from spacy.lang.char_classes import CONCAT_QUOTES, LIST_ELLIPSES, LIST_ICONS
from spacy.util import compile_infix_regex

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
abstract_corpus = textacy.Corpus(lang=textacy.load_spacy_lang('en_core_web_lg'), data=abstract_list)

for token in abstract_corpus[0]:
    print(token.text)

# TODO : add custom tokenization rule 

# default tokenizer
nlp = spacy.load("en_core_web_lg")
doc = nlp(abstract_list[257])

# modify tokenizer infix patterns
infixes = (
    LIST_ELLIPSES
    + LIST_ICONS
    + [
        r"(?<=[0-9])[+\-\*^](?=[0-9-])",
        r"(?<=[{al}{q}])\.(?=[{au}{q}])".format(
            al=ALPHA_LOWER, au=ALPHA_UPPER, q=CONCAT_QUOTES
        ),
        r"(?<=[{a}]),(?=[{a}])".format(a=ALPHA),
        # EDIT: commented out regex that splits on hyphens between letters:
        #r"(?<=[{a}])(?:{h})(?=[{a}])".format(a=ALPHA, h=HYPHENS),
        r"(?<=[{a}0-9])[:<>=/](?=[{a}])".format(a=ALPHA),
    ]
)

infix_re = compile_infix_regex(infixes)
nlp.tokenizer.infix_finditer = infix_re.finditer
doc = nlp("mother-in-law")

for token in doc:
    print(token.text)
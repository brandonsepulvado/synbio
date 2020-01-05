# this file establishes the synbio class for pubmed queries using bio.entrez

# import relevant modules
from Bio import Entrez
import pandas as pd

class SynBioPubmedQuery:
    # TODO : add a db attribute to specify pubmed
    def __init__(self, email):
        self.email = email
        self.is_ident_set = False
    
    def set_identity(self):
        Entrez.email = self.email
        self.is_ident_set = True

    def get_id_list(self, search_phrase, max_to_retain=10):
        if self.is_ident_set:
            handle = Entrez.esearch(db="pubmed", term=search_phrase, idtype="acc", retmax=max_to_retain)
            record = Entrez.read(handle)
            handle.close()
            self.pubmed_ids = record['IdList']
        else:
            print('You must first set your identity.')
    # TODO : make the count a part of the get id list 
    # so that there is no need to access server for id list and counts 
    def get_pub_count(self, search_phrase):
        if self.is_ident_set:
            handle_count = Entrez.esearch(db="pubmed", rettype='count', term=search_phrase, idtype="acc")
            record_count = Entrez.read(handle_count)
            handle_count.close()
            self.num_pubs = record_count['Count']
        else:
            print('You must first set your identity.')
    
    def get_title_abstract(self):
        handle_title_abstract = Entrez.efetch(db="pubmed", id=self.pubmed_ids, rettype="xml", retmode="text")
        records_title_abstract = Entrez.read(handle_title_abstract)
        handle_title_abstract.close()
        test_list = []
        for item in records_title_abstract['PubmedArticle']:
            dict_title_abstract = {
                'title':item['MedlineCitation']['Article']['ArticleTitle'],
                'abstract':item['MedlineCitation']['Article']['Abstract']['AbstractText']
                }
            test_list.append(dict_title_abstract)
        # save this list as a list, new attribute    
        self.title_and_abstract = test_list

    def convert_to_df(self):
        self.data = pd.DataFrame.from_dict(self.title_and_abstract)


        

test_class = SynBioPubmedQuery('sepulvado-brandon@norc.org')
print(test_class.email)
print(test_class.is_ident_set)
test_class.set_identity()
print(test_class.is_ident_set)
test_class.get_pub_count(search_phrase="synthetic biology[TIAB]")
print(test_class.num_pubs)
test_class.get_id_list(search_phrase="synthetic biology[TIAB]")
print(test_class.pubmed_ids)
test_class.get_title_abstract()
print(test_class.title_and_abstract)
test_class.convert_to_df()
print(test_class.data)
# TODO : extract abstract from list
print(test_class.data['abstract'][0])
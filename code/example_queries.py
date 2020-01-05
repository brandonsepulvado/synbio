from Bio import Entrez

Entrez.email = 'sepulvado-brandon@norc.org'
# handle = Entrez.esummary(db="pubmed", id="19304878,14630660", retmode="xml")
# records = Entrez.parse(handle)

# for record in records:
#     # each record is a Python dictionary or list.
#     print(record['Title'])
# handle.close()

# handle = Entrez.esearch(db="pubmed", term="'synthetic biological'[TIAB]", idtype="acc", retmax=10)
handle = Entrez.esearch(db="pubmed", rettype='count', term="(synthetic biology[TIAB]) AND ethic*[TIAB]", idtype="acc")
record = Entrez.read(handle)
handle.close()
print(record['Count'])

# record = Entrez.read(handle)
# handle.close()

# write function to get ids; must have previously setup Entrez info
def get_id_list(search_phrase, max_to_retain=10):
    handle = Entrez.esearch(db="pubmed", term=search_phrase, idtype="acc", retmax=max_to_retain)
    record = Entrez.read(handle)
    handle.close()
    return record['IdList']

# write a function that returns the title and abstract for a given pubmed id
def get_title_abstract(pubmed_id):
    handle = Entrez.efetch(db="pubmed", id=pubmed_id, rettype="xml", retmode="text")
    records = Entrez.read(handle)
    handle.close()
    #abstract = records['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
    dict_title_abstract = {
        'title':records['PubmedArticle'][0]['MedlineCitation']['Article']['ArticleTitle'],
        'abstract':records['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
    }
    return dict_title_abstract

# # redo for multiple pubmed ids
# def get_title_abstract(pubmed_ids):
#     handle = Entrez.efetch(db="pubmed", id=pubmed_ids, rettype="xml", retmode="text")

# test on synthetic biology[TIAB]
test_result = get_id_list("synthetic biolog*[TIAB]")
print(test_result)
first_id = test_result[0]

# test abstract function
print(get_title_abstract(first_id))
print(get_title_abstract(test_result[0:2]))

handle = Entrez.efetch(db="pubmed", id=test_result[0:5], rettype="xml", retmode="text")
records = Entrez.read(handle)
handle.close()
# dict_title_abstract = {
#     'title':records['PubmedArticle'][0]['MedlineCitation']['Article']['ArticleTitle'],
#     'abstract':records['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
# print(records['PubmedArticle'][0]['MedlineCitation']['Article']['ArticleTitle'])

# for i in enumerate(records['PubmedArticle']):
#     records['PubmedArticle'][i]['MedlineCitation']['Article']['ArticleTitle']
#     dict_title_abstract = {
#         'title':records['PubmedArticle'][i]['MedlineCitation']['Article']['ArticleTitle']
#         'abstract':records['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
#     }

test_list = []
for item in records['PubmedArticle']:
    # print(item['MedlineCitation']['Article']['ArticleTitle'])
    # print(item['MedlineCitation']['Article']['Abstract']['AbstractText'])
    dict_title_abstract = {
        'title':item['MedlineCitation']['Article']['ArticleTitle'],
        'abstract':item['MedlineCitation']['Article']['Abstract']['AbstractText']
    }
    print(dict_title_abstract)
    test_list.append(dict_title_abstract)



# # quick workflow
# handle = Entrez.esearch(db="pubmed", term="synthetic biology[TIAB] NOT 'synthetic biology'[TIAB]", idtype="acc", retmax=10)
# record = Entrez.read(handle)
# handle.close()
# print(record)



# print(Entrez.esummary(db="pubmed", id=record['IdList'][0]))

# handle1 = Entrez.efetch(db="pubmed", id=first_id)
# record1 = Entrez.read(handle1)
# info = record1[0]["TitleMainList"][0]

# handle = Entrez.efetch(db="pubmed", id=first_id,
#                        rettype="xml", retmode="text")
# records = Entrez.read(handle)
# abstracts = [pubmed_article['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
#              for pubmed_article in records['PubmedArticle']]

# # TODO : unit test to ensure that the index will always be one
# print(records['PubmedArticle'][0])
# # TODO : unit test to ensure that final index will always be zero (always only one abstract)
# print(records['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText'][0])
# print(records['PubmedArticle'][0]['MedlineCitation']['Article']['ArticleTitle'])



# print([pubmed_article['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
#              for pubmed_article in records['PubmedArticle']])
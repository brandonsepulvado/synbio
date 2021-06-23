import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms import bipartite
from networkx.algorithms import community
import matplotlib.pyplot as plt
from networkx import edge_betweenness_centrality as betweenness
import collections

# set directory
files_to_read  = list(Path('/Users/brandonsepulvado/Documents/synbio/data/web_of_science').glob('*.xlsx'))

# create list of data frames
data_list = [pd.read_excel(f, engine='openpyxl') for f in files_to_read]

# inspect
for i,df in enumerate(data_list):
    df_shape = df.shape
    if df_shape != (500, 67):
        print(df_shape, files_to_read[i])

# there are two that do not have the standard dimensions : first and last

# combine the data frames together
data_pubs = pd.concat(data_list)

# clean data ===================================================================

# which variables to keep
vars_to_keep = ['PT', 'AU', 'OA', 'PM', 'PY', 'TI', 'AB', 'DE', 'SO', 'WC', 'SC']

# subset data based upon selected variables
data_subset = data_pubs[vars_to_keep] 

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
    'DE' : 'author_keywords',
    'SO' : 'pub_title',
    'WC' : 'wos_categs',
    'SC' : 'research_areas'
}

# apply dictionary to give columns substantive names
data_subset = data_subset.rename(columns = rename_dict)

# make abstract and keywords lowercase
data_subset.loc[:,'abstract'] = data_subset.loc[:,'abstract'].str.lower()
data_subset.loc[:,'author_keywords'] = data_subset.loc[:,'author_keywords'].str.lower()
data_subset.loc[:,'pub_title'] = data_subset.loc[:,'pub_title'].str.lower()


# save data object to csv
data_subset.to_csv('/Users/brandonsepulvado/Documents/synbio/data/web_of_science/data_subset.csv', index=False)

# save ethics data =============================================================

# exclude rows with missing years
data_ethics = data_pubs[data_pubs['PY'].isnull() == False]

# exclude cases with missing abstracts
data_ethics = data_pubs[data_pubs['AB'].isnull() == False]

# keep only those cases that pertain to ethics
data_ethics = data_ethics[data_ethics['AB'].str.contains("ethic|security|safe|dilemma")]

# make all author names lowercase
data_ethics['AU'] = data_ethics['AU'].str.lower()

# save data object to csv
# data_ethics.to_csv('/Users/brandonsepulvado/Documents/synbio/data/web_of_science/data_ethics.txt', sep = '\t', index=False)

# create ethics network ========================================================

# read data
data_ethics = pd.read_excel('/Users/brandonsepulvado/Documents/synbio/data/web_of_science/data_processed/data_ethics.xlsx')

# make all author names lowercase
data_ethics['AU'] = data_ethics['AU'].str.lower()
data_ethics['AB'] = data_ethics['AB'].str.lower()

net_ethics = data_ethics[['AU', 'UT']].reset_index(drop = True)

net_ethics['AU'] = net_ethics['AU'].str.split(pat = ';')

net_ethics = net_ethics.explode('AU').reset_index(drop = True)

# remove wos:
net_ethics['UT'] = net_ethics['UT'].str.replace('WOS:', '')

# create author identifier
net_ethics['auth_id'] = range(0, len(net_ethics))


# create edge tuple column
net_ethics['edge_tuple'] = list(zip(net_ethics['auth_id'], net_ethics['UT']))

vertices_author = pd.DataFrame(net_ethics['auth_id'].unique(), columns = ['node'])

# provide type for bipartite graph
vertices_author['type'] = False
# create publication vertices
vertices_pub = pd.DataFrame(net_ethics['UT'].unique(), columns = ['node'])
# provide type for bipartite graph
vertices_pub['type'] = True
# combined
vertices = pd.concat([vertices_author[['node', 'type']], vertices_pub]).reset_index(drop = True)

# pd.merge(net_ethics, vertices_author, 
#                      left_on = 'AU', 
#                      right_on = 'node', 
#                      how='left')

# g = igraph.Graph.Bipartite(types=vertices['type'], edges=net_ethics['edge_tuple'], directed=False)

B = nx.Graph()
# Add nodes with the node attribute "bipartite"
B.add_nodes_from(net_ethics['AU'], bipartite=0)
B.add_nodes_from(net_ethics['UT'], bipartite=1)
# Add edges only between nodes of opposite node sets
B.add_edges_from(list(zip(net_ethics['AU'], net_ethics['UT'])))
nx.is_connected(B)
top_nodes = {n for n, d in B.nodes(data=True) if d['bipartite']==0}
bottom_nodes = set(B) - top_nodes
bipartite.degree_centrality(B, top_nodes)
B_onemode = bipartite.collaboration_weighted_projected_graph(B, top_nodes)
B_onemode.nodes()
B_onemode.number_of_edges()
B_onemode.number_of_nodes()
nx.number_of_isolates(B_onemode)
B_onemode.degree(weight='weight')

# save as tsv
nx.write_weighted_edgelist(B_onemode, path = Path('./output/net_collab.txt'), delimiter='\t')
B_onemode = nx.read_weighted_edgelist(path = Path('./output/net_collab.txt'), delimiter='\t')


def most_central_edge(G):
    centrality = betweenness(G, weight='weight')
    return max(centrality, key=centrality.get)

comp = community.girvan_newman(B_onemode, most_valuable_edge=most_central_edge)

nx.draw(B_onemode)
plt.show()
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

# overall network (not just ethics) ============================================

# save data pubs
data_pubs.to_csv(Path('./output/data_pubs_wos.txt'), sep='\t')

# using data_pubs

# keep only rows not missing a year
data_pubs = data_pubs[data_pubs['PY'].isnull() == False].copy()

# make year an integer
data_pubs.loc[:,'PY'] = data_pubs.loc[:,'PY'].astype(int)

# make author names lowercase
data_pubs['AU'] = data_pubs['AU'].str.lower()

# keep only author and pub id (creating bipartite edgelist)
net_full = data_pubs[['AU', 'UT']].reset_index(drop = True).copy()

# get individual authors per paper
net_full['AU'] = net_full['AU'].str.split(pat = ';')

# make one author-paper per row
net_full = net_full.explode('AU').reset_index(drop = True)

# ensure no whitespace from author unnesting
net_full['AU'] = net_full['AU'].str.strip()

# remove wos:
# net_ethics['UT'] = net_ethics['UT'].str.replace('WOS:', '')

# create author identifier
net_full['auth_id'] = range(0, len(net_full))

# create edge tuple column
net_full['edge_tuple'] = list(zip(net_full['auth_id'], net_full['UT']))

vertices_author_full = pd.DataFrame(net_full['auth_id'].unique(), columns = ['node'])

# provide type for bipartite graph
vertices_author_full['type'] = False

# create publication vertices
vertices_pub_full = pd.DataFrame(net_full['UT'].unique(), columns = ['node'])

# provide type for bipartite graph
vertices_pub_full['type'] = True

# combined
vertices_full = pd.concat([vertices_author_full[['node', 'type']], vertices_pub_full]).reset_index(drop = True)


# g = igraph.Graph.Bipartite(types=vertices['type'], edges=net_ethics['edge_tuple'], directed=False)

g_full = nx.Graph()
# Add nodes with the node attribute "bipartite"
g_full.add_nodes_from(net_full['AU'], bipartite=0)
g_full.add_nodes_from(net_full['UT'], bipartite=1)
# Add edges only between nodes of opposite node sets
g_full.add_edges_from(list(zip(net_full['AU'], net_full['UT'])))
nx.is_connected(g_full)
top_nodes = {n for n, d in g_full.nodes(data=True) if d['bipartite']==0}
bottom_nodes = set(g_full) - top_nodes
bipartite.degree_centrality(g_full, top_nodes)
g_full_onemode = bipartite.collaboration_weighted_projected_graph(g_full, top_nodes)
g_full_onemode.nodes()
g_full_onemode.number_of_edges()
g_full_onemode.number_of_nodes()
nx.number_of_isolates(g_full_onemode)
g_full_onemode.degree(weight='weight')

# get degree distribution

synbio_deg_seq_full = sorted([d for n, d in g_full_onemode.degree(weight='weight')], reverse=True)  # degree sequence
synbio_deg_count_full = collections.Counter(synbio_deg_seq_full)
deg_full, cnt_full = zip(*synbio_deg_count_full.items())

# histogram: not helpful
# fig, ax = plt.subplots()
# plt.bar(deg_full, cnt_full, width=0.80, color="b")

# plt.title("Degree Histogram")
# plt.ylabel("Count")
# plt.xlabel("Degree")
# ax.set_xticks([d + 0.4 for d in deg_full])
# ax.set_xticklabels(deg_full)
# plt.show()

# log log plot
m=3
deg_freq_full = nx.degree_histogram(g_full_onemode)
degrees_full = range(len(deg_freq_full))
plt.figure(figsize=(12, 8)) 
plt.loglog(degrees_full[m:], deg_freq_full[m:],'C1o') 
plt.xlabel('Degree')
plt.ylabel('Frequency')
plt.show()

# save as tsv
nx.write_weighted_edgelist(B_onemode, path = Path('./output/net_collab_full.txt'), delimiter='\t')
B_onemode = nx.read_weighted_edgelist(path = Path('./output/net_collab.txt'), delimiter='\t')

nx.to_pandas_adjacency(g_full_onemode)


def most_central_edge(G):
    centrality = betweenness(G, weight='weight')
    return max(centrality, key=centrality.get)

comp = community.girvan_newman(B_onemode, most_valuable_edge=most_central_edge)

nx.draw(B_onemode)
plt.show()

# size of connected components
con_comps_indiv = [len(c) for c in sorted(nx.connected_components(B_onemode), key=len, reverse=True)]

B_onemode = bipartite.weighted_projected_graph(g_full, top_nodes)
edge_weights_synbio_indiv = nx.get_edge_attributes(B_onemode, 'weight').values()
pd.Series(edge_weights_synbio_indiv).value_counts()
plt.show()

# general synbio institution co-authorship network =============================

# load data
net_synbio_inst = nx.read_weighted_edgelist(path = Path('./output/synbio_general_inst_coauth_network.txt'), delimiter='\t')

# get degree distribution
synbio_deg_seq_full = sorted([d for n, d in net_synbio_inst.degree(weight='weight')], reverse=True)  # degree sequence
synbio_deg_count_full = collections.Counter(synbio_deg_seq_full)
synbio_deg_full, synbio_cnt_full = zip(*synbio_deg_count_full.items())

# log log plot
m=3
synbio_deg_full = nx.degree_histogram(net_synbio_inst)
synbio_degrees_full = range(len(synbio_deg_full))
plt.figure(figsize=(12, 8)) 
plt.loglog(synbio_degrees_full[m:], synbio_deg_full[m:],'C1o') 
plt.xlabel('Degree')
plt.ylabel('Frequency')
plt.show()

net_synbio_inst.number_of_edges()
net_synbio_inst.number_of_nodes()

con_comps_inst = [len(c) for c in sorted(nx.connected_components(net_synbio_inst), key=len, reverse=True)]

pd.Series(net_synbio_inst).hist()
plt.show()

edge_weights_synbio_inst = nx.get_edge_attributes(net_synbio_inst, 'weight').values()
pd.Series(edge_weights_synbio_inst).hist()
plt.show()
import numpy as np
import pandas as pd
import networkx as nx
from networkx.readwrite import json_graph
import json
import run_server

# data_path = "../../data/"
# sub_dir = "reactome/67/"
# relation_file = data_path + sub_dir + "ReactomePathwaysRelation.txt"
# url_to_relation_file = "https://reactome.org/download/current/ReactomePathwaysRelation.txt"

def generate_root_node(relation_file, root_node_id = 'Homo_Sapiens', organism = "HSA"):

    rel_df = pd.read_csv(relation_file, sep = "\t", header = None)

    rel_df.columns = ['parentId', 'id']
    cut = rel_df['parentId'].str.contains(organism) & rel_df['id'].str.contains(organism)
    rel_df = rel_df.loc[cut]

    G = nx.DiGraph()
    G.add_edges_from(rel_df.values)
    roots = [n for n,d in G.in_degree() if d==0]

    roots_df = pd.DataFrame(columns = [['parentId', 'id']])
    roots_df['id'] = roots
    roots_df['parentId'] = root_node_id

    roots_df = pd.DataFrame(roots_df.values, columns = ['parentId', 'id'])
    rel_df = pd.DataFrame(rel_df.values, columns = ['parentId', 'id'])

    tree = roots_df.append(rel_df)
    return tree

def generate_sunburst_json(stats_df, relation_file, root_node_id = 'Homo_Sapiens'):
    '''
    stats_df has a very specific format:
    indexes  should be the pathway IDs
    'value' will be translated to colors in the plot.
    'ngenes' is the size of the pathwat, and will be displayed as the size of each arc.
    'descr' as the description of the pathway

    relation_file is downloaded from Reactome "ReactomePathwaysRelation.txt"
    '''

    # in_df = in_df.loc[[x for x in in_df.index if 'HSA' in x]]

    rel_df = generate_root_node(relation_file, root_node_id = root_node_id)
    topPaths = rel_df.loc[(rel_df['parentId'] == root_node_id), 'id']
    rootNgenes = np.sum(in_df.loc[[x in topPaths.tolist() for x in in_df.index],'ngenes'])
    rootNode = pd.DataFrame([[1,rootNgenes,"Homo Sapiens", 0,0,0,0]], columns = ["value", "ngenes", "descr"]).xs(0)
    rootNode.name = root_node_id

    stats_df = stats_df.append(rootNode)
    topDict = stats_df.to_dict()

    pathways = stats_df.index
    n_path = len(pathways)

    subset_vec = [x in pathways for x in rel_df.iloc[:,0]] and [x in pathways for x in rel_df.iloc[:,1]]
    sub_rel_df = rel_df[subset_vec]

    G = nx.DiGraph()

    G.add_nodes_from(pathways)
    G.add_edges_from(sub_rel_df.values)

    tree = nx.algorithms.dag.dag_to_branching(G)

    secondDict = nx.get_node_attributes(tree,'source')

    thirdDict = {'value':{}, 'ngenes':{}, 'descr':{}}
    for key, value in secondDict.items():
        thirdDict['value'].update({key : topDict['value'][value]})
        thirdDict['ngenes'].update({key : topDict['ngenes'][value]})
        thirdDict['descr'].update({key : topDict['descr'][value]})


    nx.set_node_attributes(tree, thirdDict['value'], name = 'value')
    nx.set_node_attributes(tree, thirdDict['ngenes'], name = 'ngenes')
    nx.set_node_attributes(tree, thirdDict['descr'], name = 'descr')

    root = [v for v, d in tree.in_degree() if d == 0][0]
    out_json = json_graph.tree_data(tree, root)

    return out_json

def prepare_tree_df(in_df, rectome_df):
    ngenes = reactome_df.groupby('reactome_id').count()['gene']
    descr = rectome_df.groupby('reactome_id').first()['reactome_name']
    in_df['ngenes'] = ngenes
    in_df['descr'] = descr
    in_df.columns = ['value','ngenes','descr']
    return in_df

def default(o):
     if isinstance(o, np.integer): return int(o)
     raise TypeError

def write_json(json, file_name):
    with open(file_name, 'w') as outfile:
        json.dump(json, outfile, default=default)

def generate_sunburst(values_df, reactome_df, relation_file, filename, root_node_id = 'Homo_Sapiens'):
    # relation_file is downloaded from Reactome "ReactomePathwaysRelation.txt"
    stats_df = prepare_tree_df(values_df, reactome_df)
    json = generate_sunburst_json(stats_df, relation_file, root_node_id)
    write_json(json, filename)
    run_server.run_sunburst(path='.')

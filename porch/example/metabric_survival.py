import numpy as np
import pandas as pd
import porch
from lifelines import CoxPHFitter


def load_metabric(path, relevant_columns = ['METABRIC_ID', 'age_at_diagnosis', 'last_follow_up_status','Treatment','T', 'Pam50Subtype', 'ER.Expr', 'Her2.Expr', 'PR.Expr']):
     
    discovery_ex_path = path + '/discovery_ExpressionMatrix.txt'
    validation_ex_path = path + '/validation_ExpressionMatrix.txt'
    normals_ex_path = path + '/normals_ExpressionMatrix.txt'
    
    discovery_cl_path = path + '/table_S2_revised.txt'
    validation_cl_path = path + '/table_S3_revised.txt'
    
    discovery_ex = pd.read_csv(discovery_ex_path, sep = '\t')
    validation_ex = pd.read_csv(validation_ex_path, sep = ' ')
    normals_ex = pd.read_csv(normals_ex_path, sep = '\t')
    
    discovery_cl = pd.read_csv(discovery_cl_path, sep = '\t')
    validation_cl = pd.read_csv(validation_cl_path, sep = '\t')
    
    # relevant_columns = ['METABRIC_ID', 'age_at_diagnosis', 'last_follow_up_status','Treatment','T', 'Pam50Subtype', 'ER.Expr', 'Her2.Expr', 'PR.Expr']
    
    discovery_cl = discovery_cl[relevant_columns].set_index('METABRIC_ID')
    validation_cl = validation_cl[relevant_columns].set_index('METABRIC_ID')
    
    expression = discovery_ex.join(validation_ex)
    clinical = discovery_cl.append(validation_cl)
    
    data = clinical.join(expression.T).T
     
    return data

def illumina2ensembl_dictionary(path):
    """
    This function needs a file with Illumina probe IDs and Ensembl IDs downloaded from Biomart
    """
    dictionary_pd = pd.read_csv(path, sep = '\t').dropna()
    dictionary_pd.columns = ['gene','probe']
    # dictionary = dictionary_pd.groupby(['gene','probe']).size()
    return dictionary_pd

def porch_metabric(expression_df, illumina2ensembl_path = '../../../data/reactome/illumina2ensembl.txt'):
    reactome_df = porch.get_reactome_df(organism = "HSA", gene_anot = "Ensembl")
    illumina2ensembl = illumina2ensembl_dictionary(illumina2ensembl_path)
    illumina2ensembl = illumina2ensembl.set_index('probe').to_dict()['gene']
    expression_df.index = expression_df.index.to_series().map(lambda x: illumina2ensembl.get(x,'noEnsembl'))
    activity, pathways = porch.porch_single_process(expression_df, reactome_df, gene_column = 'gene', set_column = 'reactome_id')
    return activity

def metabric_survival(activity_df, metadata_df):
    metadata_df.loc['E'] = (metadata_df.loc['last_follow_up_status'] == 'd-d.s.')*1
    result = pd.DataFrame()
    for index, row in activity_df.iterrows():
        print("Cox's regression: " + str(index))
        rowresult = porch.survival(row, metadata_df, duration_col = 'T', event_col = 'E')
        result = result.append(rowresult)
    return result




if __name__ == "__main__":
    metabric_path = '../../../../data/metabric'
    illumina2ensembl_path = '../../../../data/reactome/illumina2ensembl.txt'
    
    data = load_metabric(metabric_path )
    expression_df = data.iloc[8:,:]
    metadata_df = data.iloc[:8,:]

    activity = porch_metabric(expression_df,  illumina2ensembl_path)
    print('Writing metabric_activities.csv')
    activity.to_csv('metabric_activities.csv')
    results = metabric_survival(activity, metadata_df)
    print('Writing metabric_results.csv')
    results.to_csv('metabric_results.csv')



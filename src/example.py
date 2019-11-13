import numpy as np
import pandas as pd
import tarfile
import requests
import os
import porch
from biothings_client import get_client

def track_dl(url,tar):
    response = requests.get(url, stream=True)
    with open(tar, "wb") as handle:
        for data in response.iter_content():
            handle.write(data)

def get_tar(url,path):
    try:
        tf = tarfile.open(path)
    except:
        track_dl(url, path)
        tf = tarfile.open(path)
    return tf

def get_expression_data(path,url,file):
    df = get_data(path,url,file)
    df.dropna(axis=0, how='any', inplace=True)
#    df.set_index('Hugo_Symbol', inplace=True)
    df.set_index('Entrez_Gene_Id', inplace=True)
    #df.drop(columns=['Unnamed: 0', 'Entrez_Gene_Id'], inplace=True)
    #df.drop(columns=['Entrez_Gene_Id'], inplace=True)
    df.drop(columns=['Hugo_Symbol'], inplace=True)
    df = df.reindex(sorted(df.columns), axis=1)
    return df

def get_clinical_data(path,url,file):
    df = get_data(path,url,file).T
    df.columns = df.iloc[2]
    df.drop(columns=["A unique sample identifier.","STRING","1","SAMPLE_ID",'TCGA-BH-A1ES-01'], inplace=True)
    if 'TCGA-BH-A1ES-01' in df.columns:
        df.drop(columns=['TCGA-BH-A1ES-01'], inplace=True)
    df.drop(index=["Unnamed: 0","#Patient Identifier","Sample Identifier","Other Sample ID"], inplace=True)
    df = df.reindex(sorted(df.columns), axis=1)
    return df

def get_data(path,url,file):
    try:
        df = pd.read_csv(path, sep="\t")
    except:
        tf = get_tar(url,".porch/my.tar.gz")
        tf.extract(file,".porch/")
        df = pd.read_csv(".porch/" + file, sep="\t")
        df.to_csv(path, sep="\t")
    return df

def get_ensembl_ids(entrez_ids):
    mg = get_client('gene')
    ensembl_id_raw = mg.querymany(entrez_ids, scopes='entrezgene', fields='ensembl.gene', species='human')
    outlist,drop_list = [],[]
    for ret in ensembl_id_raw[1:]:
        if "ensembl" in ret:
            ret = ret['ensembl']
            if isinstance(ret, list):
                ret = ret[0]
            outlist.append(ret['gene'])
        else:
            drop_list.append(int(ret['query']))
    return outlist, drop_list

def tcga_example():
    print("Downloading data ...")
    try:
        os.mkdir(".porch")
    except FileExistsError:
        pass
    brca = get_expression_data(".porch/brca.tsv.gz", 'http://download.cbioportal.org/brca_tcga_pub2015.tar.gz',"data_RNA_Seq_v2_expression_median.txt")
    brca_clin_raw = get_clinical_data(".porch/brca_clin.tsv.gz", 'http://download.cbioportal.org/brca_tcga_pub2015.tar.gz',"data_clinical_sample.txt")
    print("Preprocessing data ...")
    brca.dropna(axis=0, how='any', inplace=True)
    brca = brca.loc[~(brca<=0.0).any(axis=1)]
    print(brca.index)
    ensembl_ids,drop_list = get_ensembl_ids(list(brca.index))
    # Convert all genenames we can, drop the rest
    #brca.rename(index=ensembl)
    #brca = brca.loc[~df.index.str.startswith('ENSG')]
    brca.drop(index=drop_list)
    brca = pd.DataFrame(data=np.log2(brca), index=ensembl_ids, columns=brca.columns)
    brca_clin = brca_clin_raw.loc[["PR status by ihc","ER Status By IHC","IHC-HER2"]].rename(index={"PR status by ihc" : "PR","ER Status By IHC":"ER","IHC-HER2":"HER2"})

    print("Run Porch ...")
    results_df,evaluation_df = porch.porch_reactome(brca,brca_clin,"HSA",["Pathway ~ C(PR)","Pathway ~ C(ER)","Pathway ~ C(HER2)"])

def main():
    tcga_example()

if __name__ == "__main__":
    main()
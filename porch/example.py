import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import tarfile
import requests
import os
from biothings_client import get_client
import porch
import qvalue as qv

cashe_directory = ".porch"
expression_name = "brca"
phenotype_name = "brca_clin"
preprcoc_prefix = "proc_"
significance_name  = "significance"
activity_name = "activity"

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
    df.columns = df.loc["Sample Identifier"]
    df.drop(columns=["A unique sample identifier.","STRING","1","SAMPLE_ID"], inplace=True,errors='ignore')
    if 'TCGA-BH-A1ES-01' in df.columns:
        df.drop(columns=['TCGA-BH-A1ES-01'], inplace=True)
    df.drop(index=["Unnamed: 0","#Patient Identifier","Sample Identifier","Other Sample ID"], inplace=True,errors='ignore')
    df = df.reindex(sorted(df.columns), axis=1)
    return df

def get_data(path,url,file):
    try:
        df = pd.read_csv(path, sep="\t",index_col=0)
    except:
        tf = get_tar(url,".porch/my.tar.gz")
        tf.extract(file,".porch/")
        df = pd.read_csv(".porch/" + file, sep="\t")
        df.to_csv(path, sep="\t")
    return df

def get_ensembl_ids(entrez_ids):
    mg = get_client('gene')
    ensembl_id_raw = mg.querymany(entrez_ids, scopes='entrezgene', fields='ensembl.gene', species='human')
    translate, drop_list = {},[]
    for ret in ensembl_id_raw:
        query = int(ret['query'])
        if "ensembl" in ret:
            ret = ret['ensembl']
            if isinstance(ret, list):
                ret = ret[0]
            translate[query] = ret['gene']
        else:
            drop_list.append(query)
    return translate, drop_list

def tcga_preprocess(brca, brca_clin):
    brca.dropna(axis=0, how='any', inplace=True)
    brca = brca.loc[~(brca<=0.0).any(axis=1)].copy()
    entrez2ensembl, drop_list = get_ensembl_ids(list(brca.index))
    # Convert all genenames we can, drop the rest
    brca.drop(index=drop_list, inplace=True)
    brca.rename(index=entrez2ensembl, inplace=True)
    #brca.index =  brca.index.map(str)
    #brca = pd.DataFrame(data=np.log2(brca.values), index=brca.index, columns=brca.columns)
    brca = pd.DataFrame(data=np.log2(brca.values), index=brca.index, columns=brca.columns)
    brca_clin = brca_clin.loc[["PR status by ihc","ER Status By IHC","IHC-HER2"]].rename(index={"PR status by ihc" : "PR","ER Status By IHC":"ER","IHC-HER2":"HER2"})
    relevant_col = brca_clin.T["PR"].isin(['Negative', 'Positive']) & brca_clin.T["ER"].isin(['Negative', 'Positive']) & brca_clin.T["HER2"].isin(['Negative', 'Positive'])
    brca_clin = brca_clin.T[relevant_col].T
    brca = brca.T[relevant_col].T
    # mapping = {'Negative': 0, 'Positive': 1, '[Not Available]': 1, 'Equivocal' : 1, 'Indeterminate' : 1}
    mapping = {'Negative': 0, 'Positive': 1}
    brca_clin.replace(mapping, inplace=True)
    return brca, brca_clin

def tcga_example():
    print("Downloading data ...")
    try:
        os.mkdir(cashe_directory)
    except FileExistsError:
        pass
    brca_path = cashe_directory + "/" + expression_name + ".tsv.gz"
    clin_brca_path = cashe_directory + "/" + phenotype_name + ".tsv.gz"
    processed_brca_path = cashe_directory + "/" + preprcoc_prefix + expression_name + ".tsv.gz"
    processed_clin_brca_path = cashe_directory + "/" + preprcoc_prefix + phenotype_name + ".tsv.gz"
    significance_path  = cashe_directory + "/" + significance_name + ".tsv"
    activity_path = cashe_directory + "/" + activity_name + ".tsv"
    if os.path.isfile(processed_brca_path) and os.path.isfile(processed_clin_brca_path):
        brca = pd.read_csv(processed_brca_path, sep="\t",index_col=0)
        brca_clin = pd.read_csv(processed_clin_brca_path, sep="\t",index_col=0)
    else:
        brca = get_expression_data(brca_path, 'http://download.cbioportal.org/brca_tcga_pub2015.tar.gz',"data_RNA_Seq_v2_expression_median.txt")
        brca_clin = get_clinical_data(clin_brca_path, 'http://download.cbioportal.org/brca_tcga_pub2015.tar.gz',"data_clinical_sample.txt")
        print("Preprocessing data ...")
        brca, brca_clin = tcga_preprocess(brca, brca_clin)
        brca.to_csv(processed_brca_path, sep="\t")
        brca_clin.to_csv(processed_clin_brca_path, sep="\t")
    if os.path.isfile(significance_path) and os.path.isfile(activity_path):
        significance = pd.read_csv(significance_path, sep="\t",index_col=["Test","Variable"]).T
        activity = pd.read_csv(activity_path, sep="\t",index_col=0)
    else:
        print("Run Porch ...")
        significance,activity,untested = porch.porch_reactome(brca,brca_clin,"HSA",["Pathway ~ C(PR) + C(HER2)","Pathway ~ C(PR) + C(ER) + C(PR):C(ER)","Pathway ~ C(PR) + C(ER)"])
        significance.T.to_csv(significance_path, sep="\t",index_label=["Test","Variable"])
        activity.to_csv(activity_path, sep="\t")
    ## Plot the activity of R-HSA-8931987 in TNBC vs non-TNBC
    tripple_neg = (brca_clin.T["PR"] == 0) & (brca_clin.T["ER"] == 0) & (brca_clin.T["HER2"] == 0)
    runx1 = activity.loc["R-HSA-8931987",:].T
    fig = plt.figure(figsize=(10,6))
    sns.distplot(runx1[tripple_neg], kde=False,norm_hist=True)
    sns.distplot(runx1[~tripple_neg], kde=False,norm_hist=True)
    plt.legend(labels=['TNBC','Non-TNBC'],loc='upper right')
    plt.ylabel('Density')
    plt.xlabel('Activity of \"RUNX1 regulates estrogen receptor mediated transcription\"')
#    plt.show()
    plt.savefig("rux1.png")

    # Scatterplots
    qv.qvalues(significance, ("Pathway ~ C(PR) + C(ER)","C(PR)"), ("Pathway ~ C(PR) + C(ER)","q value PR"))
    qv.qvalues(significance, ("Pathway ~ C(PR) + C(ER)","C(ER)"), ("Pathway ~ C(PR) + C(ER)","q value ER"))
    fig = plt.figure(figsize=(10,6))
    g = sns.scatterplot(data=significance,x=("Pathway ~ C(PR) + C(ER)","q value PR"),y=("Pathway ~ C(PR) + C(ER)","q value ER"))
    g.set_xscale('log')
    g.set_yscale('log')
    g.set_xlim(10e-80,0.5)
    g.set_ylim(10e-50,0.5)
    plt.savefig("PRvsER.png")

    qv.qvalues(significance, ("Pathway ~ C(PR) + C(ER) + C(PR):C(ER)","C(PR)"), ("Pathway ~ C(PR) + C(ER) + C(PR):C(ER)","q value PR"))
    qv.qvalues(significance, ("Pathway ~ C(PR) + C(ER) + C(PR):C(ER)","C(PR):C(ER)"), ("Pathway ~ C(PR) + C(ER) + C(PR):C(ER)","q value PR:ER"))
    fig = plt.figure(figsize=(10,6))
    g = sns.scatterplot(data=significance,x=("Pathway ~ C(PR) + C(ER) + C(PR):C(ER)","q value PR"),y=("Pathway ~ C(PR) + C(ER) + C(PR):C(ER)","q value PR:ER"))
    g.set_xscale('log')
    g.set_yscale('log')
   g.set_xlim(10e-80,0.5)
    g.set_ylim(10e-5,0.9)
    g.autoscale()
    plt.savefig("PRvsPRER.png")
    plt.show()



def main():
    tcga_example()

if __name__ == "__main__":
    main()

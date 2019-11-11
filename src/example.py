import pandas as pd
import tarfile
import requests
import os
import porch

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
    df.set_index('Hugo_Symbol', inplace=True)
    df.drop(columns=['Unnamed: 0', 'Entrez_Gene_Id'], inplace=True)
    #df.drop(columns=['Entrez_Gene_Id'], inplace=True)
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
        tf.extract(file)
        df = pd.read_csv(file, sep="\t")
        df.to_csv(path, sep="\t")
    return df

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
    brca = pd.DataFrame(data=np.log2(brca),index=brca.index,columns=brca.columns)
    brca_clin = brca_clin_raw[["PR status by ihc","ER Status By IHC","IHC-HER2"],:].rename(index={"PR status by ihc" : "PR","ER Status By IHC":"ER","IHC-HER2":"HER2"})

    print("Run Porch ...")
    results_df,evaluation_df = porch_reactome(brca,brca_clin,["Pathway ~ C(PR)","Pathway ~ C(ER)","Pathway ~ C(HER2)"])

def main():
    tcga_example()

if __name__ == "__main__":
    main()

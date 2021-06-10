import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import tarfile
import urllib.request
import os

from biothings_client import get_client
import porch
import porch.qvalue as qv

cache_directory = ".porch"
expression_name = "brca"
phenotype_name = "brca_clin"
preprcoc_prefix = "proc_"
significance_name  = "significance"
activity_name = "activity"

tcga_url = 'https://cbioportal-datahub.s3.amazonaws.com/brca_tcga_pub2015.tar.gz'
tcga_name = 'brca_tcga_pub2015'
expression_name = "brca_tcga_pub2015/data_RNA_Seq_v2_expression_median.txt"
clinical_name = "brca_tcga_pub2015/data_clinical_sample.txt"


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
    brca_clin = brca_clin.loc[["PR status by ihc","ER Status By IHC","IHC-HER2"]].rename(index={"PR status by ihc" : "PR","ER Status By IHC":"ER","IHC-HER2":"HER2"})
    relevant_col = brca_clin.T["PR"].isin(['Negative', 'Positive']) & brca_clin.T["ER"].isin(['Negative', 'Positive']) & brca_clin.T["HER2"].isin(['Negative', 'Positive'])
    brca_clin = brca_clin.T[relevant_col].T
    brca = brca.T[relevant_col].T
    # mapping = {'Negative': 0, 'Positive': 1, '[Not Available]': 1, 'Equivocal' : 1, 'Indeterminate' : 1}
    mapping = {'Negative': 0, 'Positive': 1}
    brca_clin.replace(mapping, inplace=True)
    return brca, brca_clin

def download_file_if_not_present(url,fileN):
    if (not os.path.exists(fileN)):
        urllib.request.urlretrieve(url, filename=fileN)

def tcga_read_data():
    """This function is quite convoluted as it downloads one tarfile and extracts two dataframes. The function also caches intermediate files"""
    print("Downloading data ...")
    try:
        os.mkdir(cache_directory)
    except FileExistsError:
        pass
    
    brca_path= os.path.join(cache_directory,  expression_name + ".tsv.gz")
    brca_clin_path = os.path.join(cache_directory,  phenotype_name + ".tsv.gz")
    try:
        brca = pd.read_csv(brca_path, sep="\t",index_col=0)
        brca_clin = pd.read_csv(brca_clin_path, sep="\t",index_col=0)
        return brca, brca_clin
    except:
        pass

    # Download TCGA data
    tcga_path = os.path.join(cache_directory,  tcga_name + ".tar.gz")
    download_file_if_not_present(tcga_url, tcga_path)
    print("Tar-file inplace, extracting tables.")

    # Decompress data into tables
    tf = tarfile.open(tcga_path)
    print(tf.getmembers())
    tf.extract(expression_name, cache_directory)
    tf.extract(clinical_name, cache_directory)

    #  def get_expression_data(self, path, file):
    df = pd.read_csv(os.path.join(cache_directory, expression_name), sep="\t")
    df.dropna(axis=0, how='any', inplace=True)
    df.set_index('Entrez_Gene_Id', inplace=True)
    # df.drop(columns=['Unnamed: 0', 'Entrez_Gene_Id'], inplace=True)
    # df.drop(columns=['Entrez_Gene_Id'], inplace=True)
    df.drop(columns=['Hugo_Symbol'], inplace=True)
    brca = df.reindex(sorted(df.columns), axis=1)

    # get_clinical_data(brca_clin_path,"data_clinical_sample.txt")
    df = pd.read_csv(os.path.join(cache_directory, clinical_name), sep="\t").T
    df.columns = df.loc["Sample Identifier"]
    df.drop(columns=["A unique sample identifier.","STRING","1","SAMPLE_ID"], inplace=True,errors='ignore')
    if 'TCGA-BH-A1ES-01' in df.columns:
        df.drop(columns=['TCGA-BH-A1ES-01'], inplace=True)
    df.drop(index=["Unnamed: 0","#Patient Identifier","Sample Identifier","Other Sample ID"], inplace=True,errors='ignore')
    brca_clin = df.reindex(sorted(df.columns), axis=1)

    # def tcga_preprocess(brca, brca_clin):
    brca.dropna(axis=0, how='any', inplace=True)
    brca = brca.loc[~(brca<=0.0).any(axis=1)].copy()
    entrez2ensembl, drop_list = get_ensembl_ids(list(brca.index))
    # Convert all genenames we can, drop the rest
    brca.drop(index=drop_list, inplace=True)
    brca.rename(index=entrez2ensembl, inplace=True)
    #brca.index =  brca.index.map(str)
    #brca = pd.DataFrame(data=np.log2(brca.values), index=brca.index, columns=brca.columns)
    brca_clin = brca_clin.loc[["PR status by ihc","ER Status By IHC","IHC-HER2"]].rename(index={"PR status by ihc" : "PR","ER Status By IHC":"ER","IHC-HER2":"HER2"})
    relevant_col = brca_clin.T["PR"].isin(['Negative', 'Positive']) & brca_clin.T["ER"].isin(['Negative', 'Positive']) & brca_clin.T["HER2"].isin(['Negative', 'Positive'])
    brca_clin = brca_clin.T[relevant_col].T
    brca = brca.T[relevant_col].T
    # mapping = {'Negative': 0, 'Positive': 1, '[Not Available]': 1, 'Equivocal' : 1, 'Indeterminate' : 1}
    mapping = {'Negative': 0, 'Positive': 1}
    brca_clin.replace(mapping, inplace=True)
    # Put the extracted matrixes to the file cashe, so we do not have to do this again if procedure is repeated.
    brca.to_csv(brca_path, sep="\t")
    brca_clin.to_csv(brca_clin_path, sep="\t")
    return brca, brca_clin

def tcga_read_data_old():
    """This function is quite convoluted as it downloads one tarfile and extracts two dataframes. The function also caches intermediate files"""
    print("Downloading data ...")
    try:
        os.mkdir(cache_directory)
    except FileExistsError:
        pass

    brca_path= os.path.join(cache_directory,  expression_name + ".tsv.gz")
    brca_clin_path = os.path.join(cache_directory,  phenotype_name + ".tsv.gz")
    proc_brca_t = FileTracker(os.path.join(cache_directory,  preprcoc_prefix + expression_name + ".tsv.gz"))
    proc_brca_clin_t = FileTracker(os.path.join(cache_directory,  preprcoc_prefix + phenotype_name + ".tsv.gz"))
    tcga = TcgaTracker('https://cbioportal-datahub.s3.amazonaws.com/brca_tcga_pub2015.tar.gz', cache_directory, "my.tar.gz")
    try:
        brca = proc_brca_t.read_file()
        brca_clin = proc_brca_clin_t.read_file()
    except Exception:
        brca = tcga.get_expression_data(brca_path,"data_RNA_Seq_v2_expression_median.txt")
        brca_clin = tcga.get_clinical_data(brca_clin_path,"data_clinical_sample.txt")
        print("Preprocessing data ...")
        brca, brca_clin = tcga_preprocess(brca, brca_clin)
        proc_brca_t.write_file(brca)
        proc_brca_clin_t.write_file(brca_clin)
    return brca, brca_clin

def tcga_example():
    brca, brca_clin = tcga_read_data()
    print("Run Porch ...")
    activity, eigen_samples, untested = porch.porch_reactome(brca)
    significance = porch.linear_model("Pathway ~ C(PR) + C(ER) + C(PR):C(ER)", activity, brca_clin)
    ## Plot the activity of R-HSA-8931987 in TNBC vs non-TNBC
    sns.set_palette("bright")
    tripple_neg = (brca_clin.T["PR"] == 0) & (brca_clin.T["ER"] == 0) & (brca_clin.T["HER2"] == 0)
    runx1 = activity.loc["R-HSA-8931987",:].T
    runx1 = runx1.drop(["annotation", "set_size"]).to_numpy()

    print(runx1)
    hist, bin_edges = np.histogram(runx1,40)
    fig = plt.figure(figsize=(10,6))
    sns.distplot(runx1[tripple_neg], bins=bin_edges, kde=False,norm_hist=True, color="red")
    sns.distplot(runx1[~tripple_neg], bins=bin_edges, kde=False,norm_hist=True, color="blue")
    plt.legend(labels=['TNBC','Non-TNBC'],loc='upper right')
    plt.ylabel('Density')
    plt.xlabel('Activity of \"RUNX1 regulates estrogen receptor mediated transcription\"')
#    plt.show()
    plt.savefig("rux1.png")

    # Scatterplots
    qv.qvalues(significance,"C(PR)", "q value PR")
    qv.qvalues(significance, "C(ER)","q value ER")
    qv.qvalues(significance, "C(PR):C(ER)", "q value PR:ER")
    fig = plt.figure(figsize=(10,6))
    g = sns.scatterplot(data=significance,x="q value PR",y="q value ER")
    g.set_xscale('log')
    g.set_yscale('log')
    g.set_xlim(10e-80,0.5)
    g.set_ylim(10e-50,0.5)
    plt.savefig("PRvsER.png")

    fig = plt.figure(figsize=(10,6))
    g = sns.scatterplot(data=significance, x="q value PR", y="q value PR:ER")
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

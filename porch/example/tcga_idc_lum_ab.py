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
expression_matrix_name = "idc"
phenotype_name = "idc_clin"
preprcoc_prefix = "proc_"
significance_name  = "significance"
activity_name = "activity"

tcga_url = 'https://cbioportal-datahub.s3.amazonaws.com/brca_tcga_pub2015.tar.gz'
tcga_name = 'brca_tcga_pub2015'
expression_name = "brca_tcga_pub2015/data_RNA_Seq_v2_expression_median.txt"
clinical_a_name = "brca_tcga_pub2015/case_lists/cases_idc_luma.txt"
clinical_b_name = "brca_tcga_pub2015/case_lists/cases_idc_lumb.txt"


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

def read_case_list(path):
    with open(path) as fh:
        for line in fh:
            if line.startswith("case_list_ids"):
                case_list = line.rstrip().split()
                return case_list[1:]

def tcga_tn_preprocess(idc, idc_clin):
    idc.dropna(axis=0, how='any', inplace=True)
    idc = idc.loc[~(idc<=0.0).any(axis=1)].copy()
    entrez2ensembl, drop_list = get_ensembl_ids(list(idc.index))
    # Convert all genenames we can, drop the rest
    idc.drop(index=drop_list, inplace=True)
    idc.rename(index=entrez2ensembl, inplace=True)
    #idc.index =  idc.index.map(str)
    #idc = pd.DataFrame(data=np.log2(idc.values), index=idc.index, columns=idc.columns)
    return idc, idc_clin

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
    
    idc_path= os.path.join(cache_directory,  expression_matrix_name + ".tsv.gz")
    idc_clin_path = os.path.join(cache_directory,  phenotype_name + ".tsv.gz")
    try:
        idc = pd.read_csv(idc_path, sep="\t",index_col=0)
        idc_clin = pd.read_csv(idc_clin_path, sep="\t",index_col=0)
        return idc, idc_clin
    except:
        pass

    # Download TCGA data
    tcga_path = os.path.join(cache_directory,  tcga_name + ".tar.gz")
    download_file_if_not_present(tcga_url, tcga_path)
    print("Tar-file inplace, extracting tables.")

    # Decompress data into tables
    tf = tarfile.open(tcga_path)
    tf.extract(expression_name, cache_directory)
    tf.extract(clinical_a_name, cache_directory)
    tf.extract(clinical_b_name, cache_directory)

    #  def get_expression_data(self, path, file):
    df = pd.read_csv(os.path.join(cache_directory, expression_name), sep="\t")
    df.dropna(axis=0, how='any', inplace=True)
    df.set_index('Entrez_Gene_Id', inplace=True)
    # df.drop(columns=['Unnamed: 0', 'Entrez_Gene_Id'], inplace=True)
    # df.drop(columns=['Entrez_Gene_Id'], inplace=True)
    df.drop(columns=['Hugo_Symbol'], inplace=True)
    idc = df.reindex(sorted(df.columns), axis=1)

    # get_clinical_data(idc_clin_path,"data_clinical_sample.txt")
    samp_a = read_case_list(os.path.join(cache_directory, clinical_a_name))
    samp_b = read_case_list(os.path.join(cache_directory, clinical_b_name))

    samp_a.remove("TCGA-BH-A1ES-01")

    idc_clin = pd.DataFrame(index=["LuminalB"],
                            columns=samp_a+samp_b,
                            data = [[0]*len(samp_a)+[1]*len(samp_b)])
    idc = idc[samp_a+samp_b]    

    idc, idc_clin = tcga_tn_preprocess(idc, idc_clin)

    # Put the extracted matrixes to the file cashe, so we do not have to do this again if procedure is repeated.
    idc.to_csv(idc_path, sep="\t")
    idc_clin.to_csv(idc_clin_path, sep="\t")
    return idc, idc_clin


def tcga_example():
    idc, idc_clin = tcga_read_data()
    print("Run Porch ...")
    activity, eigen_samples, untested = porch.porch_reactome(idc)
    significance = porch.linear_model("Pathway ~ C(LuminalB)", activity, idc_clin)
    ## Regress away proliferation
    print(idc_clin)
    print(activity)
    print(activity.loc["R-HSA-9675126"])
    idc_clin.loc["Proliferation"] = activity.loc["R-HSA-9675126"]
    significance_regressed = porch.linear_model("Pathway ~ C(LuminalB) + Proliferation", activity, idc_clin)
    significance["C(LuminalB)reg"] = significance_regressed["C(LuminalB)"]

    ## Plot the activity of R-HSA-8931987 in TNBC vs non-TNBC
    sns.set_palette("bright")

    # Scatterplots
    qv.qvalues(significance,"C(LuminalB)", "q value Lum A/B")
    qv.qvalues(significance,"C(LuminalB)reg", "q value Lum A/B-Proliferation")
    fig = plt.figure(figsize=(10,6))
    g = sns.scatterplot(data=significance,x="q value Lum A/B",y="q value Lum A/B-Proliferation")
    g.set_xscale('log')
    g.set_yscale('log')
    g.set_xlim(10e-80,0.5)
    g.set_ylim(10e-50,0.5)
    plt.savefig("qAvsB.png")

def main():
    tcga_example()

if __name__ == "__main__":
    main()

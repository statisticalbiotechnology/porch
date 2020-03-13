import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import tarfile
import requests
import os
from biothings_client import get_client
from bioservices import KEGG
import porch
import qvalue as qv

cache_directory = ".porch"
protein_expression_name = "tcell_protein"
metabolite_expression_name = "tcell_metabolite"
preprcoc_prefix = "proc_"
significance_name  = "significance"
activity_name = "activity"

proteomics_data_url = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867416313137-mmc1.xlsx"
metabolomics_data_url = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867416313137-mmc2.xlsx"

class FileCache:
    def __init__(self, file_name):
        self.file_name = file_name

    def get_file_name(self):
        return self.file_name


class UrlFileCache(FileCache):
    def __init__(self, file_name, url):
        self.file_name = file_name
        self.url = url

    def track_dl(self):
        response = requests.get(self.url, stream=True)
        with open(self.file_name, "wb") as handle:
            for data in response.iter_content():
                handle.write(data)

    def get_file_name(self):
        if not os.path.isfile(self.file_name):
            self.track_dl()
        return self.file_name

class TsvFileTracker(FileCache):
    def __init__(self, file_name, filler):
        self.file_name = file_name
        self.filler = filler

    def get_file_name(self):
        return self.file_name

    def read_file(self):
        if not os.path.isfile(self.get_file_name()):
            df = self.filler()
            self.write_file(df)
            return df
        return pd.read_csv(self.get_file_name(), sep="\t", index_col=0)

    def write_file(self,df):
        df.to_csv(self.file_name, sep="\t")

def one_row_per_proteoform(group):
    # From https://stackoverflow.com/questions/13050003/apply-function-to-pandas-dataframe-that-can-return-multiple-rows
    row = group.iloc[0]
#    proteoforms = row['Protein IDs'].split(';')
    proteoforms = str(row.name).split(';')
    copies = len(proteoforms)
    content = {'ProteinID' : proteoforms}
    row_dict = row.to_dict()
    for item in row_dict: content[item] = [row[item]] * copies
    return pd.DataFrame(content)

def one_row_per_compound_convert(group, map_kegg_chebi):
    # From https://stackoverflow.com/questions/13050003/apply-function-to-pandas-dataframe-that-can-return-multiple-rows
    row = group.iloc[0]

    chebis = []
    for kegg_id in str(row.name).split(','):
        name = 'cpd:' + kegg_id
        if name in map_kegg_chebi:
            chebis += [map_kegg_chebi[name].split(':')[1]] # The map returns compond ids like 'chebi:5292'
    copies = len(chebis)
    if copies < 1:
        return None
    content = {'MetaboliteID' : chebis}
    row_dict = row.to_dict()
    for item in row_dict: content[item] = [row[item]] * copies
    return pd.DataFrame(content)


def tcell_read_metabolomics_data():
    """This function is quite convoluted as it downloads an excelfile from a publication and extracts a dataframe, idexed by chebi. The function also caches intermediate files"""
    tcell_metabol_xls = UrlFileCache(os.path.join(cache_directory,  metabolite_expression_name + ".xlsx"), metabolomics_data_url)
    metabolomics_df = pd.read_excel(tcell_metabol_xls.get_file_name(), sheet_name = "normalized by sample mean", index_col=0, usecols="A,C:HN", skiprows = [0])
    metabolomics_df.index.name = "KEGG_ID"
    metabolomics_df = metabolomics_df.apply(np.exp2)    # The excel data is in log2 space, return it to normal
    k = KEGG(verbose=False)
    map_kegg_chebi = k.conv("chebi", "compound")
    metabolomics_df = metabolomics_df.groupby("KEGG_ID", group_keys=False).apply(lambda x: one_row_per_compound_convert(x, map_kegg_chebi)).reset_index(drop=True)
    metabolomics_df.set_index("MetaboliteID", inplace=True)
    return metabolomics_df

def tcell_read_metabolomics_frames():
    try:
        os.mkdir(cache_directory)
    except FileExistsError:
        pass
    proc_tcell_t = TsvFileTracker(os.path.join(cache_directory, metabolite_expression_name + ".tsv.gz"), tcell_read_metabolomics_data)
    metabolomics_df = proc_tcell_t.read_file()
    values = []
    for coln in metabolomics_df.columns:
        if "non act" in coln:
            values += [-1.]
        elif "ON" in coln:
            values += [0.]
        elif "act 3h" in coln:
            values += [3.]
        elif "act 12h" in coln:
            values += [12.]
        elif "act 14h" in coln:
            values += [12.]
        elif "act 24h" in coln:
            values += [24.]
        elif "act 2d" in coln:
            values += [48.]
        elif "act 3d" in coln:
            values += [72.]
        elif "act 4d" in coln:
            values += [96.]
        else:
            print(coln)
    phenotype_df = pd.DataFrame(columns=metabolomics_df.columns, data=[values], index=["Time"])
    return phenotype_df, metabolomics_df


def tcell_read_proteomics_data():
    """This function is quite convoluted as it downloads an excelfile from a publication and extracts a dataframe. The function also caches intermediate files"""

    tcell_prot_xls = UrlFileCache(os.path.join(cache_directory,  protein_expression_name + ".xlsx"),proteomics_data_url)
    proteomics_df = pd.read_excel(tcell_prot_xls.get_file_name(), sheet_name = "Data", index_col=0, usecols="A,D:U")
#    proteomics_df = pd.read_excel(tcell_prot_xls.get_file_name(), sheet_name = "Data", index_col=0, usecols="A,V:AM")
    proteomics_df = proteomics_df - proteomics_df.mean()    # Normalize by subtracting column mean
    proteomics_df = proteomics_df.apply(np.exp2)    # The excel data is in log2 space, return it to normal
    proteomics_df = proteomics_df.groupby("Protein IDs", group_keys=False).apply(one_row_per_proteoform).reset_index(drop=True)
    proteomics_df.set_index("ProteinID", inplace=True)
    return proteomics_df

def tcell_read_proteomics_frames():
    try:
        os.mkdir(cache_directory)
    except FileExistsError:
        pass
    proc_tcell_t = TsvFileTracker(os.path.join(cache_directory, protein_expression_name + ".tsv.gz"),tcell_read_proteomics_data)
    proteomics_df = proc_tcell_t.read_file()
    phenotype_df = pd.DataFrame(columns=proteomics_df.columns, data=[[0]*7+[12]*3+[24]+[48]*2+[72]*3+[96]*2], index=["Time"])
    return phenotype_df, proteomics_df

def tcell_example():
    print("Downloading data ...")
    p_phenotype_df, proteomics_df = tcell_read_proteomics_frames()
    m_phenotype_df, metabolomics_df = tcell_read_metabolomics_frames()
    print(proteomics_df.describe())
    print("Factorize data ...")
    p_activity_df,untested = porch.porch_reactome(proteomics_df, organism = "HSA", gene_anot = "UniProt")
    print(p_activity_df)
    m_activity_df,untested = porch.porch_reactome(metabolomics_df, organism = "HSA", gene_anot = "ChEBI")
    print(m_activity_df)
    print("Significance Testing ...")
    p_significance = porch.linear_model("Pathway ~ C(Time)", p_activity_df, p_phenotype_df)
    m_significance = porch.linear_model("Pathway ~ C(Time)", m_activity_df, m_phenotype_df)
    print("Multiple Hypothesis correction ...")
    qv.qvalues(p_significance,"C(Time)", "q_value_Time")
    print(m_significance)
    qv.qvalues(m_significance,"C(Time)", "q_value_Time")
    fig = plt.figure(figsize=(10,6))
    p_significance["-log10(q)"] = -np.log10(p_significance["q_value_Time"])
    g = sns.distplot(p_significance["-log10(q)"], rug=True, kde=False)
    plt.savefig("p_tcell-qtime.png")
    plt.show()
    fig = plt.figure(figsize=(10,6))
    m_significance["-log10(q)"] = -np.log10(m_significance["q_value_Time"])
    g = sns.distplot(m_significance["-log10(q)"], rug=True, kde=False)
    plt.savefig("m_tcell-qtime.png")
    plt.show()
    p_significance.sort_values(by=['q_value_Time'], inplace=True, ascending=True)
    m_significance.sort_values(by=['q_value_Time'], inplace=True, ascending=True)
    print("The most significant proteomics pathways are:")
    print(p_significance.head(n=5))
    print("The most significant metabolomics pathways are:")
    print(m_significance.head(n=10))
    most = p_significance.iloc[0:5:1].index
    p_phenotype_df = p_phenotype_df.append(p_activity_df.loc[most]).T.reset_index()
    print(p_phenotype_df)
    out_df = pd.melt(p_phenotype_df,id_vars=["Time","index"],value_vars=most, var_name='Pathway', value_name='Activity')
    sns.lineplot(data=out_df, x="Time", y="Activity", hue="Pathway")
    plt.savefig("p_tcell-qtime-top.png")
    plt.show()
    most = m_significance.iloc[0:10:1].index
    m_phenotype_df = m_phenotype_df.append(m_activity_df.loc[most]).T.reset_index()
    print(m_phenotype_df)
    out_df = pd.melt(m_phenotype_df,id_vars=["Time","index"],value_vars=most, var_name='Pathway', value_name='Activity')
    sns.lineplot(data=out_df, x="Time", y="Activity", hue="Pathway")
    plt.savefig("m_tcell-qtime-top.png")
    plt.show()

def main():
    tcell_example()

if __name__ == "__main__":
    main()

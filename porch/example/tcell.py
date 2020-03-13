import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats
import numpy as np
import pandas as pd
import tarfile
import requests
import os
from biothings_client import get_client
from bioservices import KEGG
import porch
import porch.qvalue as qv

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
    #metabolomics_df = pd.read_excel(tcell_metabol_xls.get_file_name(), sheet_name = "normalized by sample mean", index_col=0, usecols="A,C:HN", skiprows = [0])
    for col in metabolomics_df.columns:
        # Average all technical replicates (Named by trailing ".1")
        if len(col.split('.'))>1 and col.split('.')[1] == "1":
            remcol = col.split('.')[0]
            metabolomics_df[remcol] = scipy.stats.gmean(metabolomics_df[[remcol,col]],axis=1)
            metabolomics_df.drop(col, axis=1, inplace=True)
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
    values,cols = [],[]
    for coln in metabolomics_df.columns:
        if "non act" in coln:
            time = -1.
        elif "ON" in coln:
            time = 0.
        elif "act 3h" in coln:
            time = 3.
        elif "act 12h" in coln:
            time = 12.
        elif "act 14h" in coln:
            time = 14.
        elif "act 24h" in coln:
            time = 24.
        elif "act 2d" in coln:
            time = 48.
        elif "act 3d" in coln:
            time = 72.
        elif "act 4d" in coln:
            time = 96.
        else:
            print(coln)
        dish = coln.split('-')[0]
        rep = coln.split('-')[2].replace(" ","")
        values += [time]
        cols += ['_'.join([dish,str(int(time)),rep])]
    phenotype_df = pd.DataFrame(columns=cols, data=[values], index=["Time"])
    metabolomics_df.columns = cols
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
    values,cols = [],[]
    for coln in proteomics_df.columns:
        if "notact" in coln:
            time = 0.
        elif "act12h" in coln:
            time = 12.
        elif "act24h" in coln:
            time = 24.
        elif "act48h" in coln:
            time = 48.
        elif "act72h" in coln:
            time = 72.
        elif "act96h" in coln:
            time = 96.
        else:
            print(coln)
        print(coln)
        not_sure = coln.split('_')[0].replace("q","")
        rep = int(coln.split('_')[3].replace(" ",""))
        if rep<18:
            dish = int(rep/5)+2
            rep = rep%5+1
        elif rep<26:
            dish = rep%3+2
            rep = 3
        elif rep<31:
            dish = (rep-1)%3+2
            rep = 3
        else:
            dish = (rep-2)%3+2
            rep = 3
        values += [time]
        cols += ['_'.join([str(dish),str(int(time)),str(rep)])]
    proteomics_df.columns = cols
    phenotype_df = pd.DataFrame(columns=cols, data=[values], index=["Time"])
    return phenotype_df, proteomics_df

def tcell_example():
    print("Downloading data ...")
    p_phenotype_df, proteomics_df = tcell_read_proteomics_frames()
    m_phenotype_df, metabolomics_df = tcell_read_metabolomics_frames()
    #print(proteomics_df.describe())
    print(proteomics_df.columns)
    print(metabolomics_df.columns)
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
    print(p_significance.head(n=10))
    print("The most significant metabolomics pathways are:")
    print(m_significance.head(n=20))
    most = p_significance.iloc[0:5:1].index
    p_joint_df = p_phenotype_df.append(p_activity_df.loc[most]).T.reset_index()
    out_df = pd.melt(p_joint_df,id_vars=["Time","index"],value_vars=most, var_name='Pathway', value_name='Activity')
    sns.lineplot(data=out_df, x="Time", y="Activity", hue="Pathway")
    plt.savefig("p_tcell-qtime-top.png")
    plt.show()
    most = m_significance.iloc[0:10:1].index
    m_phenotype_df = m_phenotype_df.append(m_activity_df.loc[most]).T.reset_index()
    out_df = pd.melt(m_phenotype_df,id_vars=["Time","index"],value_vars=most, var_name='Pathway', value_name='Activity')
    sns.lineplot(data=out_df, x="Time", y="Activity", hue="Pathway")
    plt.savefig("m_tcell-qtime-top.png")
    plt.show()
    print("MultiOmics analysis")
    multiomics_df = pd.concat([proteomics_df,metabolomics_df],axis=0,join="inner")
    multi_phenotype_df = p_phenotype_df[multiomics_df.columns]
    multi_activity_df, untested = porch.porch_multi_reactome(multiomics_df,[["HSA","UniProt"], ["HSA","ChEBI"]])
    multi_significance = porch.linear_model("Pathway ~ C(Time)", multi_activity_df, multi_phenotype_df)
    qv.qvalues(multi_significance,"C(Time)", "q_value_Time")
    multi_significance["-log10(q)"] = -np.log10(multi_significance["q_value_Time"])
    print("The most significant multiomics pathways are:")
    print(multi_significance.head(n=500))
    for s in range(0,200,5):
        most = multi_significance.iloc[s:s+5:1].index
        multi_joint_df = p_phenotype_df.append(p_activity_df.loc[most]).T.reset_index()
        out_df = pd.melt(multi_joint_df,id_vars=["Time","index"],value_vars=most, var_name='Pathway', value_name='Activity')
        sns.lineplot(data=out_df, x="Time", y="Activity", hue="Pathway")
        plt.savefig("multi_tcell-qtime-{}{}.png".format(s,s+5))
        plt.show()

def main():
    tcell_example()

if __name__ == "__main__":
    main()

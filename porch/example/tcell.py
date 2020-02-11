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
protein_expression_name = "tcell_protein"
preprcoc_prefix = "proc_"
significance_name  = "significance"
activity_name = "activity"

proteomics_data_url = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867416313137-mmc1.xlsx"

class FileCashe:
    def __init__(self, file_name):
        self.file_name = file_name

    def get_file_name(self):
        return self.file_name


class UrlFileCashe(FileCashe):
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

class TsvFileTracker(FileCashe):
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

def tcell_read_proteomics_data():
    """This function is quite convoluted as it downloads one tarfile and extracts two dataframes. The function also cashes intermediate files"""

    tcell_prot_xls = UrlFileCashe(os.path.join(cashe_directory,  protein_expression_name + ".xlsx"),proteomics_data_url)
    proteomics_df = pd.read_excel(tcell_prot_xls.get_file_name(), sheet_name = "Data", index_col=0, usecols="A,D:U")
#    proteomics_df = pd.read_excel(tcell_prot_xls.get_file_name(), sheet_name = "Data", index_col=0, usecols="A,V:AM")
    proteomics_df = proteomics_df.groupby("Protein IDs", group_keys=False).apply(one_row_per_proteoform).reset_index(drop=True)
    proteomics_df.set_index("ProteinID", inplace=True)
    proteomics_df.apply(np.exp2, inplace=True)    # The excel data is in log2 space, return it to normal
    return proteomics_df

def tcell_read_data():
    try:
        os.mkdir(cashe_directory)
    except FileExistsError:
        pass
    proc_tcell_t = TsvFileTracker(os.path.join(cashe_directory, protein_expression_name + ".tsv.gz"),tcell_read_proteomics_data)
    proteomics_df = proc_tcell_t.read_file()
    phenotype_df = pd.DataFrame(columns=proteomics_df.columns, data=[[0]*7+[12]*3+[24]+[48]*2+[72]*3+[96]*2], index=["Time"])
    return phenotype_df, proteomics_df

def tcell_example():
    print("Downloading data ...")
    phenotype_df, proteomics_df = tcell_read_data()
    print(proteomics_df.describe())
    print("Factorize data ...")
    activity_df,untested = porch.porch_reactome(proteomics_df, organism = "HSA", gene_anot = "UniProt")
    print(activity_df)
    print("Significance Testing ...")
    significance = porch.linear_model("Pathway ~ C(Time)", activity_df, phenotype_df)
    print("Multiple Hypothesis correction ...")
    qv.qvalues(significance,"C(Time)", "q_value_Time")
    fig = plt.figure(figsize=(10,6))
    significance["-log10(q)"] = -np.log10(significance["q_value_Time"])
    g = sns.distplot(significance["-log10(q)"], rug=True, kde=False)
#   g.set_xscale('log')
#    g.set_xlim(10e-80,0.5)
#    g.set_ylim(10e-50,0.5)
    plt.savefig("tcell-qtime.png")
    plt.show()
    significance.sort_values(by=['q_value_Time'], inplace=True, ascending=True)
    print("The most significant pathways are:")
    print(significance.head(n=5))
    most = significance.iloc[0:5:1].index
    phenotype_df = phenotype_df.append(activity_df.loc[most]).T.reset_index()
    print(phenotype_df)
    out_df = pd.melt(phenotype_df,id_vars=["Time","index"],value_vars=most, var_name='Pathway', value_name='Activity')
    sns.lineplot(data=out_df, x="Time", y="Activity", hue="Pathway")
    plt.savefig("tcell-qtime-top.png")
    plt.show()

def main():
    tcell_example()

if __name__ == "__main__":
    main()

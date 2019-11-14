# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from numpy.linalg import svd
from sklearn.preprocessing import StandardScaler
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
import urllib.request
import os.path


def porch(expression_df, phenotype_df, geneset_df,
    gene_column = "gene", set_column = "pathway",
    tests = ["Pathway ~ C(Case)"]):
    """
    This is a the central routine of porch. It calculates pathway activities from the expression values of analytes,
    with a grouping given by a pathway definition. If so specified, it tests the pathway activities with a set of tests specified by the user

    Args:
        expression_df (pd.DataFrame): The DataFrame of the expression values we analyse. These values are logtransformed and subsequently standardized befor analysis
        phenotype_df (pd.DataFrame): The DataFrame of the phenotypic variables of each sample, that we want to test against.
        geneset_df (pd.DataFrame): The DataFrame of the pathway definitions.
        gene_column (str): The name of the column within geneset_df containing names of analytes.
        set_column (str): The name of the column within geneset_df containing names of pathways.
        tests (list): List of specification of tests that should be performed.

    Returns:
        A tupple of pandas DataFrames results_df, evaluation_df
        results_df, contains the output of the significance tests
        evaluation_df, contains the pathway activities
    """
    phenotype_df =  phenotype_df[[ col for col in phenotype_df.columns  if col in expression_df.columns ]]
    evaluation_df = phenotype_df.copy()
    results_df = pd.DataFrame()
    set_df = geneset_df[[gene_column, set_column]]
    set_of_all_genes = set(expression_df.index)
    results,setnames,tested_vars = [],[],[]
    for setname, geneset in set_df.groupby([set_column]):
        print("Processing " + setname)
        genes = list(set(geneset[gene_column].tolist()) & set_of_all_genes)
        expr = expression_df.loc[genes,phenotype_df.columns]
        expr.dropna(axis=0, how='any', inplace=True)
        expr = expr.loc[~(expr<=0.0).any(axis=1)]
        if expr.shape[0]>1:
            expr_stand = pd.DataFrame(index=expr.index,    columns=expr.columns,data=StandardScaler().fit_transform(np.log(expr.values.T)).T)
            U, S, Vt = svd(expr_stand.values,full_matrices=False)
            eigen_genes = (Vt.T)[:,0]
            evaluation_df.loc["Pathway"] = eigen_genes
            result, tested_vars = [],[]
            for test in tests:
                lm = ols(test, evaluation_df.T).fit()
                pvals = anova_lm(lm)["PR(>F)"].T.iloc[:-1]
                result += list(pvals.values)
                tested_vars += list(pvals.index)
            results += [result]
            setnames += [setname]
    results_cols = pd.MultiIndex.from_product([tests,tested_vars],names=["Test","Variable"])
    results_df = pd.DataFrame(data=results, columns=results_cols,index=setnames)
    return results_df,evaluation_df

def porch_reactome(expression_df, phenotype_df, organism = "HSA", tests = ["Pathway ~ C(Case)"]):
    "This is a function"
    reactome_df = get_reactome_df(organism)
    return porch(expression_df, phenotype_df, reactome_df,
        "gene", "reactome_id", tests)


def download_file(path, url):
    "This function downloads a file, path, from an url, if the file does not already exist"
    if not os.path.isfile(path):
        stream = urllib.request.urlopen(url)
        with open(path,'wb') as output:
            output.write(stream.read())
    return path

reactome_fn = "Ensembl2Reactome_All_Levels.txt"
#reactome_fn = "UniProt2Reactome_All_Levels.txt"
reactome_path = ".porch/" + reactome_fn
reactome_url = "https://reactome.org/download/current/" + reactome_fn

def get_reactome_df(organism = "HSA"):
    reactome_df = pd.read_csv(download_file(reactome_path, reactome_url),
                        sep='\t',
                        header=None,
                        usecols=[0,1,3],
                        names=["gene","reactome_id","reactome_name"])
    organism = "R-" + organism
    reactome_df = reactome_df[reactome_df["reactome_id"].str.startswith(organism) ]
    return reactome_df

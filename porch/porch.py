# -*- coding: utf-8 -*-
import multiprocessing
import numpy as np
import pandas as pd
from numpy.linalg import svd
from sklearn.preprocessing import StandardScaler
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
import urllib.request
import os.path
import sys


def porch_single_process(expression_df, geneset_df,
    gene_column = "gene", set_column = "pathway",
    test = "Pathway ~ C(Case)"):
    """
    Calculates pathway activities from the expression values of analytes,
    with a grouping given by a pathway definition.
    This call is not using parallel processing, mostly for debugging purposes.

    Args:
        expression_df (pd.DataFrame): The DataFrame of the expression values we analyse. These values are logtransformed and subsequently standardized before analysis
        geneset_df (pd.DataFrame): The DataFrame of the pathway definitions.
        gene_column (str): The name of the column within geneset_df containing names of analytes.
        set_column (str): The name of the column within geneset_df containing names of pathways.

    Returns:
        results_df, untested
        A pandas DataFrames results_df, containing the output of the significance tests
        untested, a list of the pathway that were not possible to test, due to shortage of data in expression_df.
    """
    phenotype_df =  phenotype_df[[ col for col in phenotype_df.columns  if col in expression_df.columns]]
    expression_df = expression_df[phenotype_df.columns]
    results_df = pd.DataFrame()
    set_df = geneset_df[[gene_column, set_column]]
    set_of_all_genes = set(expression_df.index)
    results, setnames, untested, activities, tested_vars = [], [], [], [], None
    for setname, geneset in set_df.groupby([set_column]):
        genes = list(set(geneset[gene_column].tolist()) & set_of_all_genes)
        proc = porch_proc
        setname, result, activity, variables = proc(setname, genes, expression_df, phenotype_df,test)
        if result:
            results += [result]
            setnames += [setname]
            activities += [activity]
            if not tested_vars:
                tested_vars = variables
        else:
            untested += [setname]
    results_df = pd.DataFrame(data=results, columns=tested_vars,index=setnames)
    activity_df = pd.DataFrame(data=activities, columns=expression_df.columns, index=setnames)
    return results_df, activity_df, untested


def porch(expression_df, geneset_df,
    gene_column = "gene", set_column = "pathway",
    test = "Pathway ~ C(Case)"):
    """
    Calculates pathway activities from the expression values of analytes,
    with a grouping given by a pathway definition.

    Args:
        expression_df (pd.DataFrame): The DataFrame of the expression values we analyse. These values are logtransformed and subsequently standardized befor analysis
        geneset_df (pd.DataFrame): The DataFrame of the pathway definitions.
        gene_column (str): The name of the column within geneset_df containing names of analytes.
        set_column (str): The name of the column within geneset_df containing names of pathways.

    Returns:
        results_df, untested
        A pandas DataFrames results_df, containing the output of the significance tests
        untested, a list of the pathway that were not possible to test, due to shortage of data in expression_df.
    """
    results_df = pd.DataFrame()
    set_df = geneset_df[[gene_column, set_column]]
    set_of_all_genes = set(expression_df.index)
    call_args = []
    for setname, geneset in set_df.groupby([set_column]):
        genes = list(set(geneset[gene_column].tolist()) & set_of_all_genes)
        call_args += [(setname, genes, expression_df, test)]
    print("Processing with {} parallel processes".format(os.cpu_count()), file=sys.stderr)
    setnames, untested, activities = [], [], []
    with multiprocessing.Pool() as executor:
        for setname, activity in executor.starmap(porch_proc,  call_args):
            if activity is None:
                untested += [setname]
            else:
                setnames += [setname]
                activities += [activity]
    activity_df = pd.DataFrame(data=activities, columns=expression_df.columns, index=setnames)
    return activity_df, untested


def porch_proc(setname, genes, expression_df, test,keep_feature_stdv=True):
    """ Core processing node of porch. Takes the analysis from expression values to significance testing. """
    print("Decomposing " + setname, file=sys.stderr)
    expr = expression_df.loc[genes]
    expr.dropna(axis=0, how='any', inplace=True)
    expr = expr.loc[~(expr<=0.0).any(axis=1)]
    if expr.shape[0]>2:
        standardizer = StandardScaler(with_std=keep_feature_stdv)
        log_data = np.log(expr.values.T.astype(float))
        standard_log_data = standardizer.fit_transform(log_data).T
        eigen_genes, _ = decomposition_method(standard_log_data)
        return setname, eigen_genes
    else:
        print("Not enough data to evaluate " + setname, file=sys.stderr)
        return setname, None

def porch_reactome(expression_df, organism = "HSA"):
    "Download the Reactome database and subsequently call porch"
    reactome_df = get_reactome_df(organism)
#    return porch_single_process(expression_df, phenotype_df, reactome_df, "gene", "reactome_id", test)
    return porch(expression_df, reactome_df,
        "gene", "reactome_id")


def svd_decomposition(data):
    U, S, Vt = svd(data, full_matrices=False)
    eigen_genes = (Vt.T)[:,0]
    eigen_samples = U[:,0]
    return eigen_genes, eigen_samples

decomposition_method = svd_decomposition

def linear_model(test,activity_df,phenotype_df):
    """
    Applies a linear model, test, that is a function of variables in phenotype_df,
    to each row in the activity_df.

    Args:
        activity_df (pd.DataFrame): The DataFrame of the pathway activity values we analyse.
        phenotype_df (pd.DataFrame): The DataFrame containing any sample oriented variables that are included in the model.
        test (str): linear model that should be tested. The model should contain the variable Pathway, that will be replaced with each pathway's activity.

    Returns:
        results_df
        A pandas DataFrames results_df, containing the output of the significance tests
    """
    expression_df = activity_df.copy()
    phenotype_df =  phenotype_df[[ col for col in phenotype_df.columns  if col in expression_df.columns]]
    expression_df = expression_df[phenotype_df.columns]
    return expression_df.apply(applicable_linear_model,
        axis=1, result_type='reduce',
        args=(test,phenotype_df))

def applicable_linear_model(row,test,phenotype_df):
    phenotype_df.loc["Pathway"] = row.values
    lm = ols(test, phenotype_df.T).fit()
    pvals = anova_lm(lm)["PR(>F)"].T.iloc[:-1]
    return pvals


def download_file(path, url):
    "This function downloads a file, path, from an url, if the file is not already cashed"
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

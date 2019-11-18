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


def porch_single_process(expression_df, phenotype_df, geneset_df,
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
        results_df, activity_df, untested
        A pandas DataFrames results_df, containing the output of the significance tests
        A pandas DataFrames activity_df, containing the pathway activities for each sample
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
        setname, result, activity, variables = porch_proc(setname, genes, expression_df, phenotype_df,tests)
        if result:
            results += [result]
            setnames += [setname]
            activities += [activity]
            if not tested_vars:
                tested_vars = variables
        else:
            untested += [setname]
    results_cols = pd.MultiIndex.from_tuples(tested_vars, names=["Test","Variable"])
    results_df = pd.DataFrame(data=results, columns=results_cols,index=setnames)
    activity_df = pd.DataFrame(data=activities, columns=expression_df.columns, index=setnames)
    return results_df, activity_df, untested


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
        results_df, activity_df, untested
        A pandas DataFrames results_df, containing the output of the significance tests
        A pandas DataFrames activity_df, containing the pathway activities for each sample
        untested, a list of the pathway that were not possible to test, due to shortage of data in expression_df.
    """
    phenotype_df =  phenotype_df[[ col for col in phenotype_df.columns  if col in expression_df.columns]]
    expression_df = expression_df[phenotype_df.columns]
    results_df = pd.DataFrame()
    set_df = geneset_df[[gene_column, set_column]]
    set_of_all_genes = set(expression_df.index)
    call_args = []
    for setname, geneset in set_df.groupby([set_column]):
        genes = list(set(geneset[gene_column].tolist()) & set_of_all_genes)
        call_args += [(setname, genes, expression_df, phenotype_df,tests)]
    print("Processing with {} parallel processes".format(os.cpu_count()), file=sys.stderr)
    results, setnames, untested, activities, tested_vars = [], [], [], [], None
    with multiprocessing.Pool() as executor:
        for setname, result, activity, variables in executor.starmap(porch_proc,  call_args):
            if result:
                results += [result]
                setnames += [setname]
                activities += [activity]
                if not tested_vars:
                    tested_vars = variables
            else:
                untested += [setname]
    results_cols = pd.MultiIndex.from_tuples(tested_vars, names=["Test","Variable"])
    results_df = pd.DataFrame(data=results, columns=results_cols,index=setnames)
    activity_df = pd.DataFrame(data=activities, columns=expression_df.columns, index=setnames)
    return results_df, activity_df, untested


def porch_proc(setname, genes, expression_df, phenotype_df,tests,keep_feature_stdv=True):
    """ Core processing node of porch. Takes the analysis from expression values to significance testing. """
    print("Processing " + setname, file=sys.stderr)
    evaluation_df = phenotype_df.copy()
    expr = expression_df.loc[genes]
    expr.dropna(axis=0, how='any', inplace=True)
    expr = expr.loc[~(expr<=0.0).any(axis=1)]
    if expr.shape[0]>2:
        standardizer = StandardScaler(with_std=keep_feature_stdv)
        log_data = np.log(expr.values.T)
        standard_log_data = standardizer.fit_transform(log_data).T
        U, S, Vt = svd(standard_log_data, full_matrices=False)
        eigen_genes = (Vt.T)[:,0]
        evaluation_df.loc["Pathway"] = eigen_genes
        result, tested_vars = [],[]
        for test in tests:
            lm = ols(test, evaluation_df.T).fit()
            pvals = anova_lm(lm)["PR(>F)"].T.iloc[:-1]
            result += list(pvals.values)
            tested_vars += [(test,var) for var in list(pvals.index)]
        return setname, result, evaluation_df.loc["Pathway"].values, tested_vars
    else:
        print("Not enough data to evaluate " + setname, file=sys.stderr)
        return setname, None, None, None


def porch_reactome(expression_df, phenotype_df, organism = "HSA", tests = ["Pathway ~ C(Case)"]):
    "This is a function"
    reactome_df = get_reactome_df(organism)
    return porch_single_process(expression_df, phenotype_df, reactome_df,
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

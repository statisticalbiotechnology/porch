# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from numpy.linalg import svd
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm

def load_dataframe(path, url, columns, indexes):
    """ This is a function"""
    try:
        tf = tarfile.open(path)
    except:
        track_dl(url, path)
        tf = tarfile.open(path)
    return tf

def porch(expression_df, phenotype_df, geneset_df, gene_column = "gene", set_column = "pathway", tests = ["Pathway ~ C(Case)"]):
    """This is a the central routine of porch. It calculates pathway activities from the expression values of analytes,
     with a grouping given by a pathway definition. If so specified, it tests the pathway activities with a set of user
     specied tests
     : param expression_df: The DataFrame of the expression values we analyse. These values are logtransformed and subsequently standardized befor analysis
     : param phenotype_df: The DataFrame of the phenotypic variables of each sample, that we want to test against.
     : param geneset_df: The DataFrame of the pathway definitions.
     : param gene_column: The name of the column within geneset_df containing names of analytes.
     : param set_column: The name of the column within geneset_df containing names of pathways.
     : param tests: List of specification of tests that should be performed.
     """
    phenotypes = phenotype_df.columns
    phenotypes_bool = geneset.columns in phenotypes
    evaluation_df = phenotype_df.copy()
    results_df = pd.DataFrame()
    for setname, geneset in df.groupby([set_column]):
        genes = list(geneset[gene_column])
        results,setnames = [],[]
        if len(genes)>1:
            expr = geneset_df.loc[genes,phenotype_bool]
            expr_stand = pd.DataFrame(index=expr.index, columns=expr.columns,data=StandardScaler().fit_transform(np.log(expr.values.T)).T)
            U, S, Vt = svd(expr_stand.values,full_matrices=False)
            eigen_genes = (Vt.T)[:,0]
            evaluation_df.loc[setname] = eigen_genes
            result, columns = [],[]
            for test in tests:
                test.replace("Pathway",setname)
                lm = ols(test, evaluation_df).fit()
                pvals = anova_lm(lm)["PR(>F)"].T.iloc[:, :-1]
                result += list(pvals.values)
                columns += list(pvals.columns)
            results += [result]
            setnames += [setname]
    results_cols = MultiIndex.from_product([tests,columns],names=["Test","Variable"])
    results_df = pd.DataFrame(data=results, columns=results_cols,index=setnames)
    return results_df,evaluation_df

def porch_reactome(expression_df, phenotype_df, tests = ["Pathway ~ C(Case)"]):
    """ This is a function"""
    pass

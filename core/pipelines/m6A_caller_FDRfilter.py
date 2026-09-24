
"""Cong Liu Yi_lab, Peking University"""
"""Feb, 2021"""
"""Email: liucong-1112@pku.edu.cn"""
"""Usage: This program is used for FDR correction"""
"""Input: [.txt]"""



import pandas as pd
import argparse
from statsmodels.stats.multitest import multipletests
import os
import sys
import time
import math
from time import strftime

OUTPUT_COLUMNS = ['Chr', 'Sites', 'Strand', 'Gene', 'CR', 'AGcov', 'Acov',
                  'Genecov', 'Ratio', 'Pvalue', 'P_adjust']


def get_data(file1,output,adjustP):
    size = os.path.getsize(file1)
    if size != 0:
        df1_t = pd.read_csv(file1, sep='\t', header=None)
        if df1_t.shape[1] != 18:
            raise ValueError('Candidate input must contain exactly 18 tab-separated columns.')
        """df1"""
        pvalue1 = df1_t.iloc[:, 14]
        if not all(math.isfinite(p) and 0 <= p <= 1 for p in pvalue1):
            raise ValueError('Candidate p-values must be finite numbers between 0 and 1.')
        adjustP1 = multipletests(pvalue1, method='fdr_bh')[1]
        x1=pd.DataFrame({'AdjustPvalue':adjustP1})
        df1_t2 = pd.concat([df1_t,x1],axis=1)
        df1 = df1_t2[df1_t2['AdjustPvalue'] < adjustP]
        df1_total_csv = get_cov_ratio(df1, 'm6A.sites')
        df1_total_csv.to_csv(output + '.csv', index=False, sep="\t")
        sys.stderr.write("[%s] FDR filtering finished successfully!\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
    else:
        pd.DataFrame(columns=OUTPUT_COLUMNS).to_csv(output + '.csv', index=False, sep='\t')
        sys.stderr.write("[%s] No candidate sites; wrote a header-only FDR table successfully.\n" %
                         (strftime("%Y-%m-%d %H:%M:%S", time.localtime())))


def get_cov_ratio(df,sam):
    chr = df.iloc[:,0]
    sites = df.iloc[:,1]
    strand = df.iloc[:,2]
    gene = df.iloc[:,3].fillna('ELSE')
    transcript = df.iloc[:,5].fillna('ELSE')
    nonCR = df.iloc[:,8]
    AGcov = df.iloc[:, 11]
    Acov = df.iloc[:, 12]
    ratio = df.iloc[:,13]
    normedRatio = df.iloc[:, 13]*(1-df.iloc[:, 8])
    Pvalue = df.iloc[:,14]
    Genecov = df.iloc[:, 16]
    p_adjust = df.iloc[:,18]
    pd_t = pd.DataFrame({'Chr':chr,'Sites':sites,'Strand':strand,'Gene':gene,'CR':1-nonCR,'AGcov':AGcov,'Acov':Acov,'Genecov':Genecov,\
     'Ratio':ratio,'Pvalue':Pvalue,'P_adjust':p_adjust})
    return pd_t


if __name__ == "__main__":
    description = """
	"""
    parser = argparse.ArgumentParser(prog="get_commonSites", fromfile_prefix_chars='@', description=description,
                                     formatter_class=argparse.RawTextHelpFormatter)
    # Require
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-i", "--file", dest="file", required=True, help="file1")
    group_required.add_argument("-o", "--output", dest="output", required=True,
                                help="output prefix; result is [prefix].csv (tab-separated)")
    group_required.add_argument("-adp", "--adjustpvalue", dest="adjustpvalue", default=0.005, type=float, help="adjustpvalue, default=0.005")

    options = parser.parse_args()
    if not 0 <= options.adjustpvalue <= 1:
        parser.error('adjustpvalue must be between 0 and 1')

    get_data(options.file,options.output,options.adjustpvalue)





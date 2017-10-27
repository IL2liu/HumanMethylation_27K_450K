import pandas as pd
import numpy as np
import sys

gene2cpg = pd.read_csv('./data/gene2cpg_unique.lst', delimiter='\t', names=['gene','CpG'])
gene2cpg.set_index(['gene'], drop=True, inplace=True)

#lung = pd.read_csv('../../TCGA_biolinks_data/data/LIHC/LIHC_mvalues.tsv',delim_whitespace=True,nrows=20)
kidney = pd.read_csv('../../TCGA_biolinks_data/data/KIRC/KIRC_mvalues.tsv',delim_whitespace=True,nrows=200)

def get_probes(df_mval, gene_symbol) :
    list1 = gene2cpg.loc[gene_symbol]['CpG'].split(' ,')
    list2 = df_mval.columns.tolist()
    list3 = list(set(list1) & set(list2))
    out = df_mval[list3]
    return(out)

def binning(col):
    if np.dtype(col) == 'float' :
        mu=col.mean()
        sigma=col.std()
        bins1 = [-np.inf,(mu-0.7*sigma),mu,(mu+0.7*sigma),np.inf]
        bins2 = sorted(set(bins1)) # in case there are duplicate arguments in bin - if coloumn has all zeroes 
        group_names = np.arange(1,len(bins2))
        out = pd.cut(col,bins=bins2, labels=group_names)
    else: out=col
    return out

kidney_bin = kidney.apply(binning , axis=0)

gene_names = gene2cpg.index.tolist()
for gene_name in gene_names[4:]:
    delta_in = get_probes(kidney_bin,gene_name)
    outfile = './delta_input/'+gene_name+'.delta_in.csv'
    delta_in.to_csv(outfile, sep='\t', index=False, header=None)
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import scipy.stats as stats
from sklearn.preprocessing import Imputer
import scipy.special as sp
import sys

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

def clean_df(meth):
    #Note: meth coloums -> cpg islands 
    #Note: meth rows -> tcga cases
    all_col = meth.columns.tolist()
    cpg_col = [i for i in all_col if i[0:2]=='cg']
    meth = meth[cpg_col]
    
    # drop coloumns & rows which have nan in more than 50% instances
    meth = meth.dropna(axis=1, thresh = 0.5*meth.shape[0]) # for columns
    meth = meth.dropna(axis=0, thresh = 0.5*meth.shape[1]) # for rows
    meth_col = meth.columns.tolist()
    meth_index = meth.index.tolist()

    # Fill nan with the median values
    fill_NaN = Imputer(missing_values=np.nan, strategy='median', axis=0)
    meth =  pd.DataFrame(fill_NaN.fit_transform(meth), index=meth_index)
    meth.columns = meth_col
    return(meth)

def featureScaling_0to1(x):
    ub = 3
    lb = -3
    x_scale = (x - lb) / (ub - lb)
    x_scale = 0 if x_scale < 0 else x_scale
    x_scale = 1 if x_scale > 1 else x_scale
    return float(x_scale)

def featureScaling(x):
    if x < -6 :
        x = -6
    elif x > 6 :
        x = 6
    else:
        x = x
    return float(x)


def m_val(meth):
    meth_mval = meth.replace(0, 0.0001)
    meth_mval = meth_mval.applymap(lambda x : sp.logit(x))
    meth_mval = meth_mval.applymap(lambda x: featureScaling(x))
    meth_mval = meth_mval.applymap(lambda x: round(x,3))
    return (meth_mval)

def binning(col):
    mu=col.mean()
    sigma=col.std()
    bins1 = [-np.inf,(mu-0.7*sigma),mu,(mu+0.7*sigma),np.inf]
    bins2 = sorted(set(bins1)) # in case there are duplicate arguments in bin - if coloumn has all zeroes 
    group_names = np.arange(1,len(bins2))
    return pd.cut(col,bins=bins2, labels=group_names)

def histogram(df, col):
    plt.hist(df[col].values)
    plt.show()
    
def patient_id(tcga_id):
    patient_id = tcga_id[:12]
    return(patient_id)

def tissue_class(meth , tissue_location):
    meth.insert(0,'tissue_class',np.nan)
    #meth['tissue_class'] = np.nan
    all_samples = meth.index.tolist()
    normal_rows = [i for i in all_samples if i[13:14]=='1']
    tumor_rows = [i for i in all_samples if i[13:14]=='0']
    meth.loc[tumor_rows,'tissue_class'] = tissue_location + "_T"
    meth.loc[normal_rows,'tissue_class'] = tissue_location + "_N"
    return(meth)



in1 = pd.read_csv(sys.argv[1], delim_whitespace=True)
in1 = in1.T
in2 = clean_df(in1)
in3 = m_val(in2)
in4 = tissue_class(in3, sys.argv[2])

outfile = sys.argv[2] + '_mvalues.tsv'

in4.to_csv(outfile , index=True, sep='\t')



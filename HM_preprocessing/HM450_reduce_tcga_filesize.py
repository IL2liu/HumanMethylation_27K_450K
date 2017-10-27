import pandas as pd
import numpy as np
import sys 


def downsize(x):
    return(round(x,3))

infile = sys.argv[1] 
in1 = pd.read_csv(infile , delim_whitespace=True)
in2 = in1.loc[:,in1.dtypes == float]
in3 = in2.applymap(downsize)

outfile = infile.split('.txt')[0]
outfile = outfile+".round.txt"
in3.to_csv(outfile, index=True, sep='\t')




            
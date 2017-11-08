import pandas as pd
import numpy as np

fin = open("data/gene_list.txt")
gene_table = pd.read_csv("data/gene_id_list.csv",sep=",",header=0,index_col=0)
gene_list = []
for n in fin:
    gene_list.append(n.strip())

n = gene_table.loc["ENSG00000106006"]["entrez"]
n = gene_table.loc["ENSG00000106006"]["symbol"]
print n
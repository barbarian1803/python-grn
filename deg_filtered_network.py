from src.Network import Network
from sets import Set
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

def calculatePearson(expr_table,gene1,gene2):
    if gene1 == gene2:
        return False

    try:
        coeff = np.corrcoef(expr_table.loc[gene1],expr_table.loc[gene2])
        return coeff[0,1]
    except:
        print("Node data are not found, edge:"+gene1+"-"+gene2)
        return False

def calculateEdgesCorrCoeff(network):
    pearsonList = []

    for n in network.getNodes():
        node = network.getNode(n)
        if node.nodeType == "gene":
            continue

        for outWard in node.getOutWardNode():
            for inWard in node.getInWardNode():
                coeff = calculatePearson(gene_expr, outWard, inWard)
                if coeff:
                    pearsonList.append(abs(coeff))

    return pearsonList

def createHistogram(list,min,max,step):
    bins = np.arange(min, max, step)
    plt.xlim(min, max)
    plt.hist(list, bins)

    plt.show()


gene_expr = pd.read_csv("data/gene_expression.csv","\t",None,0,None,0)


# read network all data
# network =  Network("data/combine_network_filtered.csv","\t")
# coeff = calculateEdgesCorrCoeff(network)
# createHistogram(coeff,0,1,0.025)

#read network filtered with logFC and DEG genes data
network2 =  Network("data/network_filter_logfc_significant.csv","\t")
coeff = calculateEdgesCorrCoeff(network2)
createHistogram(coeff,0,1,0.025)

print(stats.describe(coeff))

maxOutward = -999
nodeName = ""
nodeDegree = {}

for n in network2.getNodes():
    node = network2.getNode(n)
    numOutward = len(node.getOutWardNode())
    nodeDegree[n] = numOutward
    if numOutward>maxOutward:
        nodeName = n
        maxOutward = numOutward

print str(maxOutward)+" "+nodeName

for key, value in sorted(nodeDegree.iteritems(), key=lambda (k,v): (v,k)):
    print "%s: %s" % (key, value)


node = network2.getNode("LXRBRXRA_Q5_01")
for x in node.getOutWardNode():
    print x
import sys
import re
from sets import Set

from src.GeneNode import GeneNode
from src.NetworkUtil import NetworkUtil

import os
fins = os.listdir("lfc_longest_chain")

def generateHTML(input_file,output_file):
    ##read network data
    network =  NetworkUtil.buildNetwork(input_file,"\t")

    ##read logFC data
    logFCData = {}
    logFC_file_name = "data/log2FC_complete_significant.csv"
    logFC_file = file(logFC_file_name)
    fc_values = {}
    header = {}
    for line in logFC_file:
        #skip header line
        if "Gene" in line:
            header = re.split(',',line)[1:]
            continue

        values = re.split(',',line)
        fc_values[values[0]] = values[1:]

    ##assign number to each node
    stack = list(NetworkUtil.getStartNodesOfNetwork(network))

    while len(stack)>0:
        node = stack.pop()
        if network[node].isStartNode():
            network[node].score = 1
        else:
            inward = next(iter(network[node].getInWardNode()))
            network[node].score = network[inward].score + 1
        stack= stack + list(network[node].getOutWardNode())

    listOutput = [""]*len(network)
    listNodes = network.values()
    listNodes.sort(key=lambda x: x.score, reverse=False)

    html = file("test.html").read()
    fout = file(output_file,"w")
    strHeader = ""
    for n in header:
        strHeader = strHeader + "<th>"+n+"</th>"

    html = html.replace("{header}",strHeader)

    strOut = ""
    for n in listNodes:
        if n.nodeType=="gene":
            strOut = strOut + "<tr>"+"<td>"+str(((n.score+1)/2))+"</td>"
            strOut = strOut + "<td>"+n.geneID+"</td>"
            for fc_val in fc_values[n.geneID]:
                val = float(fc_val)
                str_val = fc_val
                css = ""
                if val >= 1:
                    css = "pos_1"
                elif val >= 0.1:
                    css = "pos_"+str_val[0:3].replace(".","_")
                elif val >= 0.05:
                    css = "pos_0_05"
                elif val <= -1:
                    css = "neg_1"
                elif val <= -0.1:
                    css = "neg_"+str_val[1:4].replace(".","_")
                elif val <= -0.05:
                    css = "neg_0_05"

                strOut = strOut + "<td class='"+css+"'>"+str("%.2f"%val)+"</td>"
            strOut = strOut + "</tr>"
    html = html.replace("{content}",strOut)
    fout.write(html)



for fin in fins:
    fin_path = "lfc_longest_chain/"+fin
    fout = "lfc_longest_chain_result/"+fin.replace(".csv",".html")
    generateHTML(fin_path,fout)

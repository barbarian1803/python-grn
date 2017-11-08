from src.GeneNode import GeneNode
from src.Network import Network
from sets import Set

##read network data
network =  Network("data/combine_network_filtered.csv","\t")

##read logFC data
logFCData = {}
logFC_file_name = "data/log2FC_complete_significant.csv"
logFC_file = file(logFC_file_name)

for line in logFC_file:
    #skip header line
    if "Gene" in line:
        continue

    values = line.split(",")

    nodeName = values[0]
    nodeLogFC = float(values[1])

    if nodeName in network.getNodes():
        network.getNode(nodeName).logFC = nodeLogFC

## remove nonLOGFC nodes

network.removeNonDiffNode(0.65)

##Filter network from short and idependent chain
network.filterNetwork()

network.printNetworkTofile(Set(),file("data/network_filter_logfc_significant.csv","w"))

##get starting nodes
start_node = Set()
for n in network.getNodes():
    if network.getNode(n).isStartNode():
        start_node.add(n)

print "No of start node: "+str(len(start_node))

##get end nodes
end_node = Set()
for n in network.getNodes():
    if network.getNode(n).isEndNode():
        end_node.add(n)

print "No of end node: "+str(len(end_node))

print network.longestNodesChain()

for n in network.getStartNodesOfNetwork():
    network.printLongestChain(n,file("lfc_longest_chain/node_."+n+".csv","w"))

exit()

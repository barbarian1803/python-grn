from src.GeneNode import GeneNode
from src.NetworkUtil import NetworkUtil
from sets import Set

##read network data
network =  NetworkUtil.buildNetwork("data/combine_network_filtered.csv","\t")

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

    if nodeName in network:
        network[nodeName].logFC = nodeLogFC

## remove nonLOGFC nodes

network = NetworkUtil.removeNonDiffNode(network,0.65)


##Filter network from short and idependent chain
network = NetworkUtil.filterNetwork(network)

NetworkUtil.printNetworkTofile(network,Set(),file("data/network_filter_logfc_significant.csv","w"))

##get starting nodes
start_node = Set()
for n in network:
    if network[n].isStartNode():
        start_node.add(n)

print "No of start node: "+str(len(start_node))

##get end nodes
end_node = Set()
for n in network:
    if network[n].isEndNode():
        end_node.add(n)

print "No of end node: "+str(len(end_node))

print NetworkUtil.longestNodesChain(network)

for n in NetworkUtil.getStartNodesOfNetwork(network):
    NetworkUtil.printLongestChain(network,n,file("lfc_longest_chain/node_."+n+".csv","w"))

exit()

from GeneNode import GeneNode
from sets import Set

class Network:
    network = {}

    ## build network data structure using edge list csv file
    def __init__(self,fileName,split):
        edge_list = open(fileName,"r")

        for l in edge_list:
            #skip header line
            if "Source" in l:
                continue

            array = l.strip().split(split)
            self.createNode(array)


    ## helper method to add node to the network from edge list element
    def createNode(self,element):
        source = element[0]
        target = element[1]

        if source not in self.network:
            self.network[source] = GeneNode(source)

        if target not in self.network:
            self.network[target] = GeneNode(target)

        self.network[source].addOutWardNode(target)
        self.network[target].addInWardNode(source)

    def getEndNodesOfNetwork(self):
        endNodeSet = Set()
        for node in self.network:
            if self.network[node].isEndNode():
                endNodeSet.add(node)
        return endNodeSet

    def getStartNodesOfNetwork(self):
        startNodeSet = Set()
        for node in self.network:
            if self.network[node].isStartNode():
                startNodeSet.add(node)
        return startNodeSet

    @staticmethod
    def load_fimo_network_result(fname):
        checked_string = []
        table = []
        fobj = open(fname)
        for line in fobj:
            line = line.rstrip()
            if line[0] == "#":
                continue
            if line in checked_string:
                continue
            checked_string.append(line)
            split = line.split("\t")
            table.append({"source":split[1],"target":split[0],"type":"regulate"})
        return table

    @staticmethod
    def load_transfac_network(fname):
        table = []
        fobj = open(fname)
        for line in fobj:
            line = line.rstrip()
            if line[0] == "#":
                continue
            split = line.split("\t")
            to = split[0].replace("$","_")
            tfs = split[3].split("|")
            for tf in tfs:
                table.append({"source":tf,"target":to,"type":"protein"})

        return table

    @staticmethod   ## input file is fimo network and transfac gene-protein network and create merged network file
    def network_builder(fimo_network, transfac_network, output):
        fimo = Network.load_fimo_network_result(fimo_network)
        transfac = Network.load_transfac_network(transfac_network)

        network = fimo+transfac

        fobj = open(output,"w")
        fobj.write("Source\tTarget\tType\n")

        for n in network:
            fobj.write(n["source"].replace("V_","")+"\t"+n["target"].replace("V_","")+"\t"+n["type"]+"\n")


    ##recursive method that is called for longest chain calculation
    def calculateLongestChain(self,node,level,visited):

        if self.network[node].isEndNode() or node in visited:
            return level+1

        visited.add(node)
        maxVal = -999
        for n in self.network[node].getOutWardNode():
            result = self.calculateLongestChain(n,level+1,visited)
            if result>maxVal:
                maxVal = result

        return maxVal

    ## calculate longest chain from every starting nodes in the network
    def longestNodesChain(self):
        startNodes = self.getStartNodesOfNetwork()
        output = {}
        for n in startNodes:
            visited = Set()
            level = self.calculateLongestChain(n,0,visited)
            output[n] = level

        return output

    ##print network into edge list file, you can skip node if you want to exclude some node
    def printNetworkTofile(self,skipNode,fName):
        fName.write("Source\tTarget\tType\n")
        toPrint = Set()
        for n in self.network:
            self.network[n].printEdgeRel(toPrint)

        for n in toPrint:
            split = n.split("\t")
            if split[0] in skipNode or split[1] in skipNode:
                continue
            fName.write(n)

    @staticmethod   ## load node value from inputted file
    def loadNodeValue(fname):
        nodeValue = {}
        for line in fname:
            if "Gene" in line:
                continue

            array = line.split("\t")

            nodeValue[array[0]] = array[1]

        return nodeValue

    ## set node value and store it in the score attribute
    def setNodeValue(self,valueList):
        for n in self.network:
            if n in valueList:
                self.network[n].score = float(valueList[n])
            else:
                self.network[n].score = float(0)

    def removeNode(self,node):   ##remove node from the network
        for n in self.network[node].getInWardNode():
            self.network[n].removeOutWardNode(node)

        for n in self.network[node].getOutWardNode():
            self.network[n].removeInWardNode(node)

        self.network.pop(node)

    ##remove edge that connect exactly same source and target
    def removeSamePath(self):
        for n in self.network:
            if self.network[n].nodeType=="protein":
                continue

            connectedNodes = Set()
            outWardEdge = self.network[n].getOutWardNode()

            for edge in outWardEdge:

                if len(self.network[edge].getInWardNode())>1:
                    continue

                edge_outward = Set(self.network[edge].getOutWardNode())
                intersect = connectedNodes.intersection(edge_outward)
                for i in intersect:
                    self.network[edge].removeOutWardNode(i)
                    self.network[i].removeInWardNode(edge)
                connectedNodes = connectedNodes.union(edge_outward)

    ##remove protein nodes that are in the starting node or end node
    def filterProteinNode(self):
        start = self.getStartNodesOfNetwork()
        for n in start:
            if self.network[n].nodeType=="protein":
                self.removeNode(n)

        end = self.getEndNodesOfNetwork()
        for n in end:
            if self.network[n].nodeType=="protein":
                self.removeNode(n)

    ##remove short sub network that is isolated, remove protein node is start/end pos, remove same path
    def filterNetwork(self):
        output = self.longestNodesChain()
        while 1 in output.values() or 2 in output.values() or 3 in output.values():
            for n in output:
                if output[n]<=3:
                    self.removeNode(n)
            self.filterProteinNode()
            output = self.longestNodesChain()

        self.removeSamePath()

    def calculateDepth(self,node,score,visited):
        if self.network[node].score<score:
            self.network[node].score = score
        if node in visited:
            return

        visited.add(node)

        for n in self.network[node].getOutWardNode():
            if n in visited:
                continue
            self.calculateDepth(n,score+1,visited)

        return

    def printChain(self,node,toPrint):
        for n in self.network[node].getInWardNode():
            if self.network[n].score==self.network[node].score-1:
                toPrint.add(n+"\t"+node+"\n")
                self.printChain(n,toPrint)

    def printLongestChain(self,startnode,fout):
        for n in self.network:
            self.network[n].score=0

        visited = Set()
        self.calculateDepth(startnode,1,visited)

        maxScore = -1
        node = Set()
        for n in self.network:
            if maxScore < self.network[n].score:
                maxScore = self.network[n].score
                node = Set()
                node.add(n)
            elif maxScore==self.network[n].score:
                node.add(n)
        toPrint = Set()
        while node:
            self.printChain(node.pop(),toPrint)

        for n in toPrint:
            fout.write(n)

    def assignLogFC(self,logFC_file_name):
        logFCData = {}
        logFC_file = file(logFC_file_name)
        for line in logFC_file:
        #skip header line
            if "Gene" in line:
                continue
            values = line.split("\t")
            if values[0] in self.network:
                self.network[values[0]].logFC = float(values[1])

    def removeNonDiffNode(self,threshold):
        toRemove = Set()
        for n in self.network:
            if self.network[n].nodeType=="protein":
                continue

            if abs(self.network[n].logFC) <= threshold:
                toRemove.add(n)

        for n in toRemove:
            self.removeNode(n)

    def getEdgeType(self,nodeSrc):
        if self.network[nodeSrc].nodeType=="protein":
            return "regulate"
        else:
            return "protein"

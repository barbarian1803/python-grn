from sets import Set

class GeneNode:
    geneID = ""
    outWardEdge = Set()
    inWardEdge = Set()
    nodeType = ""
    logFC = ""
    score = ""

    def __init__(self,geneID):
        self.geneID = geneID
        self.outWardEdge = Set()
        self.inWardEdge = Set()
        
        self.logFC = 0
        self.score = 0
        
        if "ENSG" in geneID:
            self.nodeType = "gene"
        else:
            self.nodeType = "protein"

    def printNode(self):
        print "GeneID : "+self.geneID
        print "Type : "+self.nodeType
        print "Score : "+str(self.score)
        print "OutwardEdge:"
        print self.printEdge(self.outWardEdge)
        print "InwardEdge:"
        print self.printEdge(self.inWardEdge)

    def addInWardNode(self,node):
        self.inWardEdge.add(node)

    def addOutWardNode(self,node):
        self.outWardEdge.add(node)	

    def printEdge(self,edge):
        for n in edge:
            print n

    def printEdgeRel(self,toPrint):
        for n in self.outWardEdge:
            if "ENSG" in self.geneID:
                toPrint.add(self.geneID+"\t"+n+"\t"+"protein"+"\n")
            else:
                toPrint.add(self.geneID+"\t"+n+"\t"+"regulate"+"\n")
        for n in self.inWardEdge:
            if "ENSG" in self.geneID:
                toPrint.add(n+"\t"+self.geneID+"\t"+"regulate"+"\n")
            else:
                    toPrint.add(n+"\t"+self.geneID+"\t"+"protein"+"\n")
                    
    def isEndNode(self):
        outward = [x for x in self.outWardEdge if x is not None]
        return (len(outward)==0)
        
    def isStartNode(self):
        inward = [x for x in self.inWardEdge if x is not None]
        return (len(inward)==0)
    
    def getInWardNode(self):
        return [x for x in self.inWardEdge if x is not None]
    
    def getOutWardNode(self):
        return [x for x in self.outWardEdge if x is not None]
    
    def removeInWardNode(self,nodes):
        self.inWardEdge.remove(nodes)

    def removeOutWardNode(self,nodes):
        self.outWardEdge.remove(nodes)

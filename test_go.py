from src.NetworkUtil import NetworkUtil
from src.GeneNode import GeneNode
from src.obo_parser import *
from src.obo_util import *
from src.goaf_parser import *
from src.GeneralUtil import *

all_go = Parser.get_all_ontology("data/go-basic.obo")
all_goaf = GOAFParser("data/goa_human.gaf").parse()
gene_id_list = GeneralUtil.readGeneIDDatabase()

def print_go(p,filter):
    tags = all_go[p].tags
    if (filter=="all" or filter==tags["namespace"][0]) :
        if "is_obsolete" in tags and tags["is_obsolete"][0]=="true":
            return
        print tags["id"][0]+" "+tags["namespace"][0]+" "+tags["name"][0]

def print_parent(parent,filter):
    for p in parent:
        print_go(p,filter)
        if "is_a" in all_go[p].tags:
            print_parent(all_go[p].tags["is_a"],filter)

def print_go_of_gene(gene,filter):
    uniprot_id = gene_id_list[gene]["uniprot"]
    print "Gene "+gene+" "+gene_id_list[gene]["uniprot"]+"-"+gene_id_list[gene]["symbol"]
    for goa in all_goaf[uniprot_id]:
        p = goa[4]
        print_go(p,filter)
        parents = all_go[p].tags["is_a"]
        print_parent(parents,filter)

# uniprot_id =  gene_id_list["ENSG00000212907"]["uniprot"]
# print uniprot_id

print_go_of_gene("ENSG00000129654","molecular_function")

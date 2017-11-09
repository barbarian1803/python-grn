class GeneralUtil:

    @staticmethod
    def readGeneIDDatabase():
        data = {}
        csv1 = open("data/gene_id_list.csv")
        for line in csv1:
            if "ensembl" in line:
                continue

            array=line.strip().split("\t")
            if array[0] not in data:
                data[array[0]] = {}
            data[array[0]] = {"ensembl":array[0],"symbol":array[1],"entrez":array[2],"uniprot":""}

        csv2 = open("data/gene_uniprot.csv")
        for line in csv2:
            if "ensembl" in line:
                continue
            array = line.strip().split("\t")
            if array[0] not in data:
                data[array[0]] = {}
            data[array[0]]["uniprot"]=array[1]
        return data
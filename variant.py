# Data structure for storing variant data

class Variant():

    def __init__(self, chrom, start, end, ref, alt, gene, type, af,
                  pathway=None):
        self.chrom = chrom
        self.start = start
        if end != None:
            self.end = end
        else:
            self.end = start
        self.ref = ref
        self.alt = alt
        self.gene = gene
        self.type = type
        self.af = af
        #self.pathway = pathway
        #self.gender = gender
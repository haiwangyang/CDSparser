#!/usr/bin/env python

"""
fetch data/iCDS.fa, and concatenate iCDS into complete CDS.

check if CDS has the correct stopcodon (TAA, TGA, TAG)
check if CDS is multiple of three (3N)
translate CDS into PEP
output data/CDS.fa
output data/PEP.fa
output reports of stopcodon
output reports of 3N
"""

from pyfaidx import Fasta
from Bio.Seq import Seq

def write_dct_fasta(dct, path):
    with open(path, "w") as f:
        for i in dct.keys():
            f.write(">" + i + "\n" + dct[i] + "\n")

def write_dct_table(dct, path):
    with open(path, "w") as f:
        for i in sorted(dct.keys()):
            f.write(str(i) + "\t" + str(dct[i]) + "\n")

class CDS:
    """ CDS object """
    def __init__(self, strain):
        self.strain = strain
        self.iCDSfasta = Fasta("data/" + self.strain + ".iCDS.fa", as_raw = True)
        self.get_iCDS_dct()
        self.get_CDS_dct()
        self.check_start_codon()
        self.check_stop_codon()
        self.check_3n()
        self.get_PEP_dct()
        self.get_Mpos()
        #self.get_MPEP_dct()

    def get_iCDS_dct(self):
        """ get iCDS dict """
        dct = dict()
        for name in list(self.iCDSfasta.keys()):
            position, strand, id, parent, non = name.split(";")
            scaffold, range = position.split(":")
            start, end = range.split("-")
            start = int(start)
            end = int(end)
            geneid = id.replace("ID=", "").replace("-RA:cds", "")
            seq = str(self.iCDSfasta[name])
            if not geneid in dct.keys():
                dct[geneid] = dict()
            if not strand in dct[geneid].keys():
                dct[geneid][strand] = dict()
            dct[geneid][strand][start] = seq
        self.iCDSdct = dct

    def get_CDS_dct(self):
        dct = dict()
        for geneid in self.iCDSdct.keys():
            cds = ""
            for strand in self.iCDSdct[geneid].keys():
                if strand == "+":
                    starts = sorted(self.iCDSdct[geneid][strand])
                elif strand == "-":
                    starts = sorted(self.iCDSdct[geneid][strand], reverse=True)
                for start in starts:
                    cds = cds + self.iCDSdct[geneid][strand][start]
            dct[geneid] = cds
        self.CDSdct = dct

    def check_start_codon(self):
        dct = dict()
        for geneid in self.CDSdct.keys():
            cds = self.CDSdct[geneid]
            if cds.startswith("ATG"):
                dct[geneid] = "y"
            else:
                dct[geneid] = "n"
        self.ifstartcodon = dct

    def check_stop_codon(self):
        dct = dict()
        for geneid in self.CDSdct.keys():
            cds = self.CDSdct[geneid]
            if cds.endswith("TGA") or cds.endswith("TAA") or cds.endswith("TAG"):
                dct[geneid] = "y"
            else:
                dct[geneid] = "n" 
        self.ifstopcodon = dct

    def check_3n(self):
        dct = dict()
        for geneid in self.CDSdct.keys():
            cds = self.CDSdct[geneid]
            if len(cds) % 3 == 0:
                dct[geneid] = "y"
            else:
                dct[geneid] = "n"
        self.if3n = dct

    def get_PEP_dct(self):
        dct = dict()
        for geneid in self.CDSdct.keys():
            cds = self.CDSdct[geneid]
            pep = str(Seq(cds).translate())
            #print(cds)
            #print(pep)
            #print("\n")
            dct[geneid] = pep
        self.PEPdct = dct

    def get_Mpos(self):
        """
           get M position
        """
        dct = dict()
        for geneid in sorted(self.ifstartcodon.keys()):
            PEP = self.PEPdct[geneid]
            f = PEP.find("M")
            dct[geneid] = f
        self.Mpos = dct


    def get_MPEP_dct(self):
        """
        if PEP not start with M
        find the first M
        """
        dct = dict()
        for geneid in sorted(self.ifstartcodon.keys()):            
            if self.ifstartcodon[geneid] == "y":
                dct[geneid] = self.PEPdct[geneid]
            else:
                PEP = self.PEPdct[geneid]
                f = PEP.find("M")
                if f == -1:
                    dct[geneid + " ATTENTION " + str(f)] = PEP
                else:
                    dct[geneid + " ATTENTION " + str(f)] = PEP[f:len(PEP)]
        self.MPEPdct = dct       

if __name__ == '__main__':
    for strain in ['UCSC1', 'UMSG1', 'UMSG2', 'UMSG3']:
        print(strain)
        c = CDS(strain)
        write_dct_fasta(c.CDSdct, "output/" + strain + ".CDS.fa")
        write_dct_fasta(c.PEPdct, "output/" + strain + ".PEP.fa")
        write_dct_table(c.ifstopcodon, "output/" + strain + ".ifstopcodon.txt")
        write_dct_table(c.if3n, "output/" + strain + ".if3n.txt")
        write_dct_table(c.Mpos, "output/" + strain + ".Mpos.txt")

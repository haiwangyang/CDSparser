#!/usr/bin/env python

"""
fetch data/iCDS.fa, and concatenate iCDS into complete CDS.

if neccessary:

output data/CDS.fa

translate CDS

output data/PEP.fa
"""

from pyfaidx import Fasta
from Bio.Seq import Seq

def write_dct_fasta(dct, path):
    with open(path, "w") as f:
        for i in dct.keys():
            f.write(">" + i + "\n" + dct[i] + "\n")

def write_dct_table(dct, path):
    with open(path, "w") as f:
        for i in dct.keys():
            f.write(i + "\t" + dct[i] + "\n")

class CDS:
    """ CDS object """
    def __init__(self):
        self.iCDSfasta = Fasta("data/iCDS.fa", as_raw = True)
        self.get_iCDS_dct()
        self.get_CDS_dct()
        self.check_stop_codon()
        self.check_3n()
        self.get_PEP_dct()

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
            dct[geneid] = pep
        self.PEPdct = dct

if __name__ == '__main__':
    c = CDS()
    write_dct_fasta(c.CDSdct, "data/CDS.fa")
    write_dct_fasta(c.PEPdct, "data/PEP.fa")
    write_dct_table(c.ifstopcodon, "data/ifstopcodon.txt")
    write_dct_table(c.if3n, "data/if3n.txt")

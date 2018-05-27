#!/usr/bin/env python

"""
fetch data/iCDS.fa, and concatenate iCDS into complete CDS.

if neccessary:

output data/CDS.fa

translate CDS

output data/PEP.fa
"""

from pyfaidx import Fasta

class CDS:
    """ CDS object """
    def __init__(self):
        self.iCDSfasta = Fasta("data/iCDS.fa", as_raw = True)
        self.get_iCDS_dct()
        self.get_CDS_dct()
        #self.check_stop_codon()
        self.check_3n()

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
        for geneid in self.CDSdct.keys():
            cds = self.CDSdct[geneid]
            if cds.endswith("TGA") or cds.endswith("TAA") or cds.endswith("TAG"):
                print(geneid + "\t" + "stop" + "yes")
            else:
                print(geneid + "\t" + "stop" + "no")

    def check_3n(self):
        for geneid in self.CDSdct.keys():
            cds = self.CDSdct[geneid]
            if len(cds) % 3 == 0:
                print(geneid + "\t" + "3n" + "yes")
            else:
                print(geneid + "\t" + "3n" + "no")

if __name__ == '__main__':
    c = CDS()

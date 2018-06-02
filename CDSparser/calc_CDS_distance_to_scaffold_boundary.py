#!/usr/bin/env python

"""
Two jobs for this script
(1) calculate the distance to the scaffold boundary for both start and stop codon
(2) if start codon is not M, extend to find the nearest
"""

from __future__ import print_function, division
from concatenate_iCDS import CDS
from pyfaidx import Fasta
from Bio.Seq import Seq

def get_id_from_name(name):
    """ get gene id from bed name """
    return(name.split("ID=")[1].split("-")[0])

def get_start_end_from_string(start_end):
    """ 111.222 => 111, 222 """
    start, end = start_end.split(".")
    start = int(start)
    end = int(end)
    return(start, end)

def write_dct_table(dct, path):
    with open(path, "w") as f:
        for i in sorted(dct.keys()):
            f.write(str(i) + "\t" + str(dct[i]) + "\n")

class Bed:
    """ Bed object (now only use for maker) """
    def __init__(self, strain, feature):
        self.strain = strain
        self.feature = feature
        self.get_feature_ids()
        self.get_scaffold2len()
        self.get_id2left_right()
        self.get_id2distance()
        self.scaffoldfasta = Fasta("data/" + self.strain + "_CLC_de_novo_rmhost_mod.fa", as_raw = True)

    def get_feature_ids(self):
        """ get ids set
            get id to start_ends
            get id to scaffold
            get id to strand
        """
        self.ids = set()           # gene ids
        self.id2ses = dict()       # gene start_ends
        self.id2strand = dict()    # gene strand
        self.id2scaffold = dict()
        with open("data/" + self.strain + ".A." + self.feature + ".bed", "r") as f:
            for line in f.readlines():
                scaffold, start, end, name, one, strand = line.rstrip().split("\t")
                id = get_id_from_name(name)
                self.ids.add(id)
                self.id2scaffold[id] = scaffold
                self.id2strand[id] = strand
                if id not in self.id2ses.keys():
                    self.id2ses[id] = list()
                self.id2ses[id].append(start + "." + end)

    def get_scaffold2len(self):
        self.scaffold2len = dict()
        with open("data/" + self.strain + "_CLC_de_novo_rmhost_mod.fa.fai", "r") as f:
            for line in f.readlines():
                elements = line.rstrip().split("\t")
                self.scaffold2len[elements[0]] = int(elements[1])

    def get_id2left_right(self):
        """
           get the min and max position of gene
        """
        self.id2left = dict()      # min start
        self.id2right = dict()     # max end
        self.id2span = dict()      # CDS span
        for id in self.ids:
            self.id2left[id] = 9999999999999999
            self.id2right[id] = 0
            print(self.id2left[id])
            print(self.id2right[id])        

        for id in self.ids:
            for start_end in self.id2ses[id]:
                start, end = get_start_end_from_string(start_end)
                print(self.id2scaffold[id], id, self.id2strand[id], str(start), str(end))
                print("now left is " + str(self.id2left[id]) + "; now right is " + str(self.id2right[id]))
                if start < self.id2left[id]:
                    self.id2left[id] = start
                    print("assign left: " + str(start))
                if end > self.id2right[id]:
                    self.id2right[id] = end
                    print("assign right: " + str(end))
            print("\n\n")

    def get_id2distance(self):
        """
           get distance to scaffold boudnary
        """
        self.id2five_distance = dict()
        self.id2three_distance = dict()
        for id in self.ids:
            left = self.id2left[id]
            right = self.id2right[id]
            scaffold = self.id2scaffold[id]
            strand = self.id2strand[id]
            scaffoldlen = self.scaffold2len[scaffold]
            left_distance = left
            right_distance = scaffoldlen - right
            print(id + "\t" + "ld:" + str(left_distance) + "; rd: " + str(right_distance))
            if strand == "+":
                self.id2five_distance[id] = left_distance
                self.id2three_distance[id] = right_distance
                print(id + " five: " + str(left_distance) + "; three: " + str(right_distance))
            elif strand == "-":
                self.id2five_distance[id] = right_distance
                self.id2three_distance[id] = left_distance
                print(id + " five: " + str(right_distance) + "; three: " + str(left_distance))
            print("\n\n")

if __name__ == '__main__':
    for strain in ['UCSC1', 'UMSG1', 'UMSG2', 'UMSG3']:
        b = Bed(strain, "maker")
        write_dct_table(b.id2left, "output/" + strain + ".left.txt")
        write_dct_table(b.id2right, "output/" + strain + ".right.txt")
        write_dct_table(b.id2five_distance, "output/" + strain + ".five_distance.txt")
        write_dct_table(b.id2three_distance, "output/" + strain + ".three_distance.txt")


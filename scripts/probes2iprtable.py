#!/usr/bin/python

'''
Program accepting SNP name,contig,bp in tsv format and expanding these
features by a desired number of bp to make a gff file that can be intersected
with other genomic features.
'''

#-----------------------------------------------------
# Step 0
# Define classes
#-----------------------------------------------------


class Snp_obj(object):
    def __init__(self):
        """Return a Annot_obj whose name is *transcript_id*"""
        self.snp_id = ''
        self.snp_line = ''
        self.gene_list = []
        self.gene_dict = {}
    def print_obj(self):
        for gene in self.gene_list:
            outline = "\t".join([self.snp_id, gene, ";".join(self.gene_dict[gene]), self.snp_line])
            print outline


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import numpy as np
from collections import defaultdict

ap = argparse.ArgumentParser()

# ap.add_argument('--probes',required=True,type=str,help='tsv probe file containing: name,contig,bp')
# ap.add_argument('--assembly',required=True,type=str,help='Assembly fasta')
ap.add_argument('--sig_qtl',required=True,type=str,help='table of significan qtl')
ap.add_argument('--gene_intersect',required=True,type=str,help='table of qtl and their intersected gene models')
ap.add_argument('--ipr',required=True,type=str,help='interrposcan tsv file')

conf = ap.parse_args()

with open(conf.sig_qtl) as f:
    qtl_lines = f.readlines()

with open(conf.gene_intersect) as f:
    intersect_lines = f.readlines()

with open(conf.ipr) as f:
    ipr_lines = f.readlines()



#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

qtl_dict = defaultdict(list)
qtl_list = []

for line in qtl_lines [1:]:
    line = line.rstrip()
    split_line = line.split("\t")
    qtl = split_line[4].replace('"', '').replace('.', '-')
    qtl_list.append(qtl)
    # print split_line
    snp_obj = Snp_obj()
    snp_obj.snp_id = qtl
    snp_obj.snp_line = line
    qtl_dict[qtl] = snp_obj
    # print snp_obj

ipr_dict = defaultdict(list)
for line in ipr_lines:
    line = line.rstrip()
    split_line = line.split('\t')
    # print split_line
    gene = split_line[0]
    ipr_feature = split_line[2].replace('"', "") + " (" + split_line[1] + ")"
    ipr_dict[gene].append(ipr_feature)
    # print ipr_feature

for line in intersect_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    col9_id = split_line[8].split(';')[0]
    qtl = col9_id.replace('ID=', '')
    col18_id = split_line[17].split(';')[0]
    gene = col18_id.replace('ID=', '')
    # print split_line
    # print gene

    qtl_dict[qtl].gene_list.append(gene)
    qtl_dict[qtl].gene_dict[gene] = ipr_dict[gene]
    # print qtl_dict[qtl].gene_dict[gene]

for qtl in qtl_list:
    qtl_dict[qtl].print_obj()

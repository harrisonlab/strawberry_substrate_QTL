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


def add_assembly(contig_dict, assembly_lines):
    """Add protein sequence data to each gene"""
    for line in assembly_lines:
        line = line.rstrip()
        if line.startswith('>'):
            contig_id = line.replace('>', '')
            contig_dict[contig_id] = 0
        else:
            # print line
            # len(line)
            # print len(line)
            contig_dict[contig_id] += len(line)
            # print contig_dict[contig_id]
            # quit()


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import numpy as np
from collections import defaultdict

ap = argparse.ArgumentParser()

ap.add_argument('--probes',required=True,type=str,help='tsv probe file containing: name,contig,bp')
ap.add_argument('--assembly',required=True,type=str,help='Assembly fasta')
ap.add_argument('--bp',required=True,type=str,help='Expand probe positions by +/- X bp')
conf = ap.parse_args()

with open(conf.probes) as f:
    probe_lines = f.readlines()

with open(conf.assembly) as f:
    assembly_lines = f.readlines()

bp = conf.bp

#-----------------------------------------------------
# Step 2
# Import assembly
#-----------------------------------------------------

contig_dict = defaultdict(int)
add_assembly(contig_dict, assembly_lines)

# len_dict = defaultdict(int)
# for contig_id in contig_dict.keys():
#     length = len(contig_dict[contig_id])
#     len_dict[contig_id] = length

# print contig_dict

#-----------------------------------------------------
# Step 3
# Build gff features from probe locations
#-----------------------------------------------------


for line in probe_lines:
    line = line.rstrip()
    # print line
    split_line = line.split()
    # print split_line
    name = split_line[0]
    contig = split_line[1]
    position = split_line[2]

    # position_start = str(int(position) - 1000)
    position_start = str(int(position) - int(bp))
    position_stop = str(int(position) + int(bp))
    if int(position_start) < 1:
        position_start = 1
    if int(position_start) > int(contig_dict[contig]):
        position_end = int(contig_dict[contig])

    gff_contig = contig
    gff_source = 'probes2gff'
    gff_feature = 'probe'
    gff_start = str(position_start)
    gff_stop = str(position_stop)
    gff_col6 = '.'
    gff_strand = '.'
    gff_col8 = '.'
    gff_info = 'ID=' + name + ';''Name=' + name
    gff_line = "\t".join([gff_contig, gff_source, gff_feature, gff_start, gff_stop, gff_col6, gff_strand, gff_col8, gff_info])
    print gff_line

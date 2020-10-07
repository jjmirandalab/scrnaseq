#!/usr/bin/env python

import argparse
from pathlib import Path

import HTSeq

# This script assigns the alignments from a bam file to the feature annotations in a gtf file. It then prints out a
# matrix of counts corresponding to the the features of viral genomes. Since the feature annotations in the viral
# genomes are highly overlapping the script calls the intersection or union of features for each molecule as its
# feature.

parser = argparse.ArgumentParser(description='This script will ouput feature calls for a set of cells')
parser.add_argument('bam_path', type=str, help='path to bam file')
parser.add_argument('viral_gtf_path', type=str, help='path to gtf file')
parser.add_argument('cell_barcodes_path', type=str, help='path to list of cell barcodes file')
parser.add_argument('take_feature_union', type=str, help='"true" or "false" - take the union of features for the call '
                                                         'string, intersection if false')

args = parser.parse_args()

# check arguments
if not Path(args.bam_path).is_file():
    raise RuntimeError(f'could not find bam file, {args.bam_path}')
if not Path(args.viral_gtf_path).is_file():
    raise RuntimeError(f'could not find gtf file, {args.gtf_path}')
if not Path(args.cell_barcodes_path).is_file():
    raise RuntimeError(f'could not find cell barcodes file, {args.cell_barcodes_path}')
if args.take_feature_union not in ['true', 'false']:
    raise RuntimeError(f'take feature union must be "true" or "false"')

# parse take feature union
take_feature_union = True if args.take_feature_union == "true" else False

# read in cell barcodes of predicted real cells
cell_barcodes = set()
with open(args.cell_barcodes_path) as inf:
    for line in inf:
        cell_barcodes.add(line.strip())

# read in gtf and create a Genomic Array of Sets for all exons we find
viral_gtf_path = HTSeq.GFF_Reader(args.viral_gtf_path)
exons = HTSeq.GenomicArrayOfSets('auto', stranded=True)

# get all contigs from the gtf file
viral_gtf_contigs = {f.iv.chrom for f in viral_gtf_path}

for feature in viral_gtf_path:
    if feature.type == 'exon':
        exons[feature.iv] += feature.attr['gene_id']  # add gene id to this feature's coordinates in the exons array

# get alignments by umi by cell barcode
alignments_by_umi_by_cell_barcode = dict()
umi_to_ignore_by_cell_barcode = dict()  # cell-umi barcodes that map to non-viral contigs
for almnt in HTSeq.BAM_Reader(args.bam_path):

    assert isinstance(almnt, HTSeq.SAM_Alignment)

    # ignore secondary alignments and unmapped reads
    if not almnt.aligned or almnt.not_primary_alignment:
        continue

    # ignore alignments with invalid cell barcode or umi
    tags_present = {kv[0] for kv in almnt.optional_fields}
    if 'CB' not in tags_present or 'UB' not in tags_present:
        continue

    cell_barcode = almnt.optional_field('CB').split('-')[0]  # get cell barcode

    # ignore cells not mapped to a predicted real cell
    if cell_barcode not in cell_barcodes:
        continue

    umi_barcode = almnt.optional_field('UB')  # get umi

    # if this isn't a mapping to a viral contig, record that we don't care about this cell-umi barcode and continue
    if almnt.iv.chrom not in viral_gtf_contigs:

        if cell_barcode not in umi_to_ignore_by_cell_barcode:
            umi_to_ignore_by_cell_barcode[cell_barcode] = set()

        umi_to_ignore_by_cell_barcode[cell_barcode].add(umi_barcode)

        continue

    # store alignment
    if cell_barcode not in alignments_by_umi_by_cell_barcode:
        alignments_by_umi_by_cell_barcode[cell_barcode] = dict()

    if umi_barcode not in alignments_by_umi_by_cell_barcode[cell_barcode]:
        alignments_by_umi_by_cell_barcode[cell_barcode][umi_barcode] = set()

    alignments_by_umi_by_cell_barcode[cell_barcode][umi_barcode].add(almnt)

# get all called features and features by umi by cell barcode
feature_count_by_feature_by_cell_barcode = dict()
all_called_features = set()
for cb in alignments_by_umi_by_cell_barcode:
    for ub in alignments_by_umi_by_cell_barcode[cb]:

        # assume cell bc + umi combos map to 1 molecule, ignore maps to non-viral contigs
        if cb in umi_to_ignore_by_cell_barcode and ub in umi_to_ignore_by_cell_barcode[cb]:
            continue

        features_sets = list()
        contigs = set()
        for almnt in alignments_by_umi_by_cell_barcode[cb][ub]:
            assert isinstance(almnt, HTSeq.SAM_Alignment)

            contigs.add(almnt.iv.chrom)

            features = set()
            for iv, val in exons[almnt.iv].steps():
                features.update(val)

            if len(features) > 0:
                features_sets.append(features)

        if len(features_sets) < 1:  # ignore non-feature mappings
            continue

        if len(contigs) != 1:  # assume cell bc + umi combos map to 1 molecule, ignore ambiguous mappings
            continue

        contig = contigs.pop()

        called = features_sets.pop()
        for fset in features_sets:

            if take_feature_union:
                called = called.union(fset)

            else:
                called = called.intersection(fset)  # taking intersection to avoid ambiguos mappings

        if len(called) > 0:

            call_string = list()
            call_string.extend(called)
            call_string.sort()
            call_string = contig + '_' + '_'.join(call_string)

            # store the called feature
            all_called_features.add(call_string)

            # store the feature for this cell bc + umi
            if cb not in feature_count_by_feature_by_cell_barcode:
                feature_count_by_feature_by_cell_barcode[cb] = dict()

            if call_string not in feature_count_by_feature_by_cell_barcode[cb]:
                feature_count_by_feature_by_cell_barcode[cb][call_string] = 0

            feature_count_by_feature_by_cell_barcode[cb][call_string] += 1

# print output
out_cell_barcodes = list()
out_cell_barcodes.extend(cell_barcodes)
out_features = list()
out_features.extend(all_called_features)
out_features.sort()
print('feature_id\t' + '\t'.join(out_cell_barcodes))
for feature in out_features:
    out_string = list()
    out_string.append(feature)
    for cb in out_cell_barcodes:
        if cb not in feature_count_by_feature_by_cell_barcode:
            out_string.append(str(0))
        else:
            if feature not in feature_count_by_feature_by_cell_barcode[cb]:
                out_string.append(str(0))
            else:
                out_string.append(str(feature_count_by_feature_by_cell_barcode[cb][feature]))
    print('\t'.join(out_string))

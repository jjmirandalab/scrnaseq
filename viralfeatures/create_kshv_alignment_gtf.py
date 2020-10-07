#!/usr/bin/env python

import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description='This script will cleanup a the U75698.1 gff annotation and output a gff '
                                             'for alignment')
parser.add_argument('gff_path', type=str, help='path to gtf file')

gff_path = Path(parser.parse_args().gff_path)

if not gff_path.is_file():
    raise RuntimeError(f'could not find gtf file, {gff_path}')


class Line:

    def __init__(self, lline):

        lline = lline.split('\t')

        self.fields = lline[:8]

        self.attributes = {k: v for k,v in [f.split('=') for f in lline[8].split(';')]}

    def get_type(self):

        return self.fields[2]

    def set_type(self, ttype):

        self.fields[2] = ttype


gene_ids = set()

with gff_path.open() as ifh:

    for line in ifh:

        line = line.strip()
        if line == "":
            continue

        if line[0] == "#":
            continue

        line = Line(line)

        gid = None
        output_attributes = list()

        if line.get_type() == 'exon':

            gid = line.attributes['Parent']

            output_attributes.append('gene_id "' + gid + '"')
            output_attributes.append('transcript_id "' + gid + '"')
            output_attributes.append('gene_name "' + line.attributes['product'].replace(' ', '-') + '"')

        elif line.get_type() == 'CDS':

            gid = line.attributes['Name']

            output_attributes.append('gene_id "' + gid + '"')
            output_attributes.append('transcript_id "' + gid + '"')
            output_attributes.append('gene_name "' + line.attributes['product'].replace(' ', '-') + '"')

            line.set_type('exon')

        else:
            continue

        if gid in gene_ids:
            RuntimeError(f'found a duplicate gene id')

        print('\t'.join(line.fields) + '\t' + '; '.join(output_attributes))



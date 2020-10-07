#!/usr/bin/env python

import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description='This script will cleanup the chrEBV_B_95_8_Raji.char_corrected.ann file '
                                             'for alignment')
parser.add_argument('ann_path', type=str, help='path to ann file')

ann_path = Path(parser.parse_args().ann_path)

if not ann_path.is_file():
    raise RuntimeError(f'could not find ann file, {ann_path}')


class Line:

    def __init__(self, lline):

        self.line = lline

        fields = lline.split('\t')

        self.seqname = fields[2]

        # parse attributes
        self.gene_id = fields[0]
        self.gene_name = fields[1]
        self.transcript_id = "t-" + fields[4] + ".." + fields[5]

        self.attributes = 'gene_id "' + self.gene_id + '"; transcript_id "' + self.transcript_id + \
                          '"; gene_name "' + self.gene_name + '"'

        # parse strand
        self.strand = fields[3]

        if self.strand not in {'1', '0'}:
            RuntimeError('invalid strand field, ' + self.strand )

        if self.strand == '1':
            self.strand = '+'
        else:
            self.strand = '-'

        # parse exons
        self.exon_intervals = list()

        starts = fields[9]
        stops = fields[10]

        while starts[-1] == ',':
            starts = starts[:-1]

        while stops[-1] == ',':
            stops = stops[:-1]

        starts = starts.split(',')
        stops = stops.split(',')

        if len(starts) != len(stops):
            RuntimeError('number of exon starts doesn\'t equal the number of stops')

        for i in range(len(starts)):
            self.exon_intervals.append((starts[i], stops[i]))

        self.gtf_lines = self.get_gtf_lines()
        self.first_start = int(self.gtf_lines[0].split('\t')[3])

    def get_gtf_lines(self):

        # gtf format
        # seqname source feature start end score strand frame attributes

        gtf_lines = list()
        for eiv in self.exon_intervals:

            gtf_lines.append('\t'.join([self.seqname, 'B_95_8_Raji', 'exon', eiv[0], eiv[1], '.', self.strand, '.',
                                        self.attributes]))

        return sorted(gtf_lines, key=lambda gtfline: int(gtfline.split('\t')[3]))


with ann_path.open() as ifh:

    lines_by_first_start = dict()
    for line in ifh:

        line = line.strip()
        if line == "":
            continue

        if line[0] == "#":
            continue

        line = Line(line)

        if line.first_start not in lines_by_first_start:
            lines_by_first_start[line.first_start] = list()

        lines_by_first_start[line.first_start].append(line.gtf_lines)

    for fs in sorted(lines_by_first_start.keys()):

        for lines in lines_by_first_start[fs]:
            for l in lines:
                print(l)







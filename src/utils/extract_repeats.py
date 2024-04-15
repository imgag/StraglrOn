#!/usr/bin/env python
#@author is readmanchiu (from straglr) with modifications by Thomas Braun
import argparse
import pysam
from collections import defaultdict
import re

def parse_tsv(tsv, loci=None):
    support = defaultdict(dict)
    with open(tsv, 'r') as ff:
        for line in ff:
            if line[0] == '#':
                continue
            cols = line.rstrip().split('\t')
            locus = cols[0] + ":" + cols[1] + "-" + cols[2]
            status = "does not exist"
            read_name = cols[5]
            size = cols[7]
            read_start = cols[8]
            strand = cols[9]
            
            if loci is not None and locus not in loci:
                continue
            support[locus][read_name] = int(read_start), int(size), strand
        
    return support

def extract_repeats(bam, support, flank_size=10):
    seqs = {}
    for locus in support:
        seqs[locus] = []
        chrom, start, end = re.split('[:-]', locus)
        for aln in bam.fetch(chrom, int(start), int(end)):
            
            if aln.query_name in support[locus] and not aln.query_name in seqs:
                rlen = aln.infer_read_length()
                if support[locus][aln.query_name][2] == '+':
                    start, end = support[locus][aln.query_name][0], support[locus][aln.query_name][0] + support[locus][aln.query_name][1]
                else:
                    start = rlen - (support[locus][aln.query_name][0] + support[locus][aln.query_name][1])
                    end = start + support[locus][aln.query_name][1]

                try:
                    repeat_seq = aln.query_sequence[start:end].lower()
                    left = max(0, start - flank_size), start
                    right = end, min(end + flank_size, rlen)
                    left_seq = aln.query_sequence[left[0]:left[1]].upper()
                    right_seq = aln.query_sequence[right[0]:right[1]].upper()
                    seq = left_seq + repeat_seq + right_seq
                    seqs[locus].append((aln.query_name, support[locus][aln.query_name][1], seq))
                except:
                    print('problem extracting repeat from {}'.format(aln.query_name))

    return seqs

def report(seqs, out_fa):
    with open(out_fa, 'w') as out:
        for locus in sorted(seqs.keys()):
            for read_name, repeat_size, seq in seqs[locus]:
                out.write('>{} {} {}\n{}\n'.format(read_name, locus, repeat_size, seq))

def parse_args():
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    return args
    
def fastaMaker(tsv, locus, bam, flank_size, repeat_fasta_path):
    
    support = parse_tsv(tsv, locus)
    bam = pysam.AlignmentFile(bam)
    seqs = extract_repeats(bam, support, flank_size)
    report(seqs, repeat_fasta_path)

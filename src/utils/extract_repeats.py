#!/usr/bin/env python
#@author is readmanchiu (from straglr) with modifications by Thomas Braun
import argparse
import pysam
from collections import defaultdict
import re
from .Structures import RepeatSequence, MethylationCall

def parse_tsv(tsv, loci=None):
    """Parse TSV file and return read support information"""
    support = defaultdict(dict)
    
    with open(tsv, 'r') as ff:
        for line in ff:
            if line[0] == '#':
                continue
                
            cols = line.rstrip().split('\t')
            locus = cols[0] + ":" + cols[1] + "-" + cols[2]
            status = cols[14].strip()  # was "does not exist"

            # ignore skipped/failed reads
            if status != "full":
                continue
            read_name = cols[7]  # was cols[5]
            size = cols[10]  # was cols[7]
            read_start = cols[11]  # was cols[8]
            strand = cols[12]  # was cols[9]
            
            if loci is not None and locus not in loci:
                continue
            # print(read_name)
            support[locus][read_name] = int(read_start), int(size), strand
    
    return support


def parse_bam(bam_file, support, flank_size=10):
    """
    Parse BAM file to extract sequences and methylation data
    
    Args:
        bam_file: Path to BAM file or pysam.AlignmentFile object
        support: Coordinates from supporting read tsv
        flank_size: Size of flanking regions to include
    """
    if isinstance(bam_file, str):
        bam = pysam.AlignmentFile(bam_file, "rb")
    else:
        bam = bam_file

    seqs = []
    
    for locus in support:
        chrom, start, end = re.split('[:-]', locus)
        for aln in bam.fetch(chrom, int(start), int(end)):

            if aln.query_name in support[locus]:
                rlen = aln.infer_read_length()
                if support[locus][aln.query_name][2] == '+':
                    start, end = support[locus][aln.query_name][0], support[locus][aln.query_name][0] + support[locus][aln.query_name][1]
                else:
                    start = int(rlen - (support[locus][aln.query_name][0] + support[locus][aln.query_name][1])) 
                    end = int(start + support[locus][aln.query_name][1])

                try:
                    repeat_seq = aln.query_sequence[start:end].lower()
                    left = max(0, start - flank_size), start
                    right = end, min(end + flank_size, rlen)
                    left_seq = aln.query_sequence[left[0]:left[1]].upper()
                    right_seq = aln.query_sequence[right[0]:right[1]].upper()
                    seq = left_seq + repeat_seq + right_seq

                    # Extract methylation data
                    methylation_calls = []
                    if aln.has_tag('MM') and aln.has_tag('ML'):
                        mod_string = aln.get_tag('MM')
                        mod_quals = aln.get_tag('ML')
                        
                        for mod in mod_string.split(';'):
                            if not mod or not mod.startswith('C+m'):
                                continue
                                
                            _, deltas = mod.split(',', 1)
                            deltas = list(map(int, deltas.split(',')))
                            pos = 0
                            for delta, qual in zip(deltas, mod_quals):
                                pos += delta + 1 # Delta: Number of bases to skip. Add 1 to get the position of the modified base
                                if int(left[0]) <= pos <= int(right[1]):  # Only include methylation calls in our region
                                    methylation_calls.append(
                                        MethylationCall(
                                            position=pos - start,  # Adjust position relative to sequence start
                                            is_methylated=qual >= 200,
                                            quality_score=qual
                                        )
                                    )

                    # Store results
                    seqs.append(
                        RepeatSequence(
                            locus= locus,
                            read_name= aln.query_name, 
                            repeat_size= end - start, 
                            sequence= seq,
                            start_position= aln.reference_start,
                            end_position= aln.reference_end,
                            left_flank= left_seq,
                            repeat_sequence= repeat_seq,
                            right_flank= right_seq,
                            methylation_calls= methylation_calls
                        )
                    )
                except:
                    print('problem extracting repeat from {}'.format(aln.query_name))

    if isinstance(bam_file, str):
        bam.close()
        
    return seqs

def write_fasta(seqs, out_fa):
    with open(out_fa, 'w') as out:
        for seq in seqs:
            out.write('>{} {} {}\n{}\n'.format(
                seq.read_name,
                seq.locus,
                seq.repeat_size,
                seq.sequence
            ))

def parse_args():
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    return args
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
            if status != "full" and status != "partial":
                continue
            read_name = cols[7]  # was cols[5]
            size = cols[10]  # was cols[7] (repeat length in bases)
            read_start = cols[11]  # was cols[8] (start position of repeat on read)
            strand = cols[12]  # was cols[9]
            
            if loci is not None and locus not in loci:
                continue
            # print(read_name)
            support[locus][read_name] = int(read_start), int(size), strand
    
    return support


def parse_modifications(aln, start, end, flank_size):
    """Parse modifications from alignment object into MethylationCall objects"""
    methylation_calls = {}
    
    # Process modified_bases if available
    if aln.modified_bases:
        # Iterate through all modifications
        for (base, strand, mod_type), positions in aln.modified_bases.items():
            for pos, qual in positions:
                # Check if position is within our region of interest
                if start - flank_size <= pos <= end + flank_size:
                    # Create or update MethylationCall object
                    if pos not in methylation_calls:
                        methylation_calls[pos] = MethylationCall(
                            position=pos - (start - flank_size)  # Adjust position relative to sequence start
                        )
                    # Add modification quality
                    methylation_calls[pos].modifications[mod_type] = qual
    
    return list(methylation_calls.values())

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
                    # start: read_start, end: read_start + repeat_length
                    start, end = support[locus][aln.query_name][0], support[locus][aln.query_name][0] + support[locus][aln.query_name][1]
                else:
                    # start: read_end - repeat_length, end: read_end
                    start = int(rlen - (support[locus][aln.query_name][0] + support[locus][aln.query_name][1])) 
                    end = int(start + support[locus][aln.query_name][1])

                # Check if aln.query_sequence is not None
                if aln.query_sequence is not None:
                    repeat_seq = aln.query_sequence[start:end].lower()
                    left = max(0, start - flank_size), start
                    right = end, min(end + flank_size, rlen)
                    left_seq = aln.query_sequence[left[0]:left[1]].upper()
                    right_seq = aln.query_sequence[right[0]:right[1]].upper()
                    seq = left_seq + repeat_seq + right_seq

                    # Extract methylation data
                    methylation_calls = parse_modifications(aln, start, end, flank_size)

                    # Store results
                    seqs.append(
                        RepeatSequence(
                            locus=locus,
                            read_name=aln.query_name, 
                            repeat_size=end - start, 
                            sequence=seq,
                            start_position=aln.reference_start,
                            end_position=aln.reference_end,
                            left_flank=left_seq,
                            repeat_sequence=repeat_seq,
                            right_flank=right_seq,
                            methylation_calls=methylation_calls
                        )
                    )
                else:
                    print(f"Warning: No query sequence for alignment {aln.query_name}. Skipping this alignment.")

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
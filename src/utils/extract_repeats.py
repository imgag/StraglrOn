#!/usr/bin/env python
#@author is readmanchiu (from straglr) with modifications by Thomas Braun
import argparse
import pysam
from collections import defaultdict
import re
from .Structures import RepeatSequence, MethylationCall

def parse_tsv(tsv, loci=None):
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


def parse_bam(bam_file, expansion, flank_size=10):
    """
    Parse BAM file to extract sequences and methylation data
    
    Args:
        bam_file: Path to BAM file or pysam.AlignmentFile object
        region: String in format "chr:start-end"
        flank_size: Size of flanking regions to include
    """
    if isinstance(bam_file, str):
        bam = pysam.AlignmentFile(bam_file, "rb")
    else:
        bam = bam_file

    region = f"{expansion.chr}:{expansion.start}-{expansion.end}"

  
    # Parse region
    chrom, pos = region.split(":")
    start, end = map(int, pos.split("-"))
    
    seqs = []
    
    for aln in bam.fetch(chrom, start, end):
        if aln.is_secondary or aln.is_supplementary:
            continue
            
        try:
            # Get sequence
            sequence = aln.query_sequence
            if not sequence:
                continue
                
            rlen = len(sequence)
            
            # Extract repeat region with flanks
            repeat_start = max(0, aln.reference_start - start)
            repeat_end = min(rlen, repeat_start + (end - start))
            
            repeat_seq = sequence[repeat_start:repeat_end].lower()
            
            # Extract flanking sequences
            left_start = max(0, repeat_start - flank_size)
            right_end = min(rlen, repeat_end + flank_size)
            
            left_seq = sequence[left_start:repeat_start].upper()
            right_seq = sequence[repeat_end:right_end].upper()
            
            full_seq = left_seq + repeat_seq + right_seq
            
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
                        if left_start <= pos <= right_end:  # Only include methylation calls in our region
                            methylation_calls.append(
                                MethylationCall(
                                    position=pos - left_start,  # Adjust position relative to sequence start
                                    is_methylated=qual >= 200,
                                    quality_score=qual
                                )
                            )

    
            # Store results
            seqs.append(
                RepeatSequence(
                    expansion= expansion,
                    read_name= aln.query_name, 
                    repeat_size= repeat_end - repeat_start, 
                    sequence= full_seq,
                    start_position= aln.reference_start,
                    end_position= aln.reference_end,
                    left_flank= left_seq,
                    repeat_sequence= repeat_seq,
                    right_flank= right_seq,
                    methylation_calls= methylation_calls
                )
            )
                
        except Exception as e:
            print(f'Problem extracting repeat from {aln.query_name}: {str(e)}')
            continue

    if isinstance(bam_file, str):
        bam.close()
        
    return seqs

def write_fasta(seqs, out_fa):
    with open(out_fa, 'w') as out:
        for seq in seqs:
            out.write('>{} {}:{}-{} {}\n{}\n'.format(
                seq.read_name,
                seq.expansion.chr,
                seq.expansion.start,
                seq.expansion.end,
                seq.repeat_size,
                seq.sequence
            ))

def parse_args():
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    return args
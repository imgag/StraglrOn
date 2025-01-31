import os
from collections import deque
from pathlib import Path
import argparse
from utils.Structures import Expansion
import utils.ReaderAnalyser as ra
import utils.Visualiser as vis
import utils.extract_repeats as extract_repeats
import shutil
import pysam


def is_valid_file(parser, arg, type):
    """
    Check if the provided file exists and has a .bed or .tsv extension.
    """
    if type == "bed":
        if not os.path.exists(arg):            parser.error(f"The file '{arg}' does not exist.")
        elif not arg.endswith('.bed'):
            parser.error(f"The file '{arg}' is not a .bed file.")
        else:
            return arg
        
    if type == "tsv":
        if not os.path.exists(arg):
            parser.error(f"The file '{arg}' does not exist.")
        elif not arg.endswith('.tsv'):
            parser.error(f"The file '{arg}' is not a .tsv file.")
        else:
            return arg
    

# Argument parser for command line interface including necessary file and options

parser = argparse.ArgumentParser()

# These arguments are always required

parser.add_argument("path_input_bed", type=lambda x: is_valid_file(parser, x, "bed"), help="Path to input bed file")
parser.add_argument("path_input_tsv", type=lambda x: is_valid_file(parser, x, "tsv"), help="Path to input tsv file")
parser.add_argument("loci_file", type=lambda x: is_valid_file(parser, x, "bed"), help="Path to Loci file used for straglr analysis")
parser.add_argument("-o", "--output", type=str, required=True, help="Path to output folder")

# These Arguments are optional and produce histograms for each locus and sort the output .txt file based on a normalized increase in repeats respectively

parser.add_argument("--hist", action="store_true", help="Plots histograms of pathogenic expansions")
parser.add_argument("--score", action="store_true", help="Expansion in output file is sorted by normalized size difference score")

# Activates new clustering method 

parser.add_argument("--altclust", action="store_true", help="Uses Thomas Clustering")
parser.add_argument("-c", "--cutoff", type=int, default=2, help="Sets number of reads cutoff when clustering unimodal or bimodal allele read frequencies")

# These arguments are required for producing the allele length visualization 

parser.add_argument("--alleles", action="store_true", help="Turns on the allele visualization")
parser.add_argument("--bam", type=str, help="Location of Bam file of interest")
parser.add_argument("--flank", type=int, default=25, help="flank size. Default:25")
parser.add_argument("--genome", type=str, help="location of reference genome")

args = parser.parse_args()

def main():

    # Read files        
    loci_dict = ra.lociBedReader(args.loci_file)
    expansions = ra.resultBedReader(args.path_input_bed, loci_dict)

    # Generate Histogram data and update Extensions
    vis.getHistData(args.path_input_tsv, expansions)
    
    # Alternative Clustering
    for expansion_object in expansions:
        if args.altclust:
            ra.newGenotyping(expansion_object, args.cutoff, False)
        ra.expansionScorer(expansion_object, args.altclust)

    # New scoring
    if args.score:
        #Sorts expansions by calculated normalized score
        expansions.sort(key=lambda x: float("-Inf") if x.norm_score is None else x.norm_score, reverse=True)
        
    else:
        #Sorts expansion by chromosome
        expansions.sort(key=lambda x: x.chr)
    
   
    if not os.path.exists(args.output):
        os.mkdir(args.output)
     
    # Generate histograms
    if args.hist:
        print("Generating read distribution histograms ...")
        for expansion in expansions:
            vis.plotHistogram(expansion, args.output, args.altclust)
        print("Histograms completed.")

    # Generate allele visualizations
    if args.alleles:
        print("Generating allele composition graphs ...")
        for expansion in expansions:
            # Get read support data
            support = extract_repeats.parse_tsv(
                args.path_input_tsv,
                expansion.chr + ":"+ expansion.start + "-" + expansion.end 
            )
            
            if not support:
                print(f"No reads found for {expansion.repeat_id}, skipping visualization")
                continue

            # Get sequences
            sequences = extract_repeats.parse_bam(
                args.bam,
                support,
                args.flank
            )
            
            if sequences:
                # Write FASTA and create visualization
                fasta_path = os.path.join(args.output, f"{expansion.title}.fa")
                extract_repeats.write_fasta(sequences, fasta_path)
                
                vis.alleleVisualiser(
                    expansion.repeat_unit,
                    args.flank,
                    expansion.title,
                    args.output,
                    expansion.chr,
                    expansion.start,
                    expansion.end,
                    args.genome,
                    sequences
                )

if __name__ == "__main__":
    main()
import os
from collections import deque
from pathlib import Path
import argparse
from utils.Structures import Expansion
import utils.ReaderAnalyser as ra
import utils.Visualiser as vis
import utils.extract_repeats as extract_repeats
import shutil


def is_valid_file(parser, arg, type):
    """
    Check if the provided file exists and has a .bed or .tsv extension.
    """
    if type == "bed":
        if not os.path.exists(arg):
            parser.error(f"The file '{arg}' does not exist.")
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


def outputWriter(output_expansion: "list[Expansion]", sample_id, args):
            
    # Base output here
    result_file_path = args.output + "/" + sample_id + ".tsv"
    if not (os.path.exists(args.output)):
        os.mkdir(args.output)
        
    with open(result_file_path, 'w') as f:
        print("Writing result file ...")
        f.write("#chr\tstart\tend\trepeat_id\trepeat_unit\tcopy_number\tsize\twt_size\tin_pathogenic_range\tsize_difference\tallele1_support\tallele2_support\tscore\n")
        
        if args.altclust:
            for x in output_expansion:
                f.write("{}\t{}\t{}\t{}\t{:.1f}/{:.1f}\t{}/{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(x.chr, x.start, x.end, x.repeat_id, (x.new_allele1/len(x.repeat_unit)), (x.new_allele2/len(x.repeat_unit)), x.new_allele1, x.new_allele2, x.wt_size, x.new_in_pathogenic_range, x.new_size_difference, x.new_allele1_support, x.new_allele2_support,x.new_norm_score))
        else:
            for x in output_expansion:
                f.write("{}\t{}\t{}\t{}\t{:.1f}/{:.1f}\t{}/{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(x.chr, x.start, x.end, x.repeat_id, (x.allele1_size/len(x.repeat_unit)), (x.allele2_size/len(x.repeat_unit)), x.allele1_size, x.allele2_size, x.wt_size, x.in_pathogenic_range, x.size_difference, x.allele1_support, x.allele2_support,x.norm_score))                

        print("Result file completed.")      
    
    # Base output plus Histograms
    if args.hist:
        print("Generating read distribution histograms ...")
        
        if args.altclust:
            for x in output_expansion:
                vis.plotHistogram(x, args.output, args.altclust)
        else:
            for x in output_expansion:
                vis.plotHistogram(x, args.output, args.altclust)
                
        print("Histograms completed.")

    
    # Base with Allele Visualization
    if args.alleles:
        print("Generating allele composition graphs ...")
        
        for x in output_expansion:    
            extract_repeats.fastaMaker(args.path_input_tsv, x.chr + ":"+ x.start + "-" + x.end, args.bam, args.flank, args.output + "/" + x._title + ".fa")
            vis.alleleVisualiser(args.output + "/" + x._title + ".fa", x.repeat_unit, args.flank, x._title, args.output, x.chr, x.start, x.end, args.genome)
        
        print("Allele composition graphs completed.")


def main():
        
    loci_dict = ra.lociBedReader(args.loci_file)
    expansions = ra.resultBedReader(args.path_input_bed, loci_dict)
    sample_id = Path(args.path_input_bed).stem

    vis.getHistData(args.path_input_tsv, expansions)
    
    for expansion_object in expansions:
        if __name__ == '__main__':
            expansion_object
        if args.altclust:
            ra.newGenotyping(expansion_object, args.cutoff, False)
            
        ra.expansionScorer(expansion_object, args.altclust)

    if args.score:
        #Sorts expansions by calculated normalized score
        expansions.sort(key=lambda x: float("-Inf") if x.norm_score is None else x.norm_score, reverse=True)
        
    else:
        #Sorts expansion by chromosome
        expansions.sort(key=lambda x: x.chr)
        
    outputWriter(expansions, sample_id, args)
    
if __name__ == "__main__":
    main()
from pathlib import Path
from .Structures import Expansion
import matplotlib.pyplot as plt
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
import re
import numpy as np
from Bio import SeqIO
import pysam

# set hashsalt to const value to make plots deterministic for tests
rc('svg', hashsalt="Totally_Random_String")

def getHistData(file, expansions: "list[Expansion]"):
    expansions_read_lists = {}
    for expansion in expansions:
        title=Path(file).stem + "_" + expansion.repeat_id
        genotype1=expansion.allele1_size
        genotype2=expansion.allele2_size
        if genotype1 == genotype2:
            with open(file) as readFileTSV:
                read_reader=csv.reader(readFileTSV, delimiter="\t")
                read_list = []
                read_coords_dict = {}
                for line in read_reader:
                    if not line[0].startswith('#'):
                        if expansion.start+expansion.repeat_unit == line[1]+line[3]:
                            # skip reads with read_status != full
                            if line[14] != "full":
                                continue
                            # read column 'read'/'read_name' and 'read_start'
                            read_coords_dict.update({line[7]: int(line[11])})
                            if line[11] == "NA":
                                read_coords_dict.update({line[7]: np.nan})
                            # read column 'size'
                            read_list.append(int(line[10]))

                expansion.read_dict = read_coords_dict   
                expansion.read_list = [read_list] 
                expansion.title = title 
                expansions_read_lists.update({expansion.repeat_id: (genotype1,genotype2,title,read_list)})
        else:
            with open(file) as readFileTSV:
                read_reader=csv.reader(readFileTSV, delimiter="\t")
                read_list_A1 = []
                read_list_A2 = []
                read_coords_dict = {}
                #print(locus[-2], locus[-1])
                for line in read_reader:
                    if not line[0].startswith('#'):
                        if expansion.start+expansion.repeat_unit == line[1]+line[3]:
                            # skip reads with read_status != full
                            if line[14] != "full":
                                continue
                            if expansion.copy_numberA1 == line[13]:  # read column 'allele' (not 'copy_number')
                                read_list_A1.append(int(line[10]))  # read column 'size'
                            if expansion.copy_numberA2 == line[13]:  # read column 'allele' (not 'copy_number')
                                read_list_A2.append(int(line[10]))  # read column 'size'
                            # read column 'read'/'read_name' and 'read_start'
                            read_coords_dict.update({line[7]:int(line[11])})
                expansion.read_dict = read_coords_dict   
                expansion.read_list = [read_list_A1,read_list_A2] 
                expansion.title = title         
                expansions_read_lists.update({expansion.repeat_id: (genotype1,genotype2,title,read_list_A1,read_list_A2)})

    return expansions_read_lists


def plotHistogram(expansion_object: Expansion, plotfolder, bool_altclustering):
    
    save_place = plotfolder + "/" + expansion_object.title + '_hist.svg'
    ref_size = expansion_object.wt_size
    min_pathogenic = expansion_object.pathogenic_range
    repeat_unit = expansion_object.repeat_unit
    ref_motif = expansion_object.ref_motif

    if bool_altclustering:
        allele1_size = expansion_object.new_allele1
        allele2_size = expansion_object.new_allele2
        read_size_lists = expansion_object.new_read_list
    else:    
        allele1_size = expansion_object.allele1_size
        allele2_size = expansion_object.allele2_size
        read_size_lists = expansion_object.read_list

    # get x ranges
    max_x = max(max(read_size_lists[0]), ref_size * len(ref_motif))
    min_x = min(min(read_size_lists[0]), ref_size * len(ref_motif))

    plt.figure(figsize=(16, 11), dpi=300)
    if allele1_size == allele2_size:
        plt.hist(read_size_lists[0], color='Black', density=False, bins=np.arange(min(read_size_lists[0]), max(read_size_lists[0]) + 3, 3), label='Allele 1 + 2', rwidth=0.85)
        plt.axvline(allele1_size, color='k', linestyle='dashed', linewidth=1, label=allele1_size)
        plt.axvline(allele2_size, color='k', linestyle='dashed', linewidth=1, label=allele2_size)
    else:
        plt.hist([read_size_lists[0], read_size_lists[1]], color=['Black', 'Darkgray'], label=['Allele 1', 'Allele 2'], density=False,
                 bins=np.arange(min(read_size_lists[0]+read_size_lists[1]), max(read_size_lists[0]+read_size_lists[1]) + 3, 3), rwidth=0.85)
        plt.axvline(allele1_size, color='black', linestyle='dashed', linewidth=1, label=allele1_size)
        plt.axvline(allele2_size, color='darkgray', linestyle='dashed', linewidth=1, label=allele2_size)

    plt.ylabel('Number of Reads')
    plt.xlabel('Size(bp)')
    plt.title(expansion_object.title)
    plt.axvline(ref_size * len(ref_motif), color='blue', linewidth=1, label="wt size: " + str(ref_size * len(ref_motif)))
    if min_pathogenic != 'NA':
        plt.axvline(int(min_pathogenic) * len(ref_motif), color='red', linewidth=1, label="min pathogenic: " + str(int(min_pathogenic) * len(ref_motif)))
    plt.xlim((min_x - 0.5 * (max_x - min_x)), (max_x + 0.5 * (max_x - min_x)))
    plt.legend()
    plt.savefig(save_place, format="svg")
    plt.close()


def process_motifs(seqlist, motif, flank_length):
    """Process sequence to identify motifs and their positions"""
    # Remove flanking regions
    
    # Initialize dictionary for motifs and their positions
    motifs_dict = {}
    new_seq = []

    for i,seq in enumerate(seqlist):

        if seq == motif:
            
            if motif in motifs_dict:
                motifs_dict[motif].append(len(new_seq)+flank_length)
            else:
                motifs_dict.update({motif:[len(new_seq)+flank_length]})
                
            new_seq += motif
        elif seq != "":
            
            if seq in motifs_dict:
                motifs_dict[seq].append(len(new_seq)+flank_length)
            else:
                motifs_dict.update({seq:[len(new_seq)+flank_length]})
                
            new_seq += seq
    return motifs_dict

def alleleVisualiser(repeat_unit, flank_length, title, output_folder, chrom, start, end, genome_path, sequences):
    """Create allele visualization plot"""
    # Get reference sequence
    with pysam.FastaFile(genome_path) as fasta:
        reference_sequence = fasta.fetch(chrom, int(start), int(end))
    
    motif = repeat_unit.lower()
    # Process motifs and create visualization
    motifs_dict_list = []
    # Create and save plot
    for seq in sequences:
        expansion_seq = seq.sequence[flank_length:-flank_length]
        motif_matches = re.split("(" +motif +")", expansion_seq)
        motifs_dict_list.append(process_motifs(motif_matches, motif, flank_length))

    plot_alleles(motifs_dict_list, title, output_folder, flank_length)

def rectangleMaker(motif_colors, motif_coord_list, size, flank_length, motif, methylation_data=None, read_names=None):
    patches = []
    color_list = []
    labels = {}
    
    for i, motifs_dict in enumerate(motif_coord_list):
        # Create base rectangles for sequence elements
        for motif_type, positions in motifs_dict.items():
            for pos in positions:
                rect = Rectangle((pos, 10*i), len(motif_type), 5)
                patches.append(rect)
                
                if motif_type == motif:
                    color_list.append(motif_colors["motif"])
                    labels[motif_colors["motif"]] = "Repeat Motif"
                elif len(motif_type) == len(motif):
                    color_list.append(motif_colors["non_motif"])
                    labels[motif_colors["non_motif"]] = "Non-motif Sequence"
                else:
                    color_list.append(motif_colors["flank"])
                    labels[motif_colors["flank"]] = "Flanking Sequence"
        
        # Add methylation visualization if available
        if methylation_data and read_names and i < len(read_names)-1:  # Skip reference sequence
            read_name = read_names[i]
            if read_name in methylation_data:
                for methyl_call in methylation_data[read_name]:
                    plot_pos = methyl_call.position - flank_length
                    if 0 <= plot_pos <= len(motifs_dict.get(motif, [])):
                        rect = Rectangle((plot_pos, 10*i), 2, 5)
                        patches.append(rect)
                        
                        if methyl_call.is_methylated:
                            color_list.append(motif_colors["methylated"])
                            labels[motif_colors["methylated"]] = "Methylated CpG"
                        else:
                            color_list.append(motif_colors["unmethylated"])
                            labels[motif_colors["unmethylated"]] = "Unmethylated CpG"
    
    return patches, color_list, labels


def plot_alleles(motifs_dict_list, title, output_folder, flank_length):
    """Create allele visualization plot using rectangles"""
    plt.figure(figsize=(20, len(motifs_dict_list)), dpi=300)
    
    # Define colors for different sequence types
    motif_colors = {
        "flank": "lightgray",
        "motif": "red",
        "non_motif": "gray"
    }
    
    # Calculate max sequence length
    max_size = 0
    for d in motifs_dict_list:
        size = 0
        for key, coords in d.items():
            if isinstance(coords, list):
                size += len(key) * len(coords)
        max_size = max(max_size, size + 2 * flank_length)
    
    # Create rectangles for visualization
    patches, colors, labels = rectangleMaker(
        motif_colors, 
        motifs_dict_list, 
        max_size, 
        flank_length,
        title.split('_')[1]  # Use repeat_id as motif
    )
    
    # Create patch collection and add to plot
    collection = PatchCollection(patches, facecolors=colors)
    ax = plt.gca()
    ax.add_collection(collection)
    
    # Set plot parameters
    plt.xlim(-10, max_size + 10)
    plt.ylim(-5, 10 * len(motifs_dict_list))
    plt.title(title)
    
    # Create legend
    legend_elements = [plt.Rectangle((0,0), 1, 1, facecolor=color, label=label)
                      for color, label in labels.items()]
    plt.legend(handles=legend_elements)
    
    # Save plot
    plt.savefig(f"{output_folder}/{title}_alleles.svg", format="svg", bbox_inches='tight')
    plt.close()
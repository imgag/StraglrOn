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


def process_motifs(sequence: str, motif: str, flank_length: int) -> dict:
    """
    Workflow  from here with read list:
        - right and left flanks are cut off -> list (cut_read_list)
        - motif is used in splitting (while keeping the splitting motif) the string generating list of motif -> list (read_motif_list)
        - individual list entries are used to form a new string and use the iterator to find motif start coordinates, alternate motifs and motif lengths -> motif_coord
        - dictionary for each motif is created with a list of the coordinates on the original read_list -> motif_dict_list
    """

    # Remove flanking regions
    core_seq = sequence[flank_length:-flank_length]
    
    # Split at motif boundaries while keeping the motifs
    clean_motif = motif.replace('*', '')  # Remove * from motif
    escaped_motif = re.escape(clean_motif)  # Escape any remaining special characters     
    parts = re.split("(" + escaped_motif + ")", core_seq)
    
    # Track positions and build motif dictionary
    pos = flank_length  # Start after left flank
    motifs_dict = {}
    
    for part in parts:
        if part:  # Skip empty strings from consecutive motifs
            if part == motif:
                if motif in motifs_dict:
                    motifs_dict[motif].append(pos)
                else:
                    motifs_dict[motif] = [pos]
            else:
                if part in motifs_dict:
                    motifs_dict[part].append(pos)
                else:
                    motifs_dict[part] = [pos]
            pos += len(part)
                
    return motifs_dict


def alleleVisualiser(motif, flank_length, title, output_folder, chromosome, start, end, reference_genome, sequences):
    """
    Create waterfall plot visualization including methylation data
    Args:
        sequences: List[RepeatSequence] containing sequence and methylation data
    """
    chr = chromosome
    start = int(start)
    end = int(end)
    
    refGen = pysam.FastaFile(reference_genome)
    reference_sequence = refGen.fetch(chr, start, end)
    
    # Sort sequences by repeat size
    sequences.sort(key=lambda x: x.repeat_size)
    
    # Calculate maximum sequence length
    max_seq_length = max(len(seq.sequence) for seq in sequences)
    # Also consider reference sequence length
    max_seq_length = max(max_seq_length, len(reference_sequence))
    
    motif_colors = {
        "flank": "grey",
        "motif": "green",
        "non_motif": "red",
        "methylated": "blue",
        "unmethylated": "yellow"
    }
    
    motif = motif.lower()

    # Process sequences to get motifs
    motifs_dict_list = []
    for seq in sequences:
        motifs_dict_list.append(process_motifs(seq.sequence, motif, flank_length))

    # Add reference sequence motifs
    motifs_dict_list.append(process_motifs(reference_sequence.lower(), motif, flank_length))

    # Create visualization
    fig, ax = plt.subplots(figsize=(16, 11), dpi=300)
    plot_title = title + "_" + motif
    plt.xlim([0, max_seq_length])  # Use calculated max length
    plt.ylim([0, (len(sequences)+2)*10])  # +2 for reference sequence and padding
    
    patches, color_list, labels = rectangleMaker(
        motif_colors, 
        motifs_dict_list, 
        sequences,  # Pass the sorted sequences directly
        flank_length,
        motif
    )
    
    our_cmap = ListedColormap(color_list)
    patches_collection = PatchCollection(patches, cmap=our_cmap)
    patches_collection.set_array(np.arange(len(patches)))
    ax.add_collection(patches_collection)
    handles = []

    for patch_color in set(color_list):
        patch = matplotlib.patches.Patch(color=patch_color, label=labels[patch_color])
        handles.append(patch)
        
    plt.legend(handles=handles, loc="lower right")
    plt.title(plot_title)
    plt.axis("off")
    plt.savefig((output_folder + "/" + title + ".svg").replace(":", "_"), format="svg")
    plt.close()


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
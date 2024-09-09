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


def alleleVisualiser(fasta_file, motif, flank_length, title, output_folder, chromosome, start, end, reference_genome):
    
    chr = chromosome
    start = int(start)
    end = int(end)
    
    refGen = pysam.FastaFile(reference_genome)
    
    reference_sequence = refGen.fetch(chr,start,end)
    
    fasta_sequences = list(SeqIO.parse(open(fasta_file),'fasta'))
    read_list = []
    for fasta in fasta_sequences:
        read_list.append(str(fasta.seq))
    
    motif_colors = ["grey", "green", "red"]
    
    #readlist sorted by length for waterfall visualisation
    read_list.sort(key=lambda x: len(x))

    motif = motif.lower()
    flank_length = 25

    '''
    Workflow  from here with read list:
        - right and left flanks are cut off -> list (cut_read_list)
        - motif is used in splitting (while keeping the splitting motif) the string generating list of motif -> list (read_motif_list)
        - individual list entries are used to form a new string and use the iterator to find motif start coordinates, alternate motifs and motif lengths -> motif_coord
        - dictionary for each motif is created with a list of the coordinates on the original read_list -> motif_dict_list
    '''

    cut_read_list = []
    for x in read_list:
        cut = x[flank_length:-flank_length]
        cut_read_list.append(cut)

    cut_read_list.insert(len(cut_read_list), reference_sequence.lower())
    read_motif_list = []
    for x in cut_read_list:
        read_motif = re.split("(" +motif +")", x)
        read_motif_list.append(read_motif)

    motifs_dict_list = []

    for x in read_motif_list:
        
        new_seq = []
        motif_coords = []
        motifs_dict = {}

        for i, s in enumerate(x):

            if x[i] == motif:
                
                if motif in motifs_dict:
                    motifs_dict[motif].append(len(new_seq)+flank_length)
                else:
                    motifs_dict.update({motif:[len(new_seq)+flank_length]})
                    
                new_seq += motif
            elif x[i] != "":
                
                if x[i] in motifs_dict:
                    motifs_dict[x[i]].append(len(new_seq)+flank_length)
                else:
                    motifs_dict.update({x[i]:[len(new_seq)+flank_length]})
                    
                new_seq += x[i]
                
        motifs_dict_list.append(motifs_dict)

    # For setting of the scale of the x axis as well as the positioning of the right flank patch
    size = len(max(cut_read_list, key=len))+flank_length*2

    # define Matplotlib figure and axis
    fig, ax = plt.subplots(figsize=(16, 11), dpi=300)
    plot_title = fasta_sequences[0].description.split(" ")[1] + "_" + motif
    # Plot dimension depending on number and maximum length
    plt.xlim([0,size])
    plt.ylim([0, (len(read_list)+1)*10])
    patches_list, color_list, labels = rectangleMaker(motif_colors, motifs_dict_list, size, flank_length, motif)
    our_cmap = ListedColormap(color_list)
    patches_collection = PatchCollection(patches_list, cmap=our_cmap)
    patches_collection.set_array(np.arange(len(patches_list)))
    ax.add_collection(patches_collection)
    handles=[]

    for patch_color in set(color_list):
        patch = matplotlib.patches.Patch(color=patch_color, label=labels[patch_color])
        handles.append(patch)
        
    plt.legend(handles=handles, loc="lower right")
    plt.title(plot_title)
    plt.axis("off")
    plt.show()
    plt.savefig((output_folder + "/" + title + ".svg").replace(":", "_"), format="svg")
    plt.close()


def rectangleMaker(motif_colors, motif_coord_list, size, flank_length, motif):
    patches = []
    color_list = []
    labels = {}
    
    for i, x in enumerate(motif_coord_list):
        color_chooser = 0
        color = motif_colors[color_chooser]
        patches.append(Rectangle((0, i*10), flank_length, 5, label="Flank"))
        color_list.append(color)
        labels.update({color: "Flank"})
        
        patches.append(Rectangle((size-flank_length, i*10), flank_length, 5, label="Flank"))
        color_list.append(color)
        labels.update({color:"Flank"})
            
        for key in x:
            if key == motif:
                color = motif_colors[1]
                for coord in x[key]:
                    patches.append(Rectangle((coord, 10*i), len(key), 5, label=key))
                    color_list.append(color)
                labels.update({color:key})
            
            else:
                color = motif_colors[2]
                for coord in x[key]:
                    patches.append(Rectangle((coord, 10*i), len(key), 5, label="non-motif"))
                    color_list.append(color)
                labels.update({color:"non-motif"})
        
    return patches, color_list, labels
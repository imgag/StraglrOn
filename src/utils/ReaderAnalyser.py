import csv
import pysam
from pathlib import Path
from .Structures import Expansion
from .Structures import Locus
import numpy as np
from sklearn.mixture import GaussianMixture


def resultBedReader(file, loci_dict):
    with open(file) as straglr:

        straglr_reader = csv.reader(straglr, delimiter='\t')
        expansion_list: list[Expansion] = []
        sample_id = Path(file).stem

        # loci   allele1 allele2 YES/NO  Referenz    size_difference    PathRange
        for straglr in straglr_reader:
            if straglr[0] != '#chrom':
                # Variables are assigned based on the coordinates in the straglr reader
                coords = straglr[0] + ":" + straglr[1] + "-" + straglr[2]
                chr = straglr[0]
                start = straglr[1]
                end = straglr[2]
                loci = loci_dict[coords][0].name
                reference_size = float(loci_dict[coords][0].reference_size)
                ref_motif = loci_dict[coords][0].motif
                if straglr[14].strip() == "":
                    int_path_range = 'NA'
                else:
                    int_path_range = straglr[14]
                motif = straglr[3]
                alleles = 0
                # "-" means that there was no coverage in this region -> skip this loci
                if straglr[4] == '-':
                    continue

                # "-" means that straglr assigned no different number for second allele -> allele 1 = allele 2; In the following 2 different allele sizes were assigned
                elif straglr[8] != '-':
                    allele1 = round(float(straglr[4]))
                    allele2 = round(float(straglr[7]))
                    copy_number_1 = straglr[5]
                    copy_number_2 = straglr[8]
                    allele1_support = straglr[6]
                    allele2_support = straglr[9]
                    alleles = 2

                else:
                    allele1 = round(float(straglr[4]))
                    allele2 = round(float(straglr[4]))
                    copy_number_1 = straglr[5]
                    copy_number_2 = straglr[5]
                    allele1_support = straglr[6]
                    allele2_support = straglr[6]
                    alleles = 1

                expansion_object = Expansion(chr, start, end, loci, motif, ref_motif, allele1, allele2, reference_size, int_path_range, copy_number_1, copy_number_2, allele1_support,
                                             allele2_support, sample_id)
                analyseGenotype(expansion_object, alleles)
                expansion_list.append(expansion_object)

        # List of expansions, though not pathogenic by definition, are sorted by chromosome or the normalized size difference between sample and reference length on hg38 in descending order

    return expansion_list


def analyseGenotype(expansion_object: Expansion, alleles):
    if alleles == 0:
        # repeat not covered
        expansion_object.in_pathogenic_range = 'NA'
        expansion_object.size_difference = np.nan
    elif alleles == 1:

        if expansion_object.pathogenic_range == 'NA':
            expansion_object.in_pathogenic_range = 'NA'
            if expansion_object.allele1_size > expansion_object.allele2_size:
                expansion_object.size_difference = expansion_object.allele1_size - expansion_object.wt_size
            else:
                expansion_object.size_difference = expansion_object.allele2_size - expansion_object.wt_size

                # Check for pathogenic expansion size happens here; Size of sample expansion checked against lower boundary of pathogenic range from literature; Assignment of respective loci and sizes to corresponding lists
        elif expansion_object.allele1_size > int(expansion_object.pathogenic_range) or expansion_object.allele2_size > int(expansion_object.pathogenic_range):
            expansion_object.in_pathogenic_range = 'Yes'
            if expansion_object.allele1_size > expansion_object.allele2_size:
                expansion_object.size_difference = expansion_object.allele1_size - expansion_object.wt_size
            else:
                expansion_object.size_difference = expansion_object.allele2_size - expansion_object.wt_size

        # If size of sample expansion is still below pathogenic range, they are added into general list
        else:
            expansion_object.in_pathogenic_range = 'No'
            if expansion_object.allele1_size > expansion_object.allele2_size:
                expansion_object.size_difference = expansion_object.allele1_size - expansion_object.wt_size
            else:
                expansion_object.size_difference = expansion_object.allele2_size - expansion_object.wt_size

    else:
        if expansion_object.pathogenic_range == 'NA':
            expansion_object.in_pathogenic_range = 'NA'
            # Same check and assignment as above
            if expansion_object.allele1_size > expansion_object.allele2_size:
                expansion_object.size_difference = expansion_object.allele1_size - expansion_object.wt_size
            else:
                expansion_object.size_difference = expansion_object.allele2_size - expansion_object.wt_size

        elif expansion_object.allele1_size > int(expansion_object.pathogenic_range):
            expansion_object.in_pathogenic_range = 'Yes'
            expansion_object.size_difference = expansion_object.allele1_size - expansion_object.wt_size

        else:
            expansion_object.in_pathogenic_range = 'No'
            if expansion_object.allele1_size > expansion_object.allele2_size:
                expansion_object.size_difference = expansion_object.allele1_size - expansion_object.wt_size
            else:
                expansion_object.size_difference = expansion_object.allele2_size - expansion_object.wt_size


# TODO: use pandas for parsing
def lociBedReader(bedFile):
    Loci = {}
    with open(bedFile) as lociCoords:
        lociReader_reader = csv.reader(lociCoords, delimiter='\t')
        header = []
        for locus in lociReader_reader:
            if locus[0] == '#chr':
                header = list(locus)
                header[0] = header[0][1:]
                print(header)
            else:
                key = locus[header.index("chr")] + ":" + locus[header.index("start")] + "-" + locus[header.index("end")]
                new_object = Locus(locus[header.index("repeat_id")], locus[header.index("chr")], locus[header.index("start")], locus[header.index("end")],
                                   locus[header.index("repeat_motif")], locus[header.index("repeat_id")], "", locus[header.index("ref_size")],
                                   "NA", "NA")
                if key in Loci:
                    Loci[key].append(new_object)
                else:
                    Loci.update({key: [new_object]})
    return Loci


def newGenotyping(expansion_object: Expansion, cutoff, new: bool):
    concat_reads = []
    number_of_reads = 0

    if new:
        for new_readlist in expansion_object.new_read_list:
            concat_reads += new_readlist
            number_of_reads = len(concat_reads)
            
    else: 
        for original_readlist in expansion_object.read_list:
            concat_reads += original_readlist
            number_of_reads = len(concat_reads)

    rearanged_array = np.array(concat_reads).reshape(-1,1)
    X = rearanged_array

    N = np.arange(1, 4)
    M = np.arange(1, 3)

    # fit models with 1-3 components

    if len(np.unique(X)) < 2:
        models = [GaussianMixture(1, covariance_type='full', init_params="kmeans", max_iter=500).fit(X)]

    elif len(np.unique(X)) == 2:
        models = [GaussianMixture(m, covariance_type='full', init_params="kmeans", max_iter=500).fit(X)
                  for m in M]
    else:
        models = [GaussianMixture(n, covariance_type='full', init_params="kmeans", max_iter=500).fit(X)
                  for n in N]

    BIC = [m.bic(X) for m in models]

    best_model = models[np.argmin(BIC)]

    cluster_assignment = best_model.predict(X)

    cluster1 = []
    cluster2 = []
    cluster3 = []
    clusters = []

    for i in range(len(X)):
        
        if cluster_assignment[i] == 0:
            cluster1.append(X[i][0])
        if cluster_assignment[i] == 1:
            cluster2.append(X[i][0])
        if cluster_assignment[i] == 2:
            cluster3.append(X[i][0])

    clusters.append(cluster1)
    clusters.append(cluster2)
    clusters.append(cluster3)

    clusters.sort(key=len, reverse=True)

    cluster1 = clusters[0]
    cluster2 = clusters[1]
    cluster3 = clusters[2]
    if cluster3:
        
        distance13 = abs(np.mean(cluster3) - np.mean(cluster1))
        distance12 = abs(np.mean(cluster3) - np.mean(cluster2))
        if distance12 < distance13:
            if distance12 < max(5, 2*np.var(cluster2)):
                cluster2 += cluster3
        else:
            if distance13 < max(5, 2*np.var(cluster2)):
                cluster1 += cluster3
    allele1 = cluster1
    allele2 = cluster2

    if allele1 and allele2:

        expansion_object.new_read_list = [allele1, allele2]
        expansion_object.new_allele1 = np.median(allele1)
        expansion_object.new_allele2 = np.median(allele2)

        if expansion_object.new_allele1 > expansion_object.new_allele2:
            expansion_object.new_size_difference = expansion_object.new_allele1 - expansion_object.wt_size
            if expansion_object.pathogenic_range == "NA":
                expansion_object.new_in_pathogenic_range = "NA"
            elif expansion_object.new_allele1 > int(expansion_object.pathogenic_range):
                expansion_object.new_in_pathogenic_range = "Yes"
            else:
                expansion_object.new_in_pathogenic_range = "No"

        else:
            expansion_object.new_size_difference = expansion_object.new_allele2 - expansion_object.wt_size
            if expansion_object.pathogenic_range == "NA":
                expansion_object.new_in_pathogenic_range = "NA"
            elif expansion_object.new_allele2 > int(expansion_object.pathogenic_range):
                expansion_object.new_in_pathogenic_range = "Yes"
            else:
                expansion_object.new_in_pathogenic_range = "No"

        expansion_object.new_copy_numberA1 = expansion_object.new_allele1 / len(expansion_object.repeat_unit)
        expansion_object.new_copy_numberA2 = expansion_object.new_allele2 / len(expansion_object.repeat_unit)
        expansion_object.new_allele1_support = len(allele1)
        expansion_object.new_allele2_support = len(allele2)
    else:
        expansion_object.new_read_list = [allele1 + allele2]
        expansion_object.new_allele1 = np.median(allele1)
        expansion_object.new_allele2 = np.median(allele1)
        expansion_object.new_size_difference = expansion_object.new_allele1 - expansion_object.wt_size

        if expansion_object.pathogenic_range == "NA":
            expansion_object.new_in_pathogenic_range = "NA"
        elif expansion_object.new_allele1 > int(expansion_object.pathogenic_range):
            expansion_object.new_in_pathogenic_range = "Yes"
        else:
            expansion_object.new_in_pathogenic_range = "No"

        expansion_object.new_copy_numberA1 = expansion_object.new_allele1 / len(expansion_object.repeat_unit)
        expansion_object.new_copy_numberA2 = expansion_object.new_allele2 / len(expansion_object.repeat_unit)
        expansion_object.new_allele1_support = len(allele1)

    if number_of_reads > 3 * cutoff:

        for new_list in expansion_object.new_read_list:

            if len(new_list) <= cutoff:
                expansion_object.new_read_list.remove(new_list)
                newGenotyping(expansion_object, cutoff, True)


def expansionScorer(expansion: Expansion, new_clustering: bool):
    if expansion.pathogenic_range != "NA":
        if new_clustering:
            score = (expansion.new_size_difference) / (int(expansion.pathogenic_range) - expansion.wt_size)
            expansion.new_norm_score = round(score, 4)

        else:
            score = (expansion.size_difference) / (int(expansion.pathogenic_range) - expansion.wt_size)
            expansion.norm_score = round(score, 4)
    else:
        expansion.norm_score = None
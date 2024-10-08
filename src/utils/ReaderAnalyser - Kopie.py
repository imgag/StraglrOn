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
        
        #loci   allele1 allele2 YES/NO  Referenz    size_difference    PathRange
        for straglr in straglr_reader:
            if straglr[0]!='#chrom':
                #Variables are assigned based on the coordinates in the straglr reader
                coords = straglr[0]+straglr[1]+straglr[2]
                chr = straglr[0]
                start = straglr[1]
                end = straglr[2]
                loci = loci_dict[coords][0].name
                reference_size = int(loci_dict[coords][0].reference_size)
                int_path_range = loci_dict[coords][0].pathogenic_range
                motif = straglr[3]
                alleles = 0
                # "-" means that there was no coverage in this region
                if straglr[4] == '-':
                    allele1 = np.nan
                    allele2 = np.nan
                    copy_number_1 = np.nan
                    copy_number_2 = np.nan
                    allele1_support = 0
                    allele2_support = 0
                    alleles = 0
                #"-" means that straglr assigned no different number for second allele -> allele 1 = allele 2; In the following 2 different allele sizes were assigned
                elif straglr[8]!='-':
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
                            
                expansion_object = Expansion(chr, start, end, loci, motif, allele1, allele2, reference_size, int_path_range, copy_number_1, copy_number_2, allele1_support, allele2_support, sample_id)
                analyseGenotype(expansion_object, alleles)
                expansion_list.append(expansion_object)
                
        #List of expansions, though not pathogenic by definition, are sorted by chromosome or the normalized size difference between sample and reference length on hg38 in descending order         
       
    return expansion_list

def analyseGenotype(expansion_object: Expansion, alleles):

    if alleles == 0:
        # repeat not covered
        expansion_object.in_pathogenic_range = 'NA'
        expansion_object.size_difference = np.nan
    elif alleles == 1:
        
        if expansion_object.pathogenic_range == 'NA':
            expansion_object.in_pathogenic_range ='NA'
            if expansion_object.allele1_size > expansion_object.allele2_size:
                expansion_object.size_difference = expansion_object.allele1_size-expansion_object.wt_size
            else:
                expansion_object.size_difference = expansion_object.allele2_size-expansion_object.wt_size            

        #Check for pathogenic expansion size happens here; Size of sample expansion checked against lower boundary of pathogenic range from literature; Assignment of respective loci and sizes to corresponding lists
        elif expansion_object.allele1_size > int(expansion_object.pathogenic_range) or expansion_object.allele2_size > int(expansion_object.pathogenic_range):
            expansion_object.in_pathogenic_range = 'Yes'
            if expansion_object.allele1_size > expansion_object.allele2_size:
                expansion_object.size_difference = expansion_object.allele1_size-expansion_object.wt_size
            else:
                expansion_object.size_difference = expansion_object.allele2_size-expansion_object.wt_size
                
        #If size of sample expansion is still below pathogenic range, they are added into general list
        else:   
            expansion_object.in_pathogenic_range = 'No'
            if expansion_object.allele1_size > expansion_object.allele2_size:
                expansion_object.size_difference = expansion_object.allele1_size-expansion_object.wt_size
            else:
                expansion_object.size_difference = expansion_object.allele2_size-expansion_object.wt_size
    
    else:            
        if expansion_object.pathogenic_range == 'NA':
            expansion_object.in_pathogenic_range = 'NA'
            #Same check and assignment as above
            if expansion_object.allele1_size>expansion_object.allele2_size:
                expansion_object.size_difference = expansion_object.allele1_size-expansion_object.wt_size
            else:
                expansion_object.size_difference = expansion_object.allele2_size-expansion_object.wt_size
                            
        elif expansion_object.allele1_size> int(expansion_object.pathogenic_range):
            expansion_object.in_pathogenic_range = 'Yes'
            expansion_object.size_difference = expansion_object.allele1_size-expansion_object.wt_size
            
        else:
            expansion_object.in_pathogenic_range = 'No'
            if expansion_object.allele1_size>expansion_object.allele2_size:
                expansion_object.size_difference = expansion_object.allele1_size-expansion_object.wt_size
            else:
                expansion_object.size_difference = expansion_object.allele2_size-expansion_object.wt_size
                
def lociBedReader(bedFile):
    Loci = {}
    with open(bedFile) as lociCoords:
        lociReader_reader=csv.reader(lociCoords, delimiter='\t')
        for locus in lociReader_reader:
            if locus[0]!='#chrom':
                key = locus[0] + locus[1] + locus[2]
                new_object = Locus(locus[4],locus[0],locus[1],locus[2],locus[3],locus[4],locus[5],locus[6],locus[7],locus[8])
                if key in Loci:
                    Loci[key].append(new_object)
                else:
                    Loci.update({key:[new_object]})
    return Loci

def newGenotyping(expansion_object: Expansion, cutoff, new: bool):
        
    concat_reads = []
    number_of_reads = 0
    
    if new:
    
        for new_readlist in expansion_object.new_read_list:
            #print("Oldlist: ")
            #print(original_readlist)
            concat_reads += new_readlist
            number_of_reads = len(concat_reads)
            #print(number_of_reads)
            
    else: 
    
        for original_readlist in expansion_object.read_list:
            #print("Oldlist: ")
            #print(original_readlist)
            concat_reads += original_readlist
            number_of_reads = len(concat_reads)
            #print(number_of_reads)
            
    rearanged_array = np.array(concat_reads).reshape(-1,1)
    #print(expansion_object.repeat_id, expansion_object.repeat_unit, " STARTS HERE: ")
    #print(rearanged_array)
    
    X = rearanged_array
    
    N = np.arange(1,4)
    M = np.arange(1,3)
    
    # fit models with 1-2 components
    
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
    
    #print(cluster_assignment)
    #print(X)
    #print(np.std(X))

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
            
    #print("This is cluster1: ", cluster1)
    #print("This is cluster2: ", cluster2)
    #print("This is cluster3: ", cluster3)
    
    clusters.append(cluster1)       
    clusters.append(cluster2)
    clusters.append(cluster3)
    
    clusters.sort(key=len, reverse=True)        
    
    cluster1 = clusters[0]
    cluster2 = clusters[1]
    cluster3 = clusters[2]
    
    #print("Longest cluster: ", cluster1)
    #print("Middle cluster: ", cluster2)
    #print("Shortest cluster: ",cluster3)
    
    if cluster3:
        distance13 = abs(np.mean(cluster3) - np.mean(cluster1))
        distance12 = abs(np.mean(cluster3) - np.mean(cluster2))
        if  distance12 < distance13: 
            if distance12 < max(5, 2*np.var(cluster2)):
                cluster2 += cluster3
                #print("Cluster 3 wurde cluster 2 geadded werden") 
        else:
            if distance13 < max(5, 2*np.var(cluster2)):
                cluster1 += cluster3
                #print("Cluster 3 wurde cluster 1 geadded werden")
        #print(np.mean(cluster3), " ", np.var(cluster3))
        #print(np.mean(cluster2), " ", np.var(cluster2))
        #print(np.mean(cluster1), " ", np.var(cluster1))
    
    allele1 = cluster1
    allele2 = cluster2
    
    #print(allele1)
    #print(allele2)
    
    if allele1 and allele2: 
          
        expansion_object.new_read_list = [allele1,allele2]
        expansion_object.new_allele1 = np.median(allele1)
        expansion_object.new_allele2 = np.median(allele2)
        
        if expansion_object.new_allele1 > expansion_object.new_allele2:
            expansion_object.new_size_difference = expansion_object.new_allele1-expansion_object.wt_size
            if expansion_object.pathogenic_range == "NA":
                expansion_object.new_in_pathogenic_range = "NA"
            elif expansion_object.new_allele1 > int(expansion_object.pathogenic_range):
                expansion_object.new_in_pathogenic_range = "Yes"
            else: 
                expansion_object.new_in_pathogenic_range = "No"
                
        else:
            expansion_object.new_size_difference = expansion_object.new_allele2-expansion_object.wt_size
            if expansion_object.pathogenic_range == "NA":
                expansion_object.new_in_pathogenic_range = "NA"
            elif expansion_object.new_allele2 > int(expansion_object.pathogenic_range):
                expansion_object.new_in_pathogenic_range = "Yes"
            else: 
                expansion_object.new_in_pathogenic_range = "No"
                
        expansion_object.new_copy_numberA1 = expansion_object.new_allele1/len(expansion_object.repeat_unit)
        expansion_object.new_copy_numberA2 = expansion_object.new_allele2/len(expansion_object.repeat_unit)
        expansion_object.new_allele1_support = len(allele1)
        expansion_object.new_allele2_support =  len(allele2)
        #print("New: " + expansion_object.title, expansion_object.new_read_list)
    else:
        expansion_object.new_read_list = [allele1+allele2]
        expansion_object.new_allele1 = np.median(allele1)
        expansion_object.new_allele2 = np.median(allele1)
        expansion_object.new_size_difference = expansion_object.new_allele1-expansion_object.wt_size
        
        if expansion_object.pathogenic_range == "NA":
                expansion_object.new_in_pathogenic_range = "NA"
        elif expansion_object.new_allele1 > int(expansion_object.pathogenic_range):
            expansion_object.new_in_pathogenic_range = "Yes"
        else: 
            expansion_object.new_in_pathogenic_range = "No"
            
        expansion_object.new_copy_numberA1 = expansion_object.new_allele1/len(expansion_object.repeat_unit)
        expansion_object.new_copy_numberA2 = expansion_object.new_allele2/len(expansion_object.repeat_unit)
        expansion_object.new_allele1_support = len(allele1)
        expansion_object.new_allele2_support =  len(allele1)
        
        #print(expansion_object.title, expansion_object.new_read_list)

    #print("The length of the current lists is {}".format([len(n) for n in expansion_object.new_read_list]))

    #temp_list = expansion_object.new_read_list
    
    #print(len(expansion_object.new_read_list))
    
    if number_of_reads > 3*cutoff:

        for new_list in expansion_object.new_read_list:
            
            if len(new_list) <= cutoff:
                expansion_object.new_read_list.remove(new_list)
                #expansion_object.read_list = temp_list
                newGenotyping(expansion_object, cutoff, True)
        
def expansionScorer(expansion : Expansion, new_clustering : bool):
    #score = x-xmin/xmax-xmin
    if expansion.pathogenic_range != "NA":
        if new_clustering:
            score = (expansion.new_size_difference)/(int(expansion.pathogenic_range)-expansion.wt_size)
            expansion.new_norm_score = round(score, 4)
            
        else:
            score = (expansion.size_difference)/(int(expansion.pathogenic_range)-expansion.wt_size)
            #print(expansion.in_pathogenic_range, score, expansion.size_difference, int(expansion.pathogenic_range))     
            expansion.norm_score = round(score, 4) 
    else:
        expansion.norm_score = None  
 
'''    
def bamfileReader(bamFile, pathogenics):
    samfile = pysam.AlignmentFile(bamFile, "rb")
    counter = 0
    expansion = pathogenics
    chromosome = expansion.chr
    start = int(expansion.start)
    end = int(expansion.end)
    ru = expansion.repeat_unit
    print(expansion.repeat_id, chromosome, start, end, ru)
    for key, value in expansion.read_dict.items():
        print(key, value, len(expansion.read_dict))
    for read in samfile.fetch(chromosome, start, end):
        try:
            print(read.query_name, expansion.read_dict[read.query_name])
            print(read.query_sequence[expansion.read_dict[read.query_name]:(expansion.read_dict[read.query_name]+50)])
            print()
        
        except KeyError:
            print("Key " + read.query_name + " not found!")
        counter += 1
    samfile.close()
    print(counter)


def cyclicPermutation(motif):
    size = len(motif)
    motifs = deque(motif)
    perms = []
    
    for i in range(size):
        motifs.rotate(1)
        perms.append(''.join(list(motifs)))    

    return permm
'''
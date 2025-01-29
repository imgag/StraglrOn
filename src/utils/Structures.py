#Expansion class used for variable and function storage for analysis and presentation        
class Expansion:
    def __init__(self, chr, start, end, locus, repeat_id, repeat_unit, ref_motif, allele1_size, allele2_size, wt_size, pathogenic_range, copy_numberA1, copy_numberA2, allele1_support, allele2_support, sample_id):
        # Data imported from straglr output bamfile in ResultBedReader = required
        self.chr = chr
        self.start = start
        self.end = end
        self.locus = locus
        self.repeat_id = repeat_id
        self.repeat_unit = repeat_unit
        self.ref_motif = ref_motif
        self.allele1_size = allele1_size
        self.allele2_size = allele2_size
        self.copy_numberA1 = copy_numberA1
        self.copy_numberA2 = copy_numberA2
        self.allele1_support = allele1_support
        self.allele2_support = allele2_support 
        self.wt_size = wt_size
        self.sample_id = sample_id
        
        # Data imported and inferred from loci bed file in ResultBedReader = required
        
        self.pathogenic_range = pathogenic_range
        
        # Data inferred in analyze genotype
        self.in_pathogenic_range = None
        self.size_difference = None
    
        # Data imported and inferred from tsv file in getHistData = optional
        self.title = ""  
        self.read_list = []
        self.read_dict = {}
        
        # Data inferred in expansionScorer = optional
        self.norm_score = None
        
        # Data inferred from newGenotyping = optional
        self.new_read_list = []  
        self.new_allele1 = None
        self.new_allele2 = None
        self.new_in_pathogenic_range = "None"
        self.new_size_difference = None
        self.new_allele1_support = None
        self.new_allele2_support = None
        self.new_copy_numberA1 = None
        self.new_copy_numberA2 = None
        self.new_allele1_support = None
        self.new_allele2_support =  None
        self.new_in_pathogenic_range = None
        self.new_norm_score = None
        

    @property
    def title(self):
        return self._title
    
    @title.setter
    def title(self, string: str):
        subs = string.replace(":", "_")
        self._title = subs
        
    
#Locus class imported from bed/txt file used in straglr run
class Locus:
    def __init__(self, name, chromosome, start, end, motif, locus, associated_disease, reference_size, normal_range, min_pathogenic):
        self.name = name
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.motif = motif
        self.locus = locus
        self.associated_disease = associated_disease
        self.reference_size = reference_size
        self.normal_range = normal_range
        self.min_pathogenic = min_pathogenic
        #self.pathogenic_motif = pathogenic_motiv
        #self.pathogenic_motif_range = pathogenic_motif_range

# Methylation call class used in extract_repeats.py
class MethylationCall:
    def __init__(self, position: int, is_methylated: bool, quality_score: float):
        self.position = position
        self.is_methylated = is_methylated
        self.quality_score = quality_score

class RepeatSequence:
    def __init__(self, locus: str, read_name: str, repeat_size: int, sequence: str, 
                 start_position: int, end_position: int, 
                 left_flank: str, repeat_sequence: str, right_flank: str,
                 methylation_calls=None):
        self.locus = locus
        self.read_name = read_name
        self.repeat_size = repeat_size
        self.sequence = sequence
        self.start_position = start_position
        self.end_position = end_position
        self.left_flank = left_flank
        self.repeat_sequence = repeat_sequence
        self.right_flank = right_flank
        self.methylation_calls = methylation_calls if methylation_calls else []
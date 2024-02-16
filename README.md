# StraglrON
Addon for Straglr by Thomas Braun.

Developed for better readablity and visualization. Base requirements are both files produced by Straglr (tsv and bed file) as well as the loci bed file used in the Straglr apllication. The loci bed file in this application needs to be supplied with additional information in the following format:
chr	start	end	motif	repeat_id	associated_disorder	ref_size	normal_range	pathogenic_range
If there is no information on the last column "NA" is sufficient. 

Base application will provide a .txt file with loci information sorted by chromosome. 

The optional modes "--hist" and "--alleles" will provide histograms and waterfall visualization of all locis returned by Straglr to be of interest. For the waterfall visualization you must provide the respective bam file and the used reference genome. 

The command "--altclust" activates an alternate grouping of read lengths which can be more sensitive than the default clustering Staglr employs.

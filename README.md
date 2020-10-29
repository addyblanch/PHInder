# PHInder - Pathogen Host Interaction Finder

Idenify putative interactions between bacteria and the host using Metatranscriptomic data.

Matching analysis as a single piece of code that will look for correlations between two data sets
that come from the same samples. We are particularly thinking about bacteria and gene data from
the same set of individuals.The code only looks at presence/absence rather than abundance.

Requires two data frames (output from a tool such as Salmon)

*1. Count data from bacterial RNASeq assignemnt*
*2. Count data from host RNASeq assignment*

See the example data for appropriate input. Both are data-frames with column 1
representing the bacteria/gene/whatever is being tested. The remaining columns should
represent each sample and should have the same labels to allow for matchings. Any columns
whose name is unique to one of the data frames will be deleted. min.presence states how often
an element should appear in order to be considered for analysis. The default setting of 1
will remove all non-present/ever-present elements. Increasing this will reduce the number of elements
considered but will increase the power of corellation tetss. Max presence is automatically set 
to the number of samples-min.presence for symmetry but can be changed.
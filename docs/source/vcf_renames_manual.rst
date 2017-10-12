vcf_rename_tags manual
======================

*vcf_rename_tags* is designed to update the INFO tags or Contig names in the VCF/BCF file.

::
   
   vcf_renames_tags  input_file
   Options:
     -list tags.txt     This file include two columns, they are old and new tags.
     -ncbi              Convert UCSC chromosome name to NCBI type.
     -ucsc              Convert NCBI chromosome name to UCSC type.
     -O <type>          Output type.[bzu]
     -o output_file     Output file.


 * -list parameter require a two column table file, each column seperate by a tab. First column is the old tags, and second column is the renamed tags.
   For example, "INFO/AF  dbSNP_AF" told the program to rename *AF* tag in the INFO field to *dbSNP_AF*, "CTG/chr1   1" told the program to rename
   contig name *chr1* to *1*. Please notice that capped 'CTG/' is mandantory for contig names, and if no such a capped information specifed, our program
   will only check the tags in the INFO fields.
 * -ncbi flag is used to convert UCSC contig names to NCBI names. For example, *chr1* will convert to *1* in the output file. Now our program only
   support human hg19 and hg38, for other released reference, should be specified with -list parameter.
 * -ucsc flag is used to convert NCBI contig names to UCSC names.
 * -O [bzu]  Output file format, *b* for compressed bcf, *z* for bgzipped vcf, *u* for uncompressed vcf.
 * -o specify output file.

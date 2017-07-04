# vcf2tsv manual

**vcf2tsv** was designed to convert to VCF/BCF file to tab-separated file with predefined *tags*.

About *tags*:  
* All the tags in the VCFs should be properly defined in the header.
* *tags* included CHROM, POS, ID, QUAL, Filter and all the tags in INFO and FORMAT fields.
* Besides the predefined tags, there are some buildin words for the format string,
	- **BED**, bed format of each position;
	- **SAMPLE**, sample column of this position, if this tag defined, samples will divided into lines with per sample per line;
	- **GENOTYPE**, gentype of this position, genotype should be *ref-alt, ref-ref, alt-alt*

Help message of **vcf2tsv**,
```
#vcf2tsv -h

About : Convert BCF/VCF to tsv file by selecting tags.
Usage:
	vcf2tsv -f string [Options] in.vcf.gz
Options:
	-f, --format        See man page for deatils.
	-s, --split         Split by [ALT].
	-p, --print-header  Print the header comment.
	-r, --skip-ref      Skip reference positions, when GT is "0/0"]`.
	-u, --skip-uncover  Skip uncover positions.
	-G, --no-GT         No check GT tag. For convert INFO only.
Website :
https://github.com/shiquan/vcfanno
```

About the parameters,
	- **f** : this is mandontory, selected tags must be defined with this parameter, and the columns of output file will consistant with these tags;
	- **s**: split flag, if ALT is set, one allele per record will be output. Multiple alternative alleles will be divided into several record for same position. *Notice*, tags with 'A' or 'R' will be splited by alleles;
	- **p**: print header descriptions for each column. These description lines will be comment with '#';

	- **r**: skip reference allele, if this parameter set, heterogentic position will only export alternative allele records


# vcf_rename_tags manual

```
vcf_rename_tags -list tags.txt -O <type> -o output_file input_file
```

About the parameters,
* **-list**, tags rename file, two columns consist of old tag name and new tag name;
* **-O** [b|u|v|z], output format, b for compressed bcf, v for uncompressed bcf, u for vcf, z for compressed vcf;
* **-o**, output file, standard output in default.

**vcf_rename_tags** support rename Contig and INFO tags for now, Contig tags should be inited with *CTG/*  in the tag names. 

## Example.
```
$ echo 'CTG/17\tchr17' > contig.list

$ cat contig.list
CTG/17	chr17

$ cat example/demo.vcf
##fileformat=VCFv4.2
##reference=file://17.fa.gz
##contig=<ID=17,length=81195210>
##ALT=<ID=X,Description="Represents allele(s) other than observed.">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	demo
17	41222826	.	A	C	.	.	.	GT	0/1
17	41223242	.	G	C	.	.	.	GT	0/1
17	41234451	.	G	A	.	.	.	GT	0/1
17	41258326	.	A	G	.	.	.	GT	0/1
17	41258504	.	A	C	.	.	.	GT	0/1

$ vcf_rename_tags -list contig.list example/demo.vcf
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##reference=file://17.fa.gz
##contig=<ID=chr17,length=81195210>
##ALT=<ID=X,Description="Represents allele(s) other than observed.">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	demo
chr17	41222826	.	A	C	.	.	.	GT	0/1
chr17	41223242	.	G	C	.	.	.	GT	0/1
chr17	41234451	.	G	A	.	.	.	GT	0/1
chr17	41258326	.	A	G	.	.	.	GT	0/1
chr17	41258504	.	A	C	.	.	.	GT	0/1

// Please notice the contig names have changed from `17` to `chr17`.
// You could also rename the INFO tags with predefined rename list.

```

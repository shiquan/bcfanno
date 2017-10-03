Getting Starting
================

This document will show you how to get up and running with bcfanno.
You will download the source code from GitHub in several minuters,
build your own configure file and annotate variants in the demo vcf.

build the program
------------------
Download source codes from github::

  $git clone https://github.com/shiquan/bcfanno.git
  $cd bcfanno
  $make

Following execute programs should be compiled after several miniters.

* **bcfanno** , core program to annotate genetic variants
* **tsv2vcf** (https://github.com/shiquan/bcfanno/blob/master/documents/tsv2vcf_manual.md) ,  generate VCF databases from tab-seperated file
* **vcf2tsv** (https://github.com/shiquan/bcfanno/blob/master/documents/vcf2tsv_manual.md), convert VCF file to tab-separated file with selected tags
* **vcf_rename_tags** (https://github.com/shiquan/bcfanno/blob/master/documents/vcf_rename_tags_manual.md), rename tags or contig names in the VCF file, usually used to format the databases

*Some other programs should also be install for ongoing test.*

* **`bcftools_`** 
* **`tabix_`** 



.. _bcftools:http://www.htslib.org/download/
.. _tabix:http://www.htslib.org/download/

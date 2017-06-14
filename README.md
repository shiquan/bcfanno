
VCFANNO
----------

VCFANNO was designed to annotate VCF/BCF files by using local and online databases. The early propose was design an efficient standalone program and put all the annotate informations in the seperated tags in VCF file.

Before you try to use VCFANNO, please make sure you know how to find the document of [VCF/BCF format](https://samtools.github.io/hts-specs/VCFv4.2.pdf).



**The design of vcfanno**

The VCF file are currently most used to store genomic variants and other related information. For other The INFO column of VCF/BCF file could specify all possible information related with variants in one position, and 
The core VCFANNO [htslib](http://htslib.org/)





**How to build the programs**







**Write configure file.**





**Generate databases for annotation.**



(1) Gene region databases.

(2) Allele frequency databases.

* G1000
* ExAC
* GomAD
* Local databases.

(3) Prediction databases.

* [How to generate dbNSFP and dbscsnv databases.](https://github.com/shiquan/vcfanno/blob/master/documents/database/dbNSFP.md)
* Condel

(4) Genotype and phenotype databases.

(5) Transcript related databases.

(6) SQL databases.

* How to connect online sql server.

  ​

**Convert annotated vcf file to other formats.**





**Interpret the annotations.**





**Bug report or suggestions**.



**The story behind bcfanno.**



**Reference**


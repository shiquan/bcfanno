Features
========

* **Fast : High Speed, Small Memory**
bcfanno was designed to annotate VCF/BCF files by using local and online databases. bcfanno wrote in pure C language, designed to cache the index structure instead of records in the memory. This design very small memory (usually <5M per dataset). Databases in VCF/BCF and bed format are supported, users could download and build open source and in house databases manually.

* **Flexible : Every File is a Database**
VCF is well supported and extensively used in a lot of projects. With the advances of VCF file, bcfanno was designed to put all the annotated data (dubbed *tags*) in the *INFO* and export annotated VCF. And the annotated VCF can be reused as database in next batch. For example, for one project, user could annotate all the reported genetic variants with preferred databases in the first batch, and the annotated VCF could be released as an ALL-IN-ONE database to be reused to annotate any further variant file in the same project but with very small memory and CPU resources.  For another example, to annotate a lot of VCFs, saying 1,000 WGS VCFs, merge all the VCFs into one VCF first and annotate the merged VCF in one batch, so for same genome position in different samples will be annotated just once. Compare with flat text files, manage and update VCFs could be much easier and more efficient with current stat softwares like bcftools and vcflib. bcfanno also supply two programs (**tsv2vcf** and **vcf2tsv**) to convert tab seperared text file to VCF and convert VCF to user friendly text file.

* **Accuracy : Realignment of transcript**

* **Extensible**



  
  
**Databases could be classified into three types**
-------------

bcfanno support three kinds of databases:

* *Allele Specific Databases (ASD)*, like allele frequency or any other allele related databases;
* *Region Specific Databases (RSD)*, like gene region, function regions or any other genome region related databases;
* *Transcript Specific Databases (TSD)*, like protein id or any other transcript related databases.


Some current stat-of-art databases listed below:
1. **ASD**: dbSNP, dbNSFP, dbSCSNV, dbExAC, HGVS nomenclature, GWAStrait, ...
2. **RSD**: OMIM, Geome Location, Gene, ...
3. **TSD**: Ensemble names


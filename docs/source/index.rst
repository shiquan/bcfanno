BCFANNO Documentation
===================================

Release note:
-------------

*release v1.0 :  first stable version, renamed vcfanno to bcfanno.*

**Introduction**
----------------
Flood of DNA data sequenced and the majority of the challange encounter with understanding these data is anaysis and interpret the genetic variants. Current-state-of-the-art bioinformatic tools like BWA, GATK, samtools accelerated scientific research and generated high quality genetic variants in different labs around the world. A number of annotation tools, like ANNOVAR, VAAST, SettleSeq, SNPeff, VEP, and InterVar, can predict the genetic variants affect transcript structure and combine more than 20 databases to interpret the pathogenic variants. However, to efficient annotate large scale data like genetic variants in WGS samples is still a challange. Cloud based tools like vcfanno, highly accelerate the annotation but required large clusters, huge computer memory and CPU, which is another bottleneck for small labs. Here I designed bcfanno, implement from personal PC to cloud-based server, a fast, flexible annotation tool to interpret genetic variants with very low computer resource. 

**Fast : High Speed, Small Memory**

bcfanno was designed to annotate VCF/BCF files by using local and online databases. bcfanno wrote in pure C language, designed to cache the index structure instead of records in the memory. This design very small memory (usually <5M per dataset). Databases in VCF/BCF and bed format are supported, users could download and build open source and in house databases manually.

**Flexible : Every File is a Database**

VCF is well supported and extensively used in a lot of projects. With the advances of VCF file, bcfanno was designed to put all the annotated data (dubbed *tags*) in the *INFO* and export annotated VCF. And the annotated VCF can be reused as database in next batch. For example, for one project, user could annotate all the reported genetic variants with preferred databases in the first batch, and the annotated VCF could be released as an ALL-IN-ONE database to be reused to annotate any further variant file in the same project but with very small memory and CPU resources.  For another example, to annotate a lot of VCFs, saying 1,000 WGS VCFs, merge all the VCFs into one VCF first and annotate the merged VCF in one batch, so for same genome position in different samples will be annotated just once. Compare with flat text files, manage and update VCFs could be much easier and more efficient with current stat softwares like bcftools and vcflib. bcfanno also supply two programs (**tsv2vcf** and **vcf2tsv**) to convert tab seperared text file to VCF and convert VCF to user friendly text file.



Bechmark
--------
bcfanno annotate a WGS vcf with specified 22 databases in 6 CPU hours. 

bcfanno usually do not parse the *FORMAT* of VCFs, and all the *tags* will be put into the INFO region, it is strongly suggested to merge multiple samples by `bcftools merge` before annotation.

**How to build the programs**
-----------------------------

::
   git clone https://github.com/shiquan/vcfanno.git
   cd vcfanno
   make


Following execute programs should be compiled after several miniters.

* **vcfanno** , core program to annotate genetic variants
* **tsv2vcf** (https://github.com/shiquan/vcfanno/blob/master/documents/tsv2vcf_manual.md) ,  generate VCF databases from tab-seperated file
* **vcf2tsv** (https://github.com/shiquan/vcfanno/blob/master/documents/vcf2tsv_manual.md), convert VCF file to tab-separated file with selected tags
* **vcf_rename_tags** (https://github.com/shiquan/vcfanno/blob/master/documents/vcf_rename_tags_manual.md), rename tags or contig names in the VCF file, usually used to format the databases

*Some other programs should also be install for ongoing test.*

* **bcftools** (http://www.htslib.org/download/) 
* **tabix** (http://www.htslib.org/download/)



**Databases**
-------------

bcfanno support THREE kinds of databases,

* *allele specific databases*, like allele frequency or any other allele related databases;
* *region specific databases*, like gene region, function regions or any other genome region related databases;
* *transcript specific databases*, like protein id or any other transcript related databases.

I think most of database could be classsified into above three kinds. Please email if these three don't define your databases properly.

Before you try to use any database, please try to classify the data and follow the below instruction to check or convert the database for bcfanno.


**Generate databases for bcfanno**
----------------------------------
Databases should be convert to VCF/BCF and BED-like region format. Good thing is the most databases were released in VCF or BED-like format, so we just need download them and do some updates for these kind of files, like dbsnp, EXAC and ClinVar etc. However, there are still some databases like dbNSFP were released in plain text format or other format, and we should convert them manually. For this section, we are only trying to build *clinvar* for getting start. All the details about build and convert databases could be find at section [Generate databases](https://github.com/shiquan/vcfanno/blob/master/Documentation/database/more_details.md).


*Instruction to build dbsnp for bcfanno:*

Step 1, enter the example directory and download clinvar in it.

::
   wget -c ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz`

Step 2, build index for clinvar.

::
   bcftools index clinvar.vcf.gz

Step 3, write configure file for bcfanno.

Create *clinvar.json*, and this file should be looked like

::
   {
        "vcfs":[
        	{
                	"file":"path to clinvar.vcf.gz",
                        "columns":"RS,CLNSIG",
                },
         ],
   }


Step 4, annotation. (If you are in the example directory now you can just run this command)
::
   ../vcfanno -c clinvar.json demo.vcf

The demo.vcf file would be annotated with clinvar databases. Try to compare the raw vcf and annotated vcf, see what's happened.


Write configure file.
---------------------

Configure file should be wrote in JSON format. I usually suggest my colleagues to edit the belowed copy of conifgure file and change the database path and tags accordingly.

Please notice that do not change the reserved keywords : *id*, *author*, *ref*, *hgvs*, *vcfs*, and *beds*.

.. code-block:: demo.json
   :language: json

**Convert annotated vcf file to other formats.**
-----------------------------------------------

**vcf2tsv** is a part of bcfanno package, convert selected tags from VCF/BCF to tab-seperated file.  For the usage of vcf2tsv please refer to *vcf2tsv manual* (https://github.com/shiquan/vcfanno/blob/master/documents/vcf2tsv_manual.md).

.. code-block:: demo.sh
   :language: bash


**Interpret the annotations.**
-------------------------------

+ For human genetic variants

The American College of Medical Genetics and Genomics (ACMG) supply a decision-tree roadmap and recommend using 28 criteria to help the clinical researcher to interpret genetic variants, however no computation approach could interpret the genetic variants directly, that's because gathering information for all the criteria is quite complicated and no specific algorithms for implementing this guidelines specified. See details, please refer to  **ACMG guideline** (https://www.acmg.net/docs/standards_guidelines_for_the_interpretation_of_sequence_variants.pdf).



Section  **VarType and HGVSnom** (https://github.com/shiquan/vcfanno/blob/master/Documentation/genetic_variant_types.md) introduce the genetic variant types and HGVS nomenclature.

Section **Interpret genetic variants** introduce how to use bcfanno and open distribute databases to interpret the pathogenic or regulatory variants.


+ For other species

bcfanno designed to annotate VCFs with suitable databases, not restrict to human variants. However the interpret rules may vary from different labs, there is no recommended strategy.



**Bug report or suggestions**.
------------------------------
Currently, you could report bugs from GitHub or email me directly. Please be kind to specify which exactly version you test in the report message.



Reference
----------


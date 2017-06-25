
BCFANNO
==========

*release v1.0 :  first stable version, renamed vcfanno to bcfanno.*



## **Introduction**

bcfanno was designed to annotate VCF/BCF files by using local and online databases. bcfanno support three kinds of databases which list below. And users could download and build open source or in house databases manually. Please refer to **database** section for the details.

Borrow the advances of VCF format, bcfanno was designed to put all the annotated data (dubbed *tags*) in the *INFO*, and the annotated VCF file could be convert to other formats or reused as a database. So before you try to use bcfanno, please make sure you know the structure of VCF file and how to find the document of [VCF/BCF format](https://samtools.github.io/hts-specs/VCFv4.2.pdf).



## **How to build the programs**

*prerequisite:*

* zlib

`git clone https://github.com/shiquan/vcfanno.git`

`cd vcfanno`

`make`

Following execute programs should be compiled after several miniters.

* ***vcfanno*** , core program to annotate genetic variants
* ***tsv2vcf*** ,  generate VCF databases from tab-seperated file
* ***vcf2tsv***, convert VCF file to tab-separated file with selected tags
* ***vcf_rename_tags***, rename tags or contig names in the VCF file, usually used to format the databases



*Some other programs should also be install for ongoing test.*

* [bcftools](http://www.htslib.org/download/) 

* [tabix](http://www.htslib.org/download/) (tabix is now a part of HTSlib, so download the htslib and complier it, and you will find tabix then)

  ​

## **Databases**

bcfanno support THREE kinds of databases,

* *allele specific databases*, like allele frequency or any other allele related databases
* *region specific databases*, like gene region, function regions or any other genome region related databases
* *transcript specific databases*, like protein id or any other transcript related databases

I think most of database could be classsified into above three kinds. Please email if these three don't define your databases properly.

Before you try to use any database, please try to classify the data and follow the below instruction to check or convert the database for bcfanno.



## **Generate databases for bcfanno**

Databases should be convert to VCF/BCF and BED-like region format. Good thing is the most databases were released in VCF or BED-like format, so we just need download them and do some updates for these kind of files, like dbsnp, EXAC and ClinVar etc. However, there are still some databases like dbNSFP were released in plain text format or other format, and we should convert them manually. For this section, we are only trying to build *clinvar* for getting start. All the details about build and convert databases could be find at [More details about databases](https://github.com/shiquan/vcfanno/blob/master/documents/database/more_details.md).



*Instruction to build dbsnp for bcfanno:*

Step 1, enter the example directory and download clinvar in it.

`wget -c ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz`

Step 2, build index for clinvar.

`bcftools index clinvar.vcf.gz`

Step 3, write configure file for bcfanno.

Create *clinvar.json*, and this file should be looked like

```

```

{

​	"vcfs":[

​		{

​		"file":"path to clinvar.vcf.gz",

​		"columns":"RS,CLNSIG",

​		},

​	],

}

```

```

Step 4, annotation. (If you are in the example directory now you can just run this command)

`../vcfanno -c clinvar.json demo.vcf`

The demo.vcf file would be annotated with clinvar databases. Try to compare the raw vcf and annotated vcf, see what was happened.

![](https://github.com/shiquan/vcfanno/blob/master/documents/database/raw_vcf.png)

After annotation.

![](https://github.com/shiquan/vcfanno/blob/master/documents/database/anno_vcf.png)



## **Write configure file.**

Configure file should be wrote in json format. Please remember we have some reserved keywords. Just copy the demo.json in the example directory and edit your own configure file from it.

```

```

{

"id":"configure ID and version",

"author":"author of this configure file",

"ref":"hg19",  // hg19 or hg38

"hgvs":{

​    "gene_data":"/opt/databases/refgene/hg19_refgene.tsv.gz",

​    "refseq":"/opt/databases/refgene/refMrna.fa.gz",

},

"vcfs":[    

],

"beds":[

],

}

```

```

## **Convert annotated vcf file to other formats.**





## **Interpret the annotations.**





## **Bug report or suggestions**.





## **Reference**


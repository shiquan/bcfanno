
BCFANNO
==========

*release v1.0 :  first stable version, renamed vcfanno to bcfanno.*



## **Introduction**

bcfanno was designed to annotate VCF/BCF files by using local and online databases. bcfanno support three kinds of databases which list below. And users could download and build open source or in house databases manually. Please refer to **database** section for the details.

Borrow the advances of VCF format, bcfanno was designed to put all the annotated data (dubbed *tags*) in the *INFO*, and the annotated VCF file could be convert to other formats or reused as a database. So before you try to use bcfanno, please make sure you know the structure of VCF file and how to find the document of [VCF/BCF format](https://samtools.github.io/hts-specs/VCFv4.2.pdf).



## Bechmark

bcfanno annotate a WGS vcf with specified 22 databases in 6 CPU hours. bcfanno usually do not parse the FORMAT of vcf files, and all the annotators (*tags*) will be put into the INFO region, so it is suggested to merge multiple samples by `bcftools merge` before annotation.





## **How to build the programs**

*prerequisite:*

* zlib

```

git clone https://github.com/shiquan/vcfanno.git
cd vcfanno
make
```

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
{
        "vcfs":[
        	{
                	"file":"path to clinvar.vcf.gz",
                        "columns":"RS,CLNSIG",
                },
         ],
}
```

Step 4, annotation. (If you are in the example directory now you can just run this command)

`../vcfanno -c clinvar.json demo.vcf`

The demo.vcf file would be annotated with clinvar databases. Try to compare the raw vcf and annotated vcf, see what was happened.

![](https://github.com/shiquan/vcfanno/blob/master/documents/database/raw_vcf.png)

After annotation.

![](https://github.com/shiquan/vcfanno/blob/master/documents/database/anno_vcf.png)



## **Write configure file.**

Configure file should be wrote in json format. Please remember we have some reserved keywords. Just copy the demo.json in the example directory and edit your own configure file from 

```json
{
        "id":"configure ID and version",
        "author":"author of this configure file",
        "ref":"hg19",  // hg19 or hg38
        "hgvs":{
           "gene_data":"/opt/databases/refgene/hg19_refgene.tsv.gz",
           "refseq":"/opt/databases/refgene/refMrna.fa.gz",
        },
        "vcfs":[ 
          {
            "file":"path to clinvar.vcf.gz",
       	    "columns":"RS,CLNSIG",
          },
          {
     	    "file":"path to vcf database",
            "columns":"tags",
          },
        ],
        "beds":[
          {
            "file":"path to BED-like database",
            "columns":"tags",
          },          
        ],
}
```

## **Convert annotated vcf file to other formats.**

**vcf2tsv** is a part of bcfanno package, convert selected tags in VCF/BCF to tab-seperated file.  For the usage of vcf2tsv please refer to [vcf2tsv manual]().

```
vcf2tsv -f BED,REF,ALT,GT,SAMPLE,Gene,HGVSnom,ExonIntron,VarType,HGMD_tag example/demo_anno.vcf

// results 
#CHROM	START	END	REF	ALT	GT	SAMPLE	Gene	HGVSnom	ExonIntron	VarType	HGMD_tag
17	41222825	41222826	A	C	0/1	demo	BRCA1|BRCA1|BRCA1|BRCA1|BRCA1|BRCA1	NM_007294.3:c.4987-3114A>C|NM_007297.3:c.4846-3114A>C|NM_007298.3:c.1675-3114A>C|NM_007299.3:c.1675-3114A>C|NM_007300.3:c.5050-3114A>C|NR_027676.1:n.5123-3114A>C	I15|I14|I14|I15|I16|I15intron|intron|intron|intron|intron|intron	.
17	41223241	41223242	G	C	0/1	demo	BRCA1|BRCA1|BRCA1|BRCA1|BRCA1|BRCA1	NM_007294.3:c.4689C>G(p.Tyr1563Stop/p.Y1563X)|NM_007297.3:c.4548C>G(p.Tyr1516Stop/p.Y1516X)|NM_007298.3:c.1377C>G(p.Tyr459Stop/p.Y459X)|NM_007299.3:c.1377C>G(p.Tyr459Stop/p.Y459X)|NM_007300.3:c.4752C>G(p.Tyr1584Stop/p.Y1584X)|NR_027676.1:n.4825C>G	E15/C14|E14/C12|E14/C14|E15/C14|E16/C15|E15/C15	nonsense|nonsense|nonsense|nonsense|nonsense|noncoding	DM
17	41234450	41234451	G	A	0/1	demo	BRCA1|BRCA1|BRCA1|BRCA1|BRCA1|BRCA1	NM_007294.3:c.4327C>T(p.Arg1443Stop/p.R1443X)|NM_007297.3:c.4186C>T(p.Arg1396Stop/p.R1396X)|NM_007298.3:c.1018C>T(p.Arg340Stop/p.R340X)|NM_007299.3:c.1018C>T(p.Arg340Stop/p.R340X)|NM_007300.3:c.4327C>T(p.Arg1443Stop/p.R1443X)|NR_027676.1:n.4463C>T	E12/C11|E11/C9|E11/C11|E12/C11|E12/C11|E12/C12	nonsense|nonsense|nonsense|nonsense|nonsense|noncoding	DM
17	41258325	41258326	A	G	0/1	demo	BRCA1|BRCA1|BRCA1|BRCA1|BRCA1|BRCA1	NM_007294.3:c.213-1353A>G|NM_007297.3:c.72-1353A>G|NM_007298.3:c.213-1353A>G|NM_007299.3:c.213-1353A>G|NM_007300.3:c.213-1353A>G|NR_027676.1:n.352-1353A>G	I4|I3|I3|I4|I4|I4	intron|intron|intron|intron|intron|intron	.
17	41258503	41258504	A	C	0/1	demo	BRCA1|BRCA1|BRCA1|BRCA1|BRCA1|BRCA1	NM_007294.3:c.181T>G(p.Cys61Gly/p.C61G)|NM_007297.3:c.40T>G(p.Cys14Gly/p.C14G)|NM_007298.3:c.181T>G(p.Cys61Gly/p.C61G)|NM_007299.3:c.181T>G(p.Cys61Gly/p.C61G)|NM_007300.3:c.181T>G(p.Cys61Gly/p.C61G)|NR_027676.1:n.342T>G	E4/C3|E3/C1|E3/C3|E4/C3|E4/C3|E4/C4	missense|missense|missense|missense|missense|noncoding	DM
```



## **Interpret the annotations.**





## **Bug report or suggestions**.





## **Reference**


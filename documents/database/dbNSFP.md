
Format dbNSFP and dbscSNV databases for vcfanno.
==================================

*<u>Before download and format dbNSFP and dbscSNV databases for your project, please make sure you have right to use them in your research or commerical activities. This is only a protocol to convert the standard datasets into VCF/BCF format, and I don't buy any licence to republish both of them.</u>*

 The original dbNSFP and dbscsnv databases should be download from *https://sites.google.com/site/jpopgen/dbNSFP* or other mirrors. Please notice a lot of versions for both of databases should be found at this site but I recommend  you to download the most updated version. However, this is not always true, especially when you use a old version of human reference genome in the upstream analysis, like hg19. *Find the most suitable version compared with your reference is suggested.*



Standard dbNSFP and dbscSNV databases are release in gziped tab seperated text format, here we strong suggest to convert all these kind of datasets into BCF format for performance and accuracy. Please try to download the raw datasets and README file from same place, and please follow the listed steps to generate your own dbNSFP and dbscSNV BCF databases.



**Prerequisite:**

1. Free internet connection to access the homepage. For users from mainland of China, updated datasets could be found at *ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/*.
2. **bcftools**  (could be download from *http://www.htslib.org/*)
3. **tsv2vcf**    (a sub program of vcfanno package, you will find this program after make package.)
4. any perfered text editor
5. (optional) **vcf_rename_tags** (a sub program of vcfanno package)




**Protocol (for hg38):**

1. Download dbnsfp from homepage, and unzipped it.

```
$ wget -c ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv3.4c.zip

$ unzip dbNSFPv3.4c.zip

$ cd dbNSFPv3.4c

$ ls

// belowed datasets should be found at this directory
```
   ![](https://github.com/shiquan/vcfanno/blob/master/documents/database/dbNSFP_files.png)


2. Download the related README file and generate INFO descriptions file manually, all the description information could be find from README file.

   (1) Open README file (this file could be found at dbNSFP package).

   (2) Generate the INFO descriptions for our VCF/BCF database. Please make sure you know the format of VCF header clearly. If no, please refer to *http://samtools.github.io/hts-specs/VCFv4.3.pdf* for the technical knowledge and copy my pre-defined demo INFO description file (*https://github.com/shiquan/vcfanno/blob/master/documents/database/dbNSFP_header.txt*) for your sake.

   (3)  Convert each chromosome dataset into BCF files. 

   * Check the format of each dataset. First line of dataset (dubbed *header*) should be comment and the column number of *header* should be consistent with other lines (dubbed *body*); the *header* of the row should be consistent with VCF INFO *tag*.

   ![](https://github.com/shiquan/vcfanno/blob/master/documents/database/dbNSFP_header.png)

   * Convert tablet to VCF.
```

Usage : tsv2vcf -header|-h header.txt -r reference.fa [-force -pos column -O z -o out.vcf.gz] in.tsv.gz
        -header, -h     header file
        -r              reference file
        -pos            position column, if set will skip pos,start,end column in the title
        -start          start position column, inconsistance with -pos
        -end            end position column, only set if need add a END tag in the INFO
        -chr            chr column, if set will skip first column in the title
        -force          if reference seq and fasta file are inconsistent, just give a warning
        -rename         chromosome rename file

Homepage: https://github.com/shiquan/vcfanno

```

**Note :**  the *-r* parameter is mandatory, because program will check each reference base in the datasets, if there are some inconsistance for genetic bases between dbNSFP (or other databases) and human genome reference, you must figure out how it comes and fix it by using right reference or change the dataset manually (I do *not* explicitly recommend change any database manually unless you know what exactly you do and it is suggested to report bugs to the authors).

More details about tsv2vcf, please refer to \href[tsv2vcf manual](https://github.com/shiquan/vcfanno/blob/master/documents/tsv2vcf_manual.md).

   ![](https://github.com/shiquan/vcfanno/blob/master/documents/database/tsv2vcf_con.png)

   The convert command could be:

```
tsv2vcf -chr 1 -pos 2 -h dbnsfp_header.txt -r Homo_sapiens.GRCh38.dna.toplevel.fa dbNSFP3.4c_variant.chr1 | bcftools view -O b -o chr1.bcf
```


   (4)  `bcftools concat` all the chromosomes BCF files into one big BCF, and index it.

```
$ bcftools concat chr1.bcf chr2.bcf chr3.bcf chr4.bcf chr5.bcf chr6.bcf chr7.bcf chr8.bcf chr9.bcf chr10.bcf chr11.bcf chr12.bcf chr13.bcf chr14.bcf chr15.bcf chr16.bcf chr17.bcf chr18.bcf chr19.bcf chr20.bcf chr21.bcf chr22.bcf chrX.bcf chrY.bcf chrM.bcf -O b -o dbNSFPv3.4c.bcf

$ bcftools index dbNSFPv3.4c.bcf 
```

   (5)  If you need rename chromosomes, like if you download UCSC reference and the chromosome name is different in dbNSFP, use vcf_rename_tags rename the chromosome names.

```
$ vcf_rename_tags -list contig.txt dbNSFPv3.4c.bcf -O b -o dbNSFPv3.4c.renames.bcf

$ bcftools index dbNSFPv3.4c.renames.bcf

```
   *Note:* contig.txt could be found at https://github.com/shiquan/vcfanno/blob/master/documents/database/contig.txt

   ​

3. Test database with vcfanno.

   * Edit your configure file and add new database.
```
"vcfs": [
	{
		"columns":"SIFT_score,Polyphen2_HDIV_score,Polyphen2_HVAR_score,LRT_score,LRT_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_pred,FATHMM_pred,PROVEAN_pred,CADD_raw,fathmm-MKL_coding_score,MetaSVM_score,MetaLR_score,GERP++_RS,SiPhy_29way_logOdds",
	"file":"*<u>path to bcf database</u>*",
	},
]
```

   * Debug.




***FAQ***

1. How to convert dbscsnv database?

   Follow above protocol to download dbscsnv database and build the description file and convert the database by using same programs then.

2. If I use hg19 (GRCh37), which dbNSFP version should I download?

   My suggestion is dbNSFP **v2.9.1**, but you can always download the most updated version and convert the inconsistent reference bases manually.

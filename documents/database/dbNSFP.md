
Format dbNSFP and dbscSNV databases for vcfanno.
==================================

*Before download and format dbNSFP and dbscSNV databases for your project, please make sure you have right to use them in your research or commerical activities. This is only a protocol to convert the standard datasets into VCF/BCF format, and I don't buy any licence to republish both of them.*

 The original dbNSFP and dbscsnv databases should be download from *https://sites.google.com/site/jpopgen/dbNSFP* or other mirrors. Please notice a lot of versions for both of databases should be found at this site but I recommend  you to download the most updated version. However, this is not always true, especially when you use a old version of human reference genome in the upstream analysis, like hg19. *Find the most suitable version compared with your reference is suggested.*



Standard dbNSFP and dbscSNV databases are release in gziped tab seperated text format, here we strong suggest to convert all these kind of datasets into BCF format for performance and accuracy. Please try to download the raw datasets and README file from same place, and please follow the listed steps to generate your own dbNSFP and dbscSNV BCF databases.



**Prerequisite:**

1. Free internet connection to access the homepage. For users from mainland of China, updated datasets could be found at *ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/*.
2. **bcftools**  (could be download from *http://www.htslib.org/*)
3. **tsv2vcf**    (a part of vcfanno package, you will find this program after make package.)
4. any perfered text editor



**Protocol:**

1. Download dbnsfp from homepage, and unzipped it.

   `wget -c ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv3.4c.zip`

   `unzip dbNSFPv3.4c.zip ` 

   `cd dbNSFPv3.4c`

   `ls`

   `# belowed datasets should be found at this directory`

   ![](https://github.com/shiquan/vcfanno/blob/master/documents/database/dbNSFP_files.png)

2. Download the related README file and generate header file manually, all the description information could be find from README file.

   (1) Open this site in your browser: *https://drive.google.com/file/d/0B60wROKy6OqcSWJfRk80Q1pNU1E/view*

   (2) Copy or save related information into one local text file. (Please notice that the reason we keep this README file is to generate the header file for our VCF/BCF database, so just read this file and see what information you really care about.)

   (3) Generate the header file for our VCF/BCF database. Please make sure you know the format of VCF header clearly. If no, please refer to *http://samtools.github.io/hts-specs/VCFv4.3.pdf* for the technical knowledge and copy my pre-defined demo header (*https://github.com/shiquan/vcfanno/blob/master/documents/demo_header.vcf*) for your sake.

   (4)  Convert each chromosome dataset into BCF files. 

   ​	a. check the format of each chromosome dataset

   ​	b. convert	

       Usage : tsv2vcf -header|-h header.txt -r reference.fa [-force -pos column -O z -o out.vcf.gz] in.tsv.gz
           -header, -h     header file
           -r              reference file
           -pos            position column, if set will skip pos,start,end column in the title
           -start          start position column, inconsistance with -pos
           -end            end position column, only set if need add a END tag in the INFO
           -chr            chr column, if set will skip first column in the title
           -force          if reference seq and fasta file are inconsistent, just give a warning
           -rename         chromosome rename file

   ​

   (5)  `bcftools concat` all the chromosomes BCF files into one big BCF, and index it.

   ​

3. Generate BCF file with header.txt and plain text datasets.

4. Test database with vcfanno.

5. Debug.



***FAQ***

1. ​




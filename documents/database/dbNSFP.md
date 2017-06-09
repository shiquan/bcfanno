
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



**Protocol**

1. Download dbnsfp from homepage, and unzipped it.

   `wget -c ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv3.4c.zip`

   `unzip dbNSFPv3.4c.zip ` 

   `cd dbNSFPv3.4c`

   `ls`

   `#belowed datasets should be found at this directory`

   ​

2. Download the related README file and open it.

   (1) Open this site in your browser: *https://drive.google.com/file/d/0B60wROKy6OqcSWJfRk80Q1pNU1E/view*

   (2) Copy or save related information into one local text file. (Please notice that the reason we keep this README file is to generate the header file for our VCF/BCF database, so just read this file and see what information you really care about.)

   (3) Generate the header file for our VCF/BCF database. Please make sure you know the format of VCF header clearly. If no, please refer to *http://samtools.github.io/hts-specs/VCFv4.3.pdf* for the technical knowledge and copy my predefined demo header () for your sake.

   (4) 

3. Generate header file manually, all the description information could be find from README file.

4. Generate BCF file with header.txt and plain text datasets.

5. Test database with vcfanno.

6. Debug.



***FAQ***

1. ​




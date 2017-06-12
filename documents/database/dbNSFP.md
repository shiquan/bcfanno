
Format dbNSFP and dbscSNV databases for vcfanno.
==================================

*<u>Before download and format dbNSFP and dbscSNV databases for your project, please make sure you have right to use them in your research or commerical activities. This is only a protocol to convert the standard datasets into VCF/BCF format, and I don't buy any licence to republish both of them.</u>*

 The original dbNSFP and dbscsnv databases should be download from *https://sites.google.com/site/jpopgen/dbNSFP* or other mirrors. Please notice a lot of versions for both of databases should be found at this site but I recommend  you to download the most updated version. However, this is not always true, especially when you use a old version of human reference genome in the upstream analysis, like hg19. *Find the most suitable version compared with your reference is suggested.*



Standard dbNSFP and dbscSNV databases are release in gziped tab seperated text format, here we strong suggest to convert all these kind of datasets into BCF format for performance and accuracy. Please try to download the raw datasets and README file from same place, and please follow the listed steps to generate your own dbNSFP and dbscSNV BCF databases.



**Prerequisite:**

1. Free internet connection to access the homepage. For users from mainland of China, updated datasets could be found at *ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/*.
2. **bcftools**  (could be download from *http://www.htslib.org/*)
3. **tsv2vcf**    (a part of vcfanno package, you will find this program after make package.)
4. any perfered text editor




**Protocol (for hg38):**

1. Download dbnsfp from homepage, and unzipped it.

   `wget -c ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv3.4c.zip`

   `unzip dbNSFPv3.4c.zip ` 

   `cd dbNSFPv3.4c`

   `ls`

   `# belowed datasets should be found at this directory`

   ![](https://github.com/shiquan/vcfanno/blob/master/documents/database/dbNSFP_files.png)

2. Download the related README file and generate INFO descriptions file manually, all the description information could be find from README file.

   (1) Open README file (this file could be found at dbNSFP package).

   (2) Generate the INFO descriptions for our VCF/BCF database. Please make sure you know the format of VCF header clearly. If no, please refer to *http://samtools.github.io/hts-specs/VCFv4.3.pdf* for the technical knowledge and copy my pre-defined demo INFO description file (*https://github.com/shiquan/vcfanno/blob/master/documents/demo_header.vcf*) for your sake.

   ​

   (3)  Convert each chromosome dataset into BCF files. 

   * Check the format of each dataset. First line of dataset (dubbed *header*) should be comment and the column number of *header* should be consistent with other lines (dubbed *body*); the *header* of the row should be consistent with VCF INFO *tag*.

   ![](https://github.com/shiquan/vcfanno/blob/master/documents/database/dbNSFP_header.png)

   * Convert tablet to VCF.

   ​

   **Note :**  the *-r* parameter is mandatory, because program will check each reference base in the datasets, if there are some inconsistance for genetic bases between dbNSFP (or other databases) and human genome reference, you must figure out how it comes and fix it by using right reference or change the dataset manually (I do *not* explicitly recommend change any database manually unless you know what exactly you do and bugs report to the author is suggested).

   ​

   (4)  `bcftools concat` all the chromosomes BCF files into one big BCF, and index it.

   ​

3. Test database with vcfanno.

4. Debug.



***FAQ***

1. ​




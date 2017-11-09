Preparing Clinvar database for bcfanno
===============================

Current-state-of-the-art bioinformatic tools like BWA, GATK, samtools accelerated scientific research and generated high quality genetic variants in different labs around the world. VCF is the most well supported format to store and

Databases should be convert to VCF/BCF and BED-like region format. Good thing is the most databases were released in VCF or BED-like format, so we just need download them and do some updates for these kind of files, like dbsnp, EXAC and ClinVar etc. However, there are still some databases like dbNSFP were released in plain text format or other format, and we should convert them manually. For this section, we are only trying to build *clinvar* for getting start. All the details about build and convert databases could be find at section [Generate databases](https://github.com/shiquan/vcfanno/blob/master/Documentation/database/more_details.md).


*Instruction to build dbsnp for bcfanno:*

Step 1, enter the example directory and download clinvar in it.

::

   wget -c ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz 

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

   ../bcfanno -c clinvar.json demo.vcf

   
The demo.vcf file would be annotated with clinvar. Try to compare the raw vcf and annotated vcf, new tags will be added in the INFO.


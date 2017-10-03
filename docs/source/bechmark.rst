Bechmark
--------
bcfanno annotate a WGS vcf with specified 22 databases in 6 CPU hours. 

bcfanno usually do not parse the *FORMAT* of VCFs, and all the *tags* will be put into the INFO region, it is strongly suggested to merge multiple samples by `bcftools merge` before annotation.


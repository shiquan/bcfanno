Introduction
===========

Flood of DNA data sequenced and the majority of the challange encounter with understanding these data is anaysis and interpret the genetic variants. 

A number of annotation tools, like ANNOVAR, VAAST, SettleSeq, SNPeff, VEP, and InterVar, can predict the genetic variants affect transcript structure and combine more than 20 databases to interpret the pathogenic variants. However, to efficient annotate large scale data like genetic variants in WGS samples is still a challange. Cloud based tools like vcfanno, highly accelerate the annotation but required large clusters, huge computer memory and CPU, which is another bottleneck for small labs.

Even though protein biochemistry has been used to characterize missense and nonsense coding mutations that most often underlie monogenic traits, the frequency with which loss-of-function mutations and rare coding variants are being discovered in healthy individuals suggests our understanding is far from complete. 

Here I designed bcfanno, implement from personal PC to cloud-based server, a fast, flexible annotation tool to interpret genetic variants with very low computer resource. 


**Reference**
1. Wang, K., Li, M., & Hakonarson, H. (2010). ANNOVAR: Functional annotation of genetic variants from high-throughput sequencing data. Nucleic Acids Research, 38(16). http://doi.org/10.1093/nar/gkq603
2. Hu, H., Huff, C. D., Moore, B., Flygare, S., Reese, M. G., & Yandell, M. (2013). VAAST 2.0: Improved variant classification and disease-gene identification using a conservation-controlled amino acid substitution matrix. Genetic Epidemiology, 37(6), 622–634. http://doi.org/10.1002/gepi.21743
3. Li, Q., & Wang, K. (2017). InterVar: Clinical Interpretation of Genetic Variants by the 2015 ACMG-AMP Guidelines. American Journal of Human Genetics, 100(2), 267–280. http://doi.org/10.1016/j.ajhg.2017.01.004
4. McLaren, W., Gil, L., Hunt, S. E., Riat, H. S., Ritchie, G. R. S., Thormann, A., … Cunningham, F. (2016). The Ensembl Variant Effect Predictor. Genome Biology, 17(1), 122. http://doi.org/10.1186/s13059-016-0974-4
5. Pedersen, B. S., Layer, R. M., & Quinlan, A. R. (2016). Vcfanno: fast, flexible annotation of genetic variants. Genome Biology, 17(1), 118. http://doi.org/10.1186/s13059-016-0973-5
   

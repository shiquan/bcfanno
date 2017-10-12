Data format convertor
====================

**vcf2tsv** is a part of bcfanno package, convert selected tags from VCF/BCF to tab-seperated file.  For the usage of vcf2tsv please refer to *vcf2tsv manual* (https://github.com/shiquan/vcfanno/blob/master/documents/vcf2tsv_manual.md).


::
                                
   vcf2tsv -f BED,REF,ALT,GT,SAMPLE,Gene,HGVSnom,ExonIntron,VarType,HGMD_tag example/demo_anno.vcf
   // results 
   #CHROM	START	END	REF	ALT	GT	SAMPLE	Gene	HGVSnom	ExonIntron	VarType	HGMD_tag
   17	41222825	41222826	A	C	0/1	demo	BRCA1|BRCA1|BRCA1|BRCA1|BRCA1|BRCA1	NM_007294.3:c.4987-3114A>C|NM_007297.3:c.4846-3114A>C|NM_007298.3:c.1675-3114A>C|NM_007299.3:c.1675-3114A>C|NM_007300.3:c.5050-3114A>C|NR_027676.1:n.5123-3114A>C	I15|I14|I14|I15|I16|I15intron|intron|intron|intron|intron|intron	.
   17	41223241	41223242	G	C	0/1	demo	BRCA1|BRCA1|BRCA1|BRCA1|BRCA1|BRCA1	NM_007294.3:c.4689C>G(p.Tyr1563Stop/p.Y1563X)|NM_007297.3:c.4548C>G(p.Tyr1516Stop/p.Y1516X)|NM_007298.3:c.1377C>G(p.Tyr459Stop/p.Y459X)|NM_007299.3:c.1377C>G(p.Tyr459Stop/p.Y459X)|NM_007300.3:c.4752C>G(p.Tyr1584Stop/p.Y1584X)|NR_027676.1:n.4825C>G	E15/C14|E14/C12|E14/C14|E15/C14|E16/C15|E15/C15	nonsense|nonsense|nonsense|nonsense|nonsense|noncoding	DM
   17	41234450	41234451	G	A	0/1	demo	BRCA1|BRCA1|BRCA1|BRCA1|BRCA1|BRCA1	NM_007294.3:c.4327C>T(p.Arg1443Stop/p.R1443X)|NM_007297.3:c.4186C>T(p.Arg1396Stop/p.R1396X)|NM_007298.3:c.1018C>T(p.Arg340Stop/p.R340X)|NM_007299.3:c.1018C>T(p.Arg340Stop/p.R340X)|NM_007300.3:c.4327C>T(p.Arg1443Stop/p.R1443X)|NR_027676.1:n.4463C>T	E12/C11|E11/C9|E11/C11|E12/C11|E12/C11|E12/C12	nonsense|nonsense|nonsense|nonsense|nonsense|noncoding	DM
   17	41258325	41258326	A	G	0/1	demo	BRCA1|BRCA1|BRCA1|BRCA1|BRCA1|BRCA1	NM_007294.3:c.213-1353A>G|NM_007297.3:c.72-1353A>G|NM_007298.3:c.213-1353A>G|NM_007299.3:c.213-1353A>G|NM_007300.3:c.213-1353A>G|NR_027676.1:n.352-1353A>G	I4|I3|I3|I4|I4|I4	intron|intron|intron|intron|intron|intron	.
   17	41258503	41258504	A	C	0/1	demo	BRCA1|BRCA1|BRCA1|BRCA1|BRCA1|BRCA1	NM_007294.3:c.181T>G(p.Cys61Gly/p.C61G)|NM_007297.3:c.40T>G(p.Cys14Gly/p.C14G)|NM_007298.3:c.181T>G(p.Cys61Gly/p.C61G)|NM_007299.3:c.181T>G(p.Cys61Gly/p.C61G)|NM_007300.3:c.181T>G(p.Cys61Gly/p.C61G)|NR_027676.1:n.342T>G	E4/C3|E3/C1|E3/C3|E4/C3|E4/C3|E4/C4	missense|missense|missense|missense|missense|noncoding	DM



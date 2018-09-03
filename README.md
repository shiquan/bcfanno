
Getting started with BCFANNO
================

bcfanno was designed to annotate genetic variants with various biological databases and predict the variant types and HGVS names quickly. 
Here are some easy instructions to get you up and running toy data.

## Installing bcfanno

Download the source code or binary distribution from https://github.com/shiquan/bcfanno/releases

Or cloning from github, get the most updated version :
```
git clone https://github.com/shiquan/bcfanno
cd bcfanno && make
```

## Test toy data

Run bcfanno with toy data :

```
# stdout output
./bcfanno -c toy.json example/toy.vcf.gz

# Or save results to a compressed VCF :
./bcfanno -c toy.json example/toy.vcf.gz -O z -o results.vcf.gz

# Or save as bcf :
./bcfanno -c toy.json example/toy.vcf.gz -O b -o results.bcf
```

## Select annotations from result
**vcf2tsv** is an additional program in bcfanno package, use to convert VCF/BCF file to user-friendly tab-seperated-variants file.
```
./bcfanno -c toy.json example/toy.vcf.gz -q | ./vcf2tsv -f CHROM,POS,CytoBand,REF,ALT,GT,SAMPLE,RS,MolecularConsequence,Gene,HGVSnom,ExonIntron,AAlength,HGMD_disease,HGMD_tag,HGMD_pmid

## Another usage
./bcfanno -c toy.json example/toy.vcf.gz -q | ./vcf2tsv -f BED,CytoBand,TGT,SAMPLE,RS,MolecularConsequence,Gene,HGVSnom,ExonIntron,AAlength,HGMD_disease,HGMD_tag,HGMD_pmid
```

## (Optional) View annotations with Microsoft excel
This step is optional, need to install my another program [tsv2excel](https://github.com/shiquan/tsv2excel) first.
```
./bcfanno -c toy.json example/toy.vcf.gz -q | ./vcf2tsv -f BED,CytoBand,TGT,SAMPLE,RS,MolecularConsequence,Gene,HGVSnom,ExonIntron,AAlength,HGMD_disease,HGMD_tag,HGMD_pmid | tsv2excel -o toy.xlsx
```
<img src="https://raw.githubusercontent.com/shiquan/bcfanno/master/Documentation/toy_anno.png" width="805" height="483" />

## Additional programs included in bcfanno package

Beside the core program ***bcfanno***, belowed programs will also be generated in the package.

* [***tsv2vcf***]() ,  generate VCF databases from tab-seperated file
* [***vcf2tsv***](), convert VCF file to tab-separated file with selected tags
* [***vcf_rename_tags***](), rename tags or contig names in the VCF file, usually used to format the databases
* ~~[***GenePredExtGen***]() Generate genepredext format with genome annotation and reference databases.~~


## Download databases from server
Please download the GEA databases from
* http://public.genomics.org.cn/bcfanno/databases/
* [mirror] https://pan.genomics.cn/ucdisk/s/UbeeIb

** A lot of databases could be download from different institute freely, we are not plan to redistribute these datasets (And of course, some license also require us not to do that). But if you have the database problem with bcfanno, please feel free to give me an email so that we can give you possible assistance. **

## Bug report or suggestions

Kindly report bugs and suggestions through github.


## LICENSE
The full package of bcfanno is distributed by MIT License, copyright 2016-2018 BGI Research.

Belowed package or source code used in bcfanno copyrighted by other institution.
- [htslib1.6](www.htslib.org)  The MIT/Expat License, Copyright (C) 2012-2014 Genome Research Ltd.
- thread_pool.[ch] The MIT/Expat License, Copyright (c) 2013-2017 Genome Research Ltd.

## How to cite bcfanno

We are trying to publish a paper to describe bcfanno, before that please cite https://github.com/shiquan/bcfanno in your work.


## Reference
1. [HGVS nomenclature](http://varnomen.hgvs.org/)
2. [VCF format](http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41/)

# tsv2vcf manual

**tsv2vcf** was designed to convert tab seperated file to VCF file.


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




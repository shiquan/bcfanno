Preparing annotation files for *Homo Sapiens*
===========================================

1. Go to NCBI's ftp ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/ , there will be one folder in this directory for each genome reference release. Here we use the GCF_000001405.36_GRCh38.p10 for demonstration.

::

   # Download GFF annotation file.
   wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.36_GRCh38.p10/GCF_000001405.36_GRCh38.p10_genomic.gff.gz

   # Download RNA reference sequence.
   wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.36_GRCh38.p10/GCF_000001405.36_GRCh38.p10_rna.fna.gz

   # Download genome reference sequence.
   wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.36_GRCh38.p10/GCF_000001405.36_GRCh38.p10_genomic.fna.gz

   # Decompress all files
   gzip -d GCF_000001405.36_GRCh38.p10_genomic.gff.gz
   gzip -d GCF_000001405.36_GRCh38.p10_rna.fna.gz
   gzip -d GCF_000001405.36_GRCh38.p10_genomic.fna.gz

   # Index fasta files using `samtools faidx`
   samtools faidx GCF_000001405.36_GRCh38.p10_rna.fna
   samtools faidx GCF_000001405.36_GRCh38.p10_genomic.fna

   # Convert GFF to GenePred file.
   gff3ToGenePred -rnaNameAttr=transcript_id -geneNameAttr=gene -useName GCF_000001405.36_GRCh38.p10_genomic.gff.gz GCF_000001405.36_GRCh38.p10_genomic.Genepred

   # Generated GenepredPlus file.
   comp_ref_trans -data GCF_000001405.36_GRCh38.p10_genomic.Genepred -rna GCF_000001405.36_GRCh38.p10_rna.fna -ref GCF_000001405.36_GRCh38.p10_genomic.fna -format refgene | sork -k2,2 -k4,4n -k5,5n | bgzip -c > GCF_000001405.36_GRCh38.p10_genomic.GenepredPlus.gz

   # Index GenepredPlus file using tabix.
   tabix -s 2 -b 4 -e 5 GCF_000001405.36_GRCh38.p10_genomic.GenepredPlus.gz

   # Generate configure file.
   {
        "hgvs": {
             "gene_data":"./GCF_000001405.36_GRCh38.p10_genomic.GenepredPlus.gz",
              "refseq":"./GCF_000001405.36_GRCh38.p10_rna.fna",
        },
   }


 
**Frequency asked questions:**

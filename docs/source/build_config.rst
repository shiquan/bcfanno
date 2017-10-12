
Write configure file.
====================

Configure file should be wrote in JSON format. I usually suggest my colleagues to edit the belowed copy of conifgure file and change the database path and tags accordingly.

Please notice that do not change the reserved keywords : *id*, *author*, *ref*, *hgvs*, *vcfs*, and *beds*.

::
   
   {
        "id":"configure ID and version", // optional
        "author":"author of this configure file", // optional
        "ref":"path to reference",  // fasta format and be indexed with samtools faidx, FLKSEQ tag will be generated if ref is set
        "hgvs":{
           "gene_data":"/opt/databases/refgene/hg19_refgene.tsv.gz",
           "refseq":"/opt/databases/refgene/refMrna.fa.gz",
          // "trans_list":"path to transcript list",  // this is optional
          // "gene_list":"path to gene list", // this is optional
        },
        "vcfs":[
          {
            "file":"path to clinvar.vcf.gz",
            "columns":"RS,CLNSIG",
          },
          {
            "file":"path to vcf database",
            "columns":"tags",
          },
        ],
        "beds":[
          {
            "file":"path to BED-like database",
            "columns":"tags",
          },
        ],
   }



/* This library is designed to generate and wrap HGVS nomenclature.
 * See http://www.hgvs.org/mutnomen/ for more information about HGVS nimenclature.
 *
 * Copyright (c) 2016  shiquan@link.cuhk.edu.hk
 *
 * This library is free to share, copy and modify for any people or institution.
 * Please report bugs or improvement to me.
 */
#ifndef HGVS_HEADER
#define HGVS_HEADER

/*splicing consensus regions*/

#define SPLIT5_UPSTREAM   3
#define SPLIT5_DOWNSTREAM 8
#define SPLIT3_UPSTREAM   12
#define SPLIT3_DOWNSTREAM 2
#define SPLITSITE_RANG 3

/* refgene format praser */
#define REFGENE_PRASE_BIN    1
#define REFGENE_PRASE_NAME1  (1<<1)
#define REFGENE_PRASE_REG    (1<<2)
#define REFGENE_PRASE_EXONS  (1<<3 | 1<<2)
#define REFGENE_PRASE_NAME2  (1<<4 | 1<<1)
#define REFGENE_PRASE_SUFFIX (1<<5)
#define REFGENE_PRASE_ALL   (REFGENE_PRASE_BIN | REFGENE_PRASE_EXONS | REFGENE_PRASE_NAME2 | REFGENE_PRASE_SUFFIX)

enum strand {
    GENOME_STRAND_PLUS=0,
    GENOME_STRAND_MINUS=1
};

enum cdsStat {
    CDSSTAT_NONE,
    CDSSTAT_UNKN,
    CDSSTAT_INCMPL,
    CDSSTAT_CMPL
};

enum seqType {
    TYPE_UNKN,
    TYPE_CODING,
    TYPE_NONCODING,
    TYPE_GENOMIC,
    TYPE_MITOCHONDRIAL,
    TYPE_RNA,
    TYPE_PROTEIN,
};

enum funcType {
    FUNC_UNKNOWN,
    FUNC_INTERGENIC, // intergenic region 
    FUNC_REGULATORY, // regulator, enhancer, sliencer
    FUNC_UPSTREAM, // promoter
    FUNC_5UTR,
    FUNC_NR,
    FUNC_CDS,
    FUNC_SPLITSITE,
    FUNC_SPLITSITE3,
    FUNC_SPLITSITE5,
    FUNC_INTRON,
    FUNC_SPAN,
    FUNC_INITLOSS,
    FUNC_NONCHANGE,
    FUNC_STOPRETAIN,
    FUNC_STOP,
    FUNC_STOPLOSS,
    FUNC_STOPGAIN,
    FUNC_3UTR,
    FUNC_DOWNSTREAM,
};

enum varType {
    VARTYPE_REF,
    VARTYPE_SUBSITITUTION,
    VARTYPE_MISSENSE,
    VARTYPE_NONSENSE,
    VARTYPE_CONVERSION,
    VARTYPE_DELETION,
    VARTYPE_INDEL,
    VARTYPE_INSERTION,
    VARTYPE_INVERSION,
    VARTYPE_DUPLICATION,
    VARTYPE_TRANSLOCATION,
    VARTYPE_TRANSPOSITION,
    VARTYPE_UNKN,
};

enum fsType {
    IS_FRAMESHIFT,
    NOT_FRAMESHIFT,
};

/* UCSC refgene format :
  `bin` smallint(5) unsigned NOT NULL,
  `name` varchar(255) NOT NULL,
  `chrom` varchar(255) NOT NULL,
  `strand` char(1) NOT NULL,
  `txStart` int(10) unsigned NOT NULL,
  `txEnd` int(10) unsigned NOT NULL,
  `cdsStart` int(10) unsigned NOT NULL,
  `cdsEnd` int(10) unsigned NOT NULL,
  `exonCount` int(10) unsigned NOT NULL,
  `exonStarts` longblob NOT NULL,
  `exonEnds` longblob NOT NULL,
  `score` int(11) DEFAULT NULL,
  `name2` varchar(255) NOT NULL,
  `cdsStartStat` enum('none','unk','incmpl','cmpl') NOT NULL,
  `cdsEndStat` enum('none','unk','incmpl','cmpl') NOT NULL,
  `exonFrames` longblob NOT NULL,
 */

struct refgene_entry {
    int empty;
    int flag; // bitwise flag of refgene prase
    int bin;
    char *name1;
    int tid; //    char *chr;
    enum strand strand; // GENOME_STRAND
    int txStart;
    int txEnd;
    int cdsStart;
    int cdsEnd;
    int exonCount;
    int *exonStarts;
    int *exonEnds;
    int score;
    char *name2;
    enum cdsStat cdsStartStat;
    enum cdsStat cdsEndStat;
    int *exonFrames;
    int nfields;
    int *splits;
    kstring_t buffer;
};

struct varCompact {
    enum varType var;
    enum fsType fs;
};
    
struct hgvs_ale {
    struct varCompact type;
    char *cols;
    int fs; // -1 for non-fs
};

struct hgvs1 {
    int empty; // skip if empty == 1
    char *trans;
    char *gene;
    enum strand strand;
    enum seqType type;
    enum funcType func;
    //int n_allele;
    //struct hgvs_ale *als;
    int cpos;	
    int offset; // 0 for cds region
    int exId;
    int cdsId;
};

struct hgvs_record {
    int ntrans, mtrans;
    struct hgvs1 *trans;
};

extern void fill_mempool(bcf_hdr_t *h, int tid, int start, int end, int flag);
extern void release_refgene_mempool();
extern struct hgvs_record *hgvs_generate(bcf_hdr_t *h, bcf1_t *line);

#endif

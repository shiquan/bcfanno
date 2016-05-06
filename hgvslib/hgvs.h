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

#define UTR3_REG  10
#define UTR5_REG  10

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

enum sequenceType {
    TYPE_UNKN,
    TYPE_CODING,
    TYPE_NONCODING,
    TYPE_GENOMIC,
    TYPE_MITOCHONDRIAL,
    TYPE_RNA,
    TYPE_PROTEIN,
};

enum variant_type {
    VARTYPE_REF,
    VARTYPE_SUBSITITUTION,
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

struct hgvs1 {
    enum variant_type type;
    char *name;
    char *alt1;
    int cds_pos;	
    int offset;
    int exon_id;
    int cds_id;
};

struct hgvs_record {
    int tid;
    int start;
    int end;
    int lref, lalt;
    char *ref;
    char *alt;
    int ltrans;
    struct hgvs1 *trans;
};

extern void fill_mempool(int tid, int start, int end, int flag);
extern void release_refgene_mempool();

#endif

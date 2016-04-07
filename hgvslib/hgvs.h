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

/* refgene format praser */
#define REFGENE_PRASE_BIN    1
#define REFGENE_PRASE_NAME1  2
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
    CDSSTAT_INCOMPL,
    CDSSTAT_CMPL
};

enum sequenceType {
    TYPE_UNKN, TYPE_CODING, TYPE_NONCODING,
    TYPE_GENOMIC, TYPE_MITOCHONDRIAL, TYPE_RNA,
    TYPE_PROTEIN
};

enum hgvs_function { FUC};

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
    int extract_flag; // bitwise flag of refgene prase
    int bin;
    char *name;
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
};

struct refgene_mempools {
    int n, m;
    struct refgene_entry *entrys;
    int begin;
    int end;
    int tid;
    int lastbegin;
    int lastend;
};

struct hgvs {
    int tid;
    int start;
    int end;
    char *ref;
    char *alt;
    char *alt1;
    int cds_pos;
    int offset;
    int exon_id;
    int cds_id;
};

#endif

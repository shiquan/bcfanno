// This is a updated version of bedutil.h from bedutils program.
// Copyright shiquan@link.cuhk.edu.hk
//

#ifndef BED_UTILS_HEADER
#define BED_UTILS_HEADER
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "htslib/kstring.h"
#include "htslib/tbx.h"
#include "htslib/kseq.h"
#include "htslib/bgzf.h"

#ifndef KSTRINT_INIT
#define KSTRINT_INIT { 0, 0, 0}
#endif

#define MEMPOOL_LINE 10000 // todo: memory management

// read file handler
KSTREAM_INIT(BGZF*, bgzf_read, 8193);

// flag of bed_file struct
// bits offset rule : right first
//                           bed file is empty or not
//                          /
// uint8_t : | | | | | | | |
//                 | | |  \
//                 | | |   part of file is cached
//                 | | |_  has extra data (more than 3 cols in this bed file), the extra data will get lost in design
//                 | |___  sorted
//                 |_____  merged
#define bed_bit_empty  1
#define bed_bit_cached 2 // when finished, should clear this bit
#define bed_bit_extra  4
#define bed_bit_sorted (1<<3)
#define bed_bit_merged (1<<4)
#define bed_bit_backup (1<<5)

struct bed_line {
    int chrom_id;
    int32_t start;
    int32_t end;    
};

#define BED_LINE_INIT { -1, 0, 0 }

struct bed_chrom {
    int cached; // cached size
    int max; // max allocated memory size
    int i; // for loop, i will be used by bed_read_line()
    uint64_t *a; 
    int id; // name id, usually for chromosomes or contigs
    uint32_t length;    
};

// typedef struct bed_chrom* bed_chrom_point;
// KHASH_MAP_INIT_STR(reg, struct bed_chrom_point)
// typedef kh_reg_t reghash_t;

struct bedaux {
    char *fname;
    uint8_t flag;
    int l_names, m_names;
    char **names;
    // iterator for loop names, used by bed_read_line
    int i;
    // For big file, read first part into memory first and merge and read other parts.
    BGZF *fp; 
    kstream_t *ks;
    uint32_t line;
    // used by bed_fill_bigdata(), if regions are greater than block size, merge cached regions and increase block_size, 
    uint32_t block_size;
    void *hash;
    // original lines|regions
    uint32_t regions_ori;
    // gapped regions in this bed file after operations    
    uint32_t regions;
    uint64_t length_ori;
    // total length of all these chromosomes
    uint64_t length; 
};

extern void set_based_0();
extern void set_based_1();
extern struct bedaux *bedaux_init();

extern void bed_destroy(struct bedaux *bed);

extern struct bed_chrom * get_chrom(struct bedaux *bed, const char *name);
// fork a bedaux structure from bed_chrom
extern struct bedaux *bed_fork(struct bed_chrom *, const char *name, int flag);
extern struct bedaux *bed_dup(struct bedaux *bed);
// read line from chrom structure, return 1 if reach the end, -1 for error, 0 for normal
extern int bed_getline_chrom(struct bed_chrom *chrom, struct bed_line *line);
extern int bed_getline(struct bedaux *bed, struct bed_line *line);
// read a bed file
extern int bed_read_bigfile(struct bedaux *bed, const char *fname);
extern int bed_read(struct bedaux *bed, const char *fname);

// sort
extern int bed_sort(struct bedaux *bed);
// merge
extern int bed_merge(struct bedaux *bed);
extern struct bedaux *bed_merge_several_files(struct bedaux **beds, int n);
// flank | trim
extern void bed_flktrim(struct bedaux *bed, int left, int right);
extern void bed_round(struct bedaux *bed, int length);
// uniq
extern struct bedaux *bed_overlap(struct bedaux *bed);
extern struct bedaux *bed_uniq_several_files(struct bedaux **beds, int n);
extern struct bedaux *bed_uniq_bigfile(struct bedaux *bed, tbx_t *tbx);

// bed_find_rough_bigfile() is an experimental function to find uniq regions and if no uniq region then find
// most nearby regions.
// this function is used to design oligos. required two parameters, gap size and region length
// gap_size for check the nearby regions, if the closest region of target is far than gap size, ignore it.
// region_limit for generate the length of nearby regions, if find a close enough region, the length of this region
// will cap to region_limit.
extern struct bedaux *bed_find_rough_bigfile(struct bedaux *bed, htsFile *fp, tbx_t *tbx, int gap_size, int region_limit);
// diff
extern struct bedaux *bed_diff(struct bedaux *bed1, struct bedaux *bed2);
extern struct bedaux *bed_diff_bigfile(struct bedaux *bed, tbx_t *tbx);

// if bed is raw, just add new line at the end of it
// if bed is sorted, new line will kept in cooridinate,
// if bed is merged, new line will merge into it. 
extern void push_newline(struct bedaux *bed, const char *name, int start, int end);
extern void push_newline1(struct bedaux *bed, struct bed_line *l);

extern int bed_save(struct bedaux *bed, const char *fname);

extern int bed_region_covered(struct bedaux *bed, char *name, int start, int end);
extern int bed_position_covered(struct bedaux *bed, char *chr, int pos, int *_start, int *_end);
#endif

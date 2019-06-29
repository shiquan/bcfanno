#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include "utils.h"
#include "number.h"
#include "bed_utils.h"
#include "htslib/hts.h"
#include "htslib/khash.h"
#include "htslib/ksort.h"

// for very large file, there might be a memory overflow problem to keep all raw data, so here design a read-and-hold
// structure to read some parts of bed file into memory pool, sort and merge cached data first and then load remain 
// data to reduce the memory cost
// #define MEMPOOL_MAX_LINES  10000
static uint32_t mempool_max_lines = 10000;
// if file is greater than FILE_SIZE_LIMIT, just hold the fp handler for bed_read()
// #define FILE_SIZE_LIMIT 100000000
static uint32_t file_size_limit = 10000000; // 10m


void set_memory_max_lines(uint32_t n_lines)
{
    mempool_max_lines = n_lines;
}

void set_file_size_limit(uint32_t limit)
{
    file_size_limit = limit;
}

KSORT_INIT_GENERIC(uint64_t)

// hash structure, chromosome is key, struct bed_chrom is value
typedef struct bed_chrom* bed_chrom_point;
KHASH_MAP_INIT_STR(reg, bed_chrom_point)
typedef kh_reg_t reghash_type;

#ifndef KSTRING_INIT
#define KSTRING_INIT { 0, 0, 0}
#endif
static int based_1 = 0;

void set_based_0(void)
{
    based_1 = 0;
}
void set_based_1(void)
{
    based_1 = 1;
}

#define is_base_1 (based_1 == 1)
#define is_base_0 (based_1 == 0)

struct bedaux *bedaux_init()
{
    struct bedaux *bed = (struct bedaux*)malloc(sizeof(struct bedaux));
    bed->flag = bed_bit_empty;
    bed->l_names = bed->m_names = 0;
    bed->names = 0;
    bed->i = 0;
    bed->hash = kh_init(reg);
    bed->regions_ori = 0;
    bed->regions = 0;
    bed->length_ori = 0;
    bed->length = 0;
    bed->line = 0;
    bed->fname = NULL;
    bed->block_size = mempool_max_lines;
    return bed;
}
struct bed_chrom *bedchrom_init()
{
    struct bed_chrom *chrom = (struct bed_chrom*)malloc(sizeof(struct bed_chrom));
    chrom->cached = chrom->max = 0;
    chrom->a = 0;
    chrom->id = -1;
    chrom->length = 0;
    chrom->i = 0;
    return chrom;
}

void bed_destroy(struct bedaux *file)
{
    if (file == NULL) return;
    khiter_t k;
    int i;
    reghash_type *hash = (reghash_type*)file->hash;
    for (i = 0; i < file->l_names; ++i) {
	char *name = file->names[i];
	k = kh_get(reg, hash, name);
	free(name);
	if (k == kh_end(hash)) {
	    continue;
	} else {
	    struct bed_chrom * chrom = kh_val(hash, k);
	    free(chrom->a);
	    free(chrom);
	    kh_del(reg, hash, k);
	}
    }
    kh_destroy(reg, hash);
    
    free(file);    
}
int get_name_id(struct bedaux *bed, const char *name)
{
    int i;
    for (i = 0; i < bed->l_names; ++i) {
	if ( strcmp(bed->names[i], name) == 0) return i;
    }
    return -1;
}
// for read only
struct bed_chrom *get_chrom(struct bedaux *bed, const char *name)
{
    if (bed->flag & bed_bit_empty) {
	warnings("Try to get chrom structure from empty bed.");
	return NULL;
    }
    khint_t k;
    reghash_type *hash = (reghash_type*)bed->hash;
    k = kh_get(reg, hash, name);
    if (k == kh_end(hash))
        return NULL;
    struct bed_chrom *chrom = kh_val(hash, k);
    return chrom;
}

// return
// 0 for normal
// 1 for empty
// 2 for malformed line
static int parse_string(struct bedaux *bed, kstring_t *string, struct bed_line *line)
{
    int nfields = 0;
    int *splits = ksplit(string, 0, &nfields);
    if ( splits == NULL ) return 1;
    if ( nfields < 2)
	goto malformed_line;
    reghash_type * hash = (reghash_type*)bed->hash;
    khiter_t k;
    k = kh_get(reg, hash, string->s);
    char *name = string->s + splits[0];
    char *temp = string->s + splits[1];
    if ( check_num_likely(temp) )
	line->start = str2int(temp);
    else
        goto malformed_line;

    if ( nfields > 2 && check_num_likely(string->s + splits[2])) {
	line->end = str2int(string->s + splits[2]);
    } else {
	line->end = line->start;
	line->start = line->start < 1 ? 0 : line->start -1;
    }
    k = kh_get(reg, hash, name);
    int id = get_name_id(bed, name);
    if ( k == kh_end(hash) ) {
	if (id == -1) {
	    if ( bed->l_names == bed->m_names ) {
		bed->m_names = bed->m_names == 0 ? 2 : bed->m_names << 1;
		bed->names = (char**)realloc(bed->names, bed->m_names*sizeof(char*));
	    }
	    id = bed->l_names;
	    bed->names[bed->l_names++] = strdup(name);
	}
	int ret;	    
	struct bed_chrom *chrom = bedchrom_init();
	chrom->id = id;
	k = kh_put(reg, hash, bed->names[id], &ret);
	kh_val(hash, k) = chrom;
    }
    line->chrom_id = id;
    free(splits);
    return 0;

  malformed_line:
    free(splits);
    return 2;
}

static int bed_fill(struct bedaux *bed)
{
  //if (bed->flag & bed_bit_empty) return 1;
  //if (bed->flag ^ bed_bit_cached) return 1;

    // reghash_type *hash = (reghash_type*)bed->hash;
    kstring_t string = KSTRING_INIT;
    int dret;
    struct bed_line line = BED_LINE_INIT;
    while ( ks_getuntil(bed->ks, 2, &string, &dret) >= 0) {
	//int start = -1;
	//int end = -1;
	bed->line++;
	if ( string.l == 0 || string.s[0] == '\n' ) {
	    warnings("%s : line %d is empty. skip ..", bed->fname, bed->line);
	    continue;
	}
	if ( string.s[0] == '#' ) continue;
	if ( parse_string(bed, &string, &line) )
	    goto clean_string;

	if ( line.start == line.end && is_base_0 ) {
	    warnings("line %d looks like a 1-based region. Please make sure you use right parameters.", bed->line);
	    --line.start;
	}

	push_newline1(bed, &line);
      clean_string:
	string.l = 0;
    }
    if ( string.m ) free(string.s);
    bgzf_close(bed->fp);
    ks_destroy(bed->ks);
    return 0;
}
static void chrom_sort(struct bed_chrom *chrom)
{
    ks_introsort(uint64_t, chrom->cached, chrom->a);
}
static void chrom_merge(struct bed_chrom *chrom)
{    
    chrom_sort(chrom);
    
    int i;
    uint64_t *b = (uint64_t*)malloc(chrom->cached * sizeof(uint64_t));
    uint32_t start_last = 0;
    uint32_t end_last = 0;
    int l = 0;
    int length = 0;
    for ( i = 0; i < chrom->cached; ++i ) {
	uint32_t start = chrom->a[i]>>32;
	uint32_t end = (uint32_t)chrom->a[i];
	if ( end_last < 1 ) {
	    start_last = start;
	    end_last = end;
	    continue;
	}
	if ( end_last >= start) {
	    if ( end_last < end )
		end_last = end;	    
	} else {
	    b[l++] = (uint64_t) start_last<<32| end_last;
	    length += end_last - start_last;
	    start_last = start;
	    end_last = end;
	}
    }
    // tail region
    if (end_last > 0) {
	length += end_last - start_last;
	b[l++] = (uint64_t) start_last<<32| end_last;
    }
    memset(chrom->a, 0, chrom->cached * sizeof(uint64_t));
    memcpy(chrom->a, b, l * sizeof(uint64_t));
    chrom->cached = l;
    free(b);
    chrom->length = length;
}
void bed_cache_update(struct bedaux *bed)
{
    int i;
    khiter_t k;
    reghash_type *hash = (reghash_type*)bed->hash;
    bed->regions = 0;
    bed->length = 0;
    for (i = 0; i < bed->l_names; ++i) {
	char *name = bed->names[i];
	k = kh_get(reg, hash, name);
	if (k == kh_end(hash)) continue;
	struct bed_chrom *chrom = kh_val(hash, k);
	chrom_merge(chrom);
	bed->regions += chrom->cached;
	bed->length += chrom->length;
    }
    bed->block_size = bed->regions + mempool_max_lines;    
}
int bed_fill_bigdata(struct bedaux *bed)
{
    // reghash_type *hash = (reghash_type*)bed->hash;
    kstring_t string = KSTRING_INIT;
    int dret;
    struct bed_line line = BED_LINE_INIT;
    while ( ks_getuntil(bed->ks, 2, &string, &dret) >= 0) {
	// int start = -1;
	// int end = -1;
	bed->line++;
	if ( string.l == 0 || string.s[0] ) {
	    warnings("%s : line %d is empty. skip ..", bed->fname, bed->line);
	    continue;
	}
	if ( string.s[0] == '#' ) continue;

	if ( parse_string(bed, &string, &line) )
	    goto clean_string;
	
	if ( line.start == line.end && is_base_0 ) {
	    warnings("line %d looks like a 1-based region. Please make sure you use right parameters.", bed->line);
	    line.start--;
	}
	push_newline1(bed, &line);
	if ( bed->regions == bed->block_size) {
	    bed_cache_update(bed);
	}
	
      clean_string:
	string.l = 0;
    }
    if ( string.m ) free(string.s);
    bgzf_close(bed->fp);
    ks_destroy(bed->ks);    
    return 0;
}
int bed_read(struct bedaux *bed, const char *fname)
{
    // if bed size is greater than 100M, or sort already, hold the file handle
    bed->fp = bgzf_open(fname, "r");
    if (bed->fp == 0)
	error("failed to open %s : %s.", fname, strerror(errno));
    bgzf_seek(bed->fp, 0L, SEEK_END);
    uint64_t size = bgzf_tell(bed->fp);

    // go back to file begin    
    bgzf_seek(bed->fp, 0L, SEEK_SET);
    bed->ks = ks_init(bed->fp);
    bed->fname = (char*)fname;
    // remove empty flag
    bed->flag &= ~bed_bit_empty;
   // small file, cached whole file
    if ( size < file_size_limit ) {
	bed_fill(bed);
	bed->flag &= ~bed_bit_cached;
	// file is empty, set empty flag
	if ( bed->length == 0 )
	    bed->flag |= bed_bit_empty;	
	return 1;
    }

    // for huge file, just hold the handler
    bed->flag |= bed_bit_cached;
    return 0;
}
struct bedaux *bed_fork(struct bed_chrom *chrom, const char *name, int flag)
{
    struct bedaux *bed = bedaux_init();
    bed->flag = flag;
    bed->l_names = bed->m_names = 1;
    bed->names = (char**)malloc(sizeof(char*));
    bed->names[0] = strdup(name);
    reghash_type *hash = (reghash_type*)bed->hash;
    khiter_t k;
    int ret;
    k = kh_put(reg, hash, name, &ret);
    kh_val(hash, k) = chrom;
    return bed;
}
struct bed_chrom *bed_chrom_dup(struct bed_chrom *_chm)
{
    struct bed_chrom * chm = bedchrom_init();
    chm->id = _chm->id;
    chm->cached = _chm->cached;
    chm->max = _chm->cached;
    chm->a = (uint64_t*)malloc(sizeof(uint64_t)*chm->max);
    memcpy(chm->a, _chm->a, chm->cached * sizeof(uint64_t));
    chm->length = _chm->length;
    return chm;
}
struct bedaux *bed_dup(struct bedaux *_bed)
{
    if ( _bed->flag & bed_bit_cached )
	error("[bed_dup]bedaux  should be filled. Trying to fork a cached bed struct ..");
    struct bedaux *bed = bedaux_init();
    bed->flag = _bed->flag;
    bed->fname = _bed->fname;
    bed->l_names = _bed->l_names;
    bed->m_names = _bed->m_names;

    int i;
    if ( bed->m_names ) {
	bed->names = (char**)malloc(bed->m_names*sizeof(char*));
	for (i = 0; i < bed->l_names; ++i)
	    bed->names[i] = strdup(_bed->names[i]);	
    }
    bed->hash = kh_init(reg);
    reghash_type *hash = (reghash_type*)bed->hash;
    reghash_type *hash1 = (reghash_type*)_bed->hash;
    khiter_t k;
    khiter_t k1;
    for (i = 0; i < bed->l_names; ++i) {
	k = kh_get(reg, hash, bed->names[i]);
	k1 = kh_get(reg, hash1, bed->names[i]);
	if ( k1 == kh_end(hash1) ) continue;
	if ( k == kh_end(hash)) {
	    int ret;
	    k = kh_put(reg, hash, bed->names[i], &ret);	    
	}
	kh_val(hash, k) = bed_chrom_dup(kh_val(hash1, k1));	    	
    }
    return bed;
}
// 1 for end, 0 for success, -1 for no found
int bed_getline_chrom(struct bed_chrom *chm, struct bed_line *line)
{
    if ( chm == NULL )
        return -1;
    if ( chm->i >= chm->cached )
	return 1;
    
    line->chrom_id = chm->id;
    line->start = chm->a[chm->i] >> 32;
    line->end = (uint32_t)chm->a[chm->i];
    chm->i++;
    return 0;
}
int bed_getline(struct bedaux *bed, struct bed_line *line)
{
    for ( ; bed->i < bed->l_names; bed->i++ ) {
	struct bed_chrom *chm = get_chrom(bed, bed->names[bed->i]);
	if ( bed_getline_chrom(chm, line) == 0)
	  break;
    }
    if (bed->i == bed->l_names)
        return 1;
    
    return 0;    
}
int bed_sort(struct bedaux *bed)
{
    // sorted already
    if ( bed->flag & bed_bit_sorted ) return 1;
    
    int i;
    for (i = 0; i < bed->l_names; ++i) {
	struct bed_chrom *chm = get_chrom(bed, bed->names[i]);
	if (chm == NULL)
	    continue;
	chrom_sort(chm);
    }
    bed->flag |= bed_bit_sorted;
    return 0;
}
int bed_merge(struct bedaux *bed)
{
    if ( bed->flag & bed_bit_merged)
	return 1;
    // bed_sort(bed);
    int i;
    for (i = 0; i < bed->l_names; ++i) {
	struct bed_chrom *chm = get_chrom(bed, bed->names[i]);
	if (chm == NULL)
	    continue;
	chrom_merge(chm);
    }
    bed->flag |= bed_bit_sorted;
    bed->flag |= bed_bit_merged;
    return 0;
}
// todo: 
// require all bed files sorted
struct bedaux *bed_merge_several_bigdata(struct bedaux **beds, int n)
{
    return NULL;
}
struct bedaux *bed_merge_several_files(struct bedaux **beds, int n)
{
    return NULL;
}
void bed_flktrim(struct bedaux *bed, int left, int right)
{
    int i, j;
    uint64_t length = 0;
    for (i = 0; i < bed->l_names; ++i) {
	struct bed_chrom *chm = get_chrom(bed, bed->names[i]);
	if (chm == NULL)
	    continue;
	for (j = 0; j < chm->cached; ++j) {
	    uint32_t start = chm->a[j]>>32;
	    uint32_t end = (uint32_t)chm->a[j];
	    // if region is too short to trim, skip it without a warning
	    if ((int)(end -start) <= (left + right) * -1)
		continue;
	    start -= left;
	    end += right;
	    chm->a[j] = (uint64_t)start << 32|end;
	    chm->length += right + left;
	}
	length += chm->length;
    }
    bed->length = length;
}
void bed_round(struct bedaux *bed, int round_length)
{
    uint64_t length = 0;
    int i, j;
    for (i = 0; i < bed->l_names; ++i) {
	struct bed_chrom *chm = get_chrom(bed, bed->names[i]);
	if (chm == NULL)
	    continue;
	for (j = 0; j < chm->cached; ++j) {
	    uint32_t start = chm->a[j]>>32;
	    uint32_t end = (uint32_t)chm->a[j];
	    if (end - start >= round_length) continue;
	    int offset = round_length - (end -start);
	    start = start - offset/2;
	    end = end + offset/2 + (offset & 1);  // add the extra base to the end
	    chm->a[j] = (uint64_t)start << 32|end;
	    chm->length += offset;
	}
	length += chm->length;
    }
    bed->length = length;
}
struct bedaux *bed_overlap(struct bedaux *bed)
{
    return NULL;
}
struct bedaux *bed_uniq_several_files(struct bedaux **beds, int n)
{
    return NULL;
}
struct bedaux *bed_uniq_bigfile(struct bedaux *bed, tbx_t *tbx)
{
    return NULL;
}
static void copy_line(struct bed_line *dest, struct bed_line *line)
{
    dest->chrom_id = line->chrom_id;
    dest->start = line->start;
    dest->end = line->end;
}
// bed_find_rough_bigfile() is a function to retrieve most nearest or covered regions in the tbx databases for target regions
// for probe design programs, gap_size is usually slightly smaller than the fragement size.
struct bedaux *bed_find_rough_bigfile(struct bedaux *target, htsFile *fp, tbx_t *data, int gap_size, int region_limit)
{
    bed_merge(target);
    struct bed_line line = BED_LINE_INIT;

    kstring_t string = KSTRING_INIT;
    struct bedaux *design = bedaux_init();
    design->flag &= ~bed_bit_empty;
    while ( bed_getline(target, &line) == 0 ) {
	// retrieve target in dataset
	int tid = tbx_name2id(data, target->names[line.chrom_id]);
	if (tid == -1) {
	    warnings("Chromosome %s is not found data.", target->names[line.chrom_id]);
	}
	hts_itr_t *itr = tbx_itr_queryi(data, tid, line.start, line.end);
	int n_regions = 0;
	int left = 0;
	int right = 0;
	struct bed_line dl = BED_LINE_INIT;
	while ( tbx_itr_next(fp, data, itr, &string) >= 0) {	    
	    parse_string(design, &string, &dl);
	    if ( dl.start < line.start ) dl.start = line.start;
	    if ( dl.end > line.end ) dl.end = line.end;
	    if ( left == 0)
		left = dl.start;
	    if ( right < line.end )
		right = line.end;
	    push_newline1(design, &dl);
	    n_regions++;
	    string.l = 0;
	}

	// if there are too much gaps in the edges, or	
	// if no regions in dataset, find nearby regions	

	// find nearest left side regions
        if (n_regions == 0 || left - line.start > gap_size) {
	    uint32_t start = line.start - gap_size > 0 ? line.start - gap_size : 0;
	    uint32_t end = line.start;
	    itr = tbx_itr_queryi(data, tid, start, end);
	    while ( tbx_itr_next(fp, data, itr, &string) >= 0) {
		parse_string(design, &string, &dl);
		if (dl.start < start) dl.start = start;
		push_newline1(design, &dl);
	    }
	    string.l = 0;
	}
	// find nearest right side regions
	if (n_regions == 0 ||  line.end - right > gap_size) {
	    uint32_t start = line.end;
	    uint32_t end = line.end + gap_size;
	    itr = tbx_itr_queryi(data, tid, start, end);
	    while ( tbx_itr_next(fp, data, itr, &string) >= 0) {
		parse_string(design, &string, &dl);
		if (dl.end > end) dl.end = end;
		push_newline1(design, &dl);
	    }
	    string.l = 0;
	}
    }
    bed_merge(design);
    return design;
}
struct bedaux *bed_diff(struct bedaux *bed1, struct bedaux *bed2)
{
    return NULL;
}
struct bedaux *bed_diff_bigfile(struct bedaux *bed, tbx_t *tbx)
{
    return NULL;
}
void push_newline1(struct bedaux *bed, struct bed_line *l)
{    
    if (l->chrom_id == -1 || l->chrom_id > bed->l_names) 
	error("[push_newline1] chrom is not found, id : %d, lname : %d", l->chrom_id, bed->l_names);
    if (l->start > l->end) { int temp = l->end; l->end = l->start; l->start = temp; }
    khiter_t k;
    reghash_type *hash = (reghash_type*)bed->hash;    
    k = kh_get(reg, hash, bed->names[l->chrom_id]);
    struct bed_chrom *chm = kh_val(hash, k);    
    
    if (chm->cached == chm->max) {			   
	chm->max = chm->max == 0 ? 10 : chm->max << 1; 
	chm->a = (uint64_t*)realloc(chm->a, chm->max * sizeof(uint64_t));
    }
    chm->a[chm->cached++] = (uint64_t)l->start << 32 | l->end;
    bed->flag &= ~bed_bit_merged;
    bed->flag &= ~bed_bit_sorted;
    bed->regions_ori++;
    bed->length_ori += l->end - l->start;
    bed->length += l->end - l->start;
    bed->regions++;
}
void push_newline(struct bedaux *bed, const char *name, int start, int end)
{    
}

int bed_save(struct bedaux *bed, const char *fname)
{
#ifdef _DEBUG_MODE
    debug_print("[%s]", __func__);
#endif
    if ( bed == NULL) return 1;
    if ( bed->flag & bed_bit_empty ) return 1;
    if ( bed->flag & bed_bit_cached ) bed_fill(bed);
    
    FILE *fp = fopen(fname, "w");
    khiter_t k;
    int i, j;
    reghash_type * hash = (reghash_type*)bed->hash;
    for (i = 0; i < bed->l_names; ++i) {
	k = kh_get(reg, hash, bed->names[i]);
	if ( k != kh_end(hash) ) {
	    struct bed_chrom * chrom = kh_val(hash, k);
	    if ( chrom == NULL)
		continue;
	    for (j = 0; j < chrom->cached; ++j)
		fprintf(fp, "%s\t%u\t%u\n", bed->names[i], (uint32_t)(chrom->a[j] >> 32), (uint32_t)chrom->a[j]);
	}
    }
    fclose(fp);
    return 0;
}

int bed_region_covered(struct bedaux *bed, char *chr, int start, int end)
{
    struct bed_chrom *c =  get_chrom(bed, chr);
    if ( c == NULL ) return 0; // not found this chromosome

    if ( start < 0 || end < 0 ) error("Trying to find an unreasonable genomic locations.");
    
    uint64_t a = (uint64_t)start<<32|end;
    int min = 0, max = c->cached-1;
    int mid;
    for (;;) {
        mid = (min + max)/2;
        if (c->a[mid] == a || c->a[min] == a || c->a[max] == a) return 1; // found
        if (c->a[mid] > a) max = mid;
        else min = mid;

        if ( max == min ) {
            uint32_t s = c->a[max] >> 32;
            uint32_t e = (uint32_t)c->a[max];
            if ( start >= s && start <= e ) return 1;
            if ( end >= s && end <= e ) return 1;
            break;
        }
        else if (max - min == 1 ) {
            uint32_t e = c->a[max] >> 32;
            uint32_t s = (uint32_t)c->a[min];
            if ( start >= s && end <= e ) break;
            return 1;
        }
    }
    return 0;
}

static int rloc(uint64_t a, int pos, int *_s, int *_e)
{
    uint32_t s = a >>32;
    uint32_t e = (uint32_t)a;
    if ( pos >= s && pos <= e ) {
        *_s = s;
        *_e = e;
        return 1;
    }
    return 0;
}
int bed_position_covered(struct bedaux *bed, char *chr, int pos, int *_start, int *_end)
{
    // int start, end;
    struct bed_chrom *c = get_chrom(bed, chr);
    if ( c == NULL ) return 0;

    if ( pos < 1 ) error("Trying to find an unreasonable genomic location.");

    int mid, min = 0, max = c->cached -1;    
    for ( ;; ) {
        mid = (min + max)/2;
        if ( rloc(c->a[mid], pos, _start, _end) ||
             rloc(c->a[min], pos, _start, _end) ||
             rloc(c->a[max], pos, _start, _end) ) return 1;
        if ( min == max ) 
            return rloc(c->a[min], pos, _start, _end);
        else if ( max - min == 1 ) {
            if ( rloc(c->a[min], pos, _start, _end) ) return 1;
            if ( rloc(c->a[max], pos, _start, _end) ) return 1;
            return 0;
        }
        if ((uint32_t)c->a[mid] < pos ) min = mid;
        else max = mid;
    }
    return 0;
}

#ifdef _MAIN_BED
#include "utils.h"

int main(int argc, char **argv)
{
    if (argc != 2) {
	error("%s in.bed", argv[0]);
    }
    LOG_print("read %s ..", argv[1]);
    struct bedaux *bed = bedaux_init();
    bed_read(bed, argv[1]);
    bed_merge(bed);
    LOG_print("save merged target file target.bed ..");
    bed_save(bed, "target.bed");

    LOG_print("flank 100b of target file to flank.bed ..");
    bed_flktrim(bed, 100, 100);
    bed_save(bed, "flank.bed");
    LOG_print("trim 100b of flanked file to trim.bed ..");
    bed_flktrim(bed, -100, -100);
    bed_save(bed, "trim.bed");
    LOG_print("round trimed file by 10000b to round.bed ..");
    bed_round(bed, 10000);
    bed_save(bed, "round.bed");
    LOG_print("clean ..");
    bed_destroy(bed);
    return 0;
}
#endif

#include "utils.h"
#include "genepred.h"
#include "number.h"
#include "sort_list.h"
#include "zlib.h"
//#include "file.h"
#include "htslib/tbx.h"
#include "htslib/khash.h"
#include "htslib/hfile.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "htslib/bgzf.h"

#define KSTRING_INIT {0, 0, 0 }

KHASH_MAP_INIT_STR(list, char*)
typedef kh_list_t list_hash_t;

#ifdef HTS3

struct stream_lite {
    int begin, end;
    int is_eof:2, bufsize:30;
    uint64_t seek_pos;
    BGZF* f;
    unsigned char *buf;
};

int file_seek(htsFile *fp, long offset, int where)
{
    if ( fp->is_bin ) {
        return bgzf_seek(fp->fp.bgzf, offset, where);
    } else {
        ks_rewind((struct stream_lite*)fp->fp.voidp);
        ((struct stream_lite*)fp->fp.voidp)->seek_pos = offset;
        return bgzf_seek(((struct stream_lite*)fp->fp.voidp)->f, offset, where);
    }
}

#else

int file_seek(htsFile *fp, long offset, int where)
{
    if ( fp->is_bgzf ) {
        return bgzf_seek(fp->fp.bgzf, offset, where);
    } else {        
        return hseek(fp->fp.hfile, offset, where);
    }
}

#endif

#define GENEPRED_CIGAR_UNKNOWN_BASE '*'
#define GENEPRED_CIGAR_MATCH_BASE   'M'
#define GENEPRED_CIGAR_DELETE_BASE  'D'
#define GENEPRED_CIGAR_INSERT_BASE  'I'

// The fields cigar packed counts and cigar type, as follows:
//
//   offset   bits  value
//   0        3     cigar_type
//   4        32    cigar_counts
//
#define GENEPRED_CIGAR_PACKED_FIELD 3
#define GENEPRED_CIGAR_UNKNOWN_TYPE 0
#define GENEPRED_CIGAR_MATCH_TYPE   1
#define GENEPRED_CIGAR_DELETE_TYPE  2
#define GENEPRED_CIGAR_INSERT_TYPE  4
#define GENEPRED_CIGAR_MASK_TYPE    7

struct list *init_list(const char *fn)
{
    if (fn == NULL)
	return NULL;
    int i;
    khiter_t k;
    int ret;
    struct list *list = (struct list*)malloc(sizeof(struct list));
        
    list->reads = hts_readlines(fn, &list->lines);
    if ( list->reads == NULL) {
	fprintf(stderr, "%s : %s\n", fn, strerror(errno));
        free(list);
        return NULL;
    }
    
    if (list->lines == 0)
        goto empty_list;

    list->hash = kh_init(list);
    for ( i=0; i< list->lines; ++i ) {
	char *name = list->reads[i];
        if (name == NULL)
            continue;
        if (*name == '#' || *name == '/')
            continue;

	k = kh_get(list, (kh_list_t*)list->hash, name);
	if (k == kh_end((kh_list_t*)list->hash))
	    k = kh_put(list, (kh_list_t*)list->hash, name, &ret);	
    }    
    return list;
    
  empty_list:
    free(list->reads);
    free(list);
    return NULL;
}

// 1 on failure, 0 on success found.
int check_list(struct list *list, char *name)
{
    if ( list == NULL || list->lines == 0 )
        return 1;
    khint_t k;
    k = kh_get(list, (kh_list_t*)list->hash, name);
    if ( k == kh_end((kh_list_t*)list->hash) )
        return 1;
    return 0;
}
void destroy_list(struct list* list)
{
    if (list == NULL)
        return;
    if ( list->reads ) {
        int i;
        for (i=0; i<list->lines; ++i) 
	    free(list->reads[i]);
	free(list->reads);
    }
   
    if ( list->hash )
	kh_destroy(list, (kh_list_t*)list->hash);
    free(list);
}

struct genepred_spec *genepred_spec_init()
{
    struct genepred_spec *spec = (struct genepred_spec*)malloc(sizeof(struct genepred_spec));
    memset(spec, 0, sizeof(struct genepred_spec));
    return spec;
}

void genepred_spec_destroy(struct genepred_spec *spec)
{
    if ( spec->idx )
        tbx_destroy(spec->idx);
    if ( spec->fai )
        fai_destroy(spec->fai);
    
    if ( spec->genes )
        destroy_list(spec->genes);
    if ( spec->trans )
        destroy_list(spec->trans);

    hts_close(spec->fp);
    free(spec);
    spec = NULL;
}
int genepred_load_data(struct genepred_spec *spec, const char *fname)
{
    if ( spec == NULL )
        error("Invoke genepred_spec_init() before load data.");

    if ( fname == NULL )
        fname = getenv("REFGENE");

    if ( fname == NULL )
        error("No genepred or refgene database specified.");

    spec->data_fname = fname;

    // Load the database file.
    spec->fp = hts_open(fname, "r");
    if (spec->fp == NULL)
        error("%s : %s.", fname, strerror(errno));

    // Load the tabix index.
    spec->idx = tbx_index_load(fname);
    if ( spec->idx == NULL )
        error("Failed to load index of %s", fname);

    return 0;
}
int genepred_load_fasta(struct genepred_spec *spec, const char *fname)
{
    if ( spec == NULL )
        error("Invoke genepred_spec_init() before load fasta file.");

    if ( fname == NULL )
        error("No fasta file specified.");

    spec->fai_fname = fname;

    // Load the faidx of fasta file.
    spec->fai = fai_load(fname);
    if (spec->fai == NULL)
        error("%s : %s.", fname, strerror(errno));

    return 0;
}

int genepred_load_genes(struct genepred_spec *spec, const char *fname)
{
    spec->genes = init_list(fname);
    if ( spec->genes == NULL )
        return 1;
    return 0;
}
int genepred_load_trans(struct genepred_spec *spec, const char *fname)
{
    spec->trans = init_list(fname);
    if ( spec->trans == NULL )
        return 1;
    return 0;
}
struct genepred_line *genepred_line_create()
{
    struct genepred_line *line = (struct genepred_line*)malloc(sizeof(struct genepred_line));
    memset(line, 0, sizeof(struct genepred_line));    
    return line;
}

void clear_genepred_line(struct genepred_line *line)
{
    struct genepred_line **pp = &line;
    if ( line == NULL ) {
        *pp = (struct genepred_line*)malloc(sizeof(struct genepred_line));
    } else {
        if ( line->chrom != NULL)
            free(line->chrom);
        if ( line->name1 != NULL)
            free(line->name1);
        if ( line->name2 != NULL)
            free(line->name2);
        free(line->exons[0]);
        free(line->exons[1]);
        if ( line->loc_parsed ) {
            free(line->loc[0]);
            free(line->loc[1]);
        }
        if ( line->realn )
            free(line->realn);
        if ( line->cigars )
            free(line->cigars);           
    }
    line->loc_parsed = 0;
    memset(*pp, 0, sizeof(struct genepred_line));
}
void genepred_line_destroy(void *line)
{
    if ( line != NULL ) {
        clear_genepred_line((struct genepred_line*)line);
        free(line);
    }
}

static struct genepred_format refgene_formats = {
    .chrom = 2,
    .name1 = 1,
    .name2 = 12,
    .strand = 3,
    .txstart = 4,
    .txend = 5,
    .cdsstart = 6,
    .cdsend = 7,
    .exon_count = 8,
    .exonstarts = 9,
    .exonends = 10,
    .name_version = -1,
    .realn = -1,
};
    
static struct genepred_format genepred_formats = {    
    .name1 = 0,
    .chrom = 1,
    .strand = 2,
    .txstart = 3,
    .txend = 4,
    .cdsstart = 5,
    .cdsend = 6,
    .exon_count = 7,
    .exonstarts = 8,
    .exonends = 9,
    .name2 = 10,
    .name_version = 11,
    .realn = 12,
};

static struct genepred_format refflat_formats = {    
    .name1 = 0,
    .name2 = 1,
    .chrom = 2,
    .strand = 3,
    .txstart = 4,
    .txend = 5,
    .cdsstart = 6,
    .cdsend = 7,
    .exon_count = 8,
    .exonstarts = 9,
    .exonends = 10,
    .name_version = -1,
    .realn = -1,
};

// The type defined the format of database. Default is genepred format. Access genepred files of
// species from UCSC table browsers.
static struct genepred_format *type = &genepred_formats;

void set_format_refgene()
{
    type = &refgene_formats;
}
void set_format_genepred()
{
    type = &genepred_formats;
}
void set_format_refflat()
{
    type = &refflat_formats;
}
char check_strand(char *strand)
{
    if (memcmp(strand, "+", 1) == 0 )
        return '+';
    else if (memcmp(strand, "-", 1) == 0)
        return '-';
    error("Unknow strand type, %s", strand);
}
static int parse_line_core(kstring_t *string, struct genepred_line *line)
{
    clear_genepred_line(line);
    int *splits;
    int nfields;
    splits = ksplit(string, 0, &nfields);

#define BRANCH(_node, _col, _func) do { \
        if ( _col >= nfields ) { \
            _node = 0; \
        } else { \
            _node = _func(string->s + splits[_col]);\
        }\
    } while(0)

    // Parse chromosome name.
    BRANCH(line->chrom, type->chrom, strdup);
    if ( line->chrom == NULL )
        return 1;

    // Parse transcripts name or ensembl gene ID.
    BRANCH(line->name1, type->name1, strdup);

    // Parse strand, char '+' or '-'.
    BRANCH(line->strand, type->strand, check_strand);

    // Parse transcripts start.
    BRANCH(line->txstart, type->txstart, str2int);
    // Convert 0 based txstart to 1 based.
    line->txstart++;

    // Parse transcripts end.
    BRANCH(line->txend, type->txend, str2int);
    if ( type->txend == 0 )
        return 1;

    // Parse cds start, for mRNA cds start should be the downstream of txstart; for noncoding RNA
    // cds start should be equal to txend.
    BRANCH(line->cdsstart, type->cdsstart, str2int);
    
    // Parse cds end. for noncoding RNA cdsend == cdsstart == txend.
    BRANCH(line->cdsend, type->cdsend, str2int);

    // Convert 0 based cdsstart to 1 based. Only for coding transcript.
    if ( line->cdsstart < line->cdsend )
        line->cdsstart++;

    // Parse exon number.
    BRANCH(line->exon_count, type->exon_count, str2int);

    // Parse gene name, for some genepred file, no name2 specified.
    BRANCH(line->name2, type->name2, strdup);

    if (type->name_version != -1 && nfields > type->name_version ) {
        BRANCH(line->name_version, type->name_version, str2int);
    }
    else {
        line->name_version = 0;
    }
    
    if (type->realn != -1 && nfields > type->realn) {
        BRANCH(line->realn, type->realn, strdup);
    }
    else {
        line->realn = NULL;
    }
    
#undef BRANCH

    // The genepred file from UCSC define exon regions like 1,2,3,4,5. Here init exon_pair[] and
    // exon_offset_pair[].
    if ( type->exonstarts >= nfields )
        return 1;
    if ( type->exonends >= nfields )
        return 1;

    char *ss, *se, *ss1, *se1;
    
    ss = string->s + splits[type->exonstarts];
    se = string->s + splits[type->exonends];
    
    int i;
    // Alloc memory for exons[], exon_offset_pair[], and loc[].
    for ( i = 0; i < 2; i++ ) {
        line->exons[i] = (int*)calloc(line->exon_count, sizeof(int));
    }
    // Release splits.
    free(splits);

    for ( i = 0; i < line->exon_count; ++i ) {
        // Parse start.
        ss1 = ss;
        while ( ss1 && *ss1 != ',' )
            ss1++;
        ss1[0] = '\0';
        // Convert 0 based start to 1 based ?
        line->exons[BLOCK_START][i] = str2int(ss) +1;
        ss = ++ss1;
        
        // Parse end.
        se1 = se;
        while ( se1 && *se1 != ',')
            se1++;
        se1[0] = '\0';
        line->exons[BLOCK_END][i] = str2int(se);
        se = ++se1;   
    }
    
    return 0;
}
int parse_line(kstring_t *string, struct genepred_line *line)
{
    char *temp = strdup(string->s);
    int ret = 0;
    if ( temp == NULL || temp[0] == '#' )
	ret = 1;
    else if ( parse_line_core(string, line) == 1 ) {
        error_print("Format error. Failed to parse line, %s", temp);
	ret = 1;
    }
    // Release memory.
    free(temp);
    return ret;
}
int parse_line_locs(struct genepred_line *line)
{
    if (line->loc_parsed == 1 )
        error("Double parsed.");
    
    line->loc_parsed = 1;
    
    // parse cigar string
    if (line->realn != NULL ) {
        int i, j, l, m = 0;
        l = strlen(line->realn);

        for ( i = 0; i < l;) {
            j = i;
            for (; isdigit(line->realn[i]) && i < l; i++);
            if ( i == l ) {
                error_print("Failed to parse CIGAR column. %s %d %d %s", line->chrom, line->txstart, line->txend, line->realn);
                return 1;
            }
            if ( line->n_cigar == m ) {
                m = m == 0 ? 2: m + 2;
                line->cigars = (int*)realloc(line->cigars, m * sizeof(int));
            }

#define BRANCH(x, type) do {                                            \
                int c = str2int_l(line->realn+j, i-j);                  \
                x = (c<<GENEPRED_CIGAR_PACKED_FIELD) | (type & GENEPRED_CIGAR_MASK_TYPE); \
            } while(0)
              
            switch ( line->realn[i] ) {
                case GENEPRED_CIGAR_UNKNOWN_BASE: {
                    if ( m > 0 ) 
                        free(line->cigars);
                    line->n_cigar = 0;
                    line->cigars = NULL;              
                    m = 0;
                    break;
                }
                    
                case GENEPRED_CIGAR_MATCH_BASE:
                    // Symbol M denote match and mismatch, no gaps between transcript and genome.
                    BRANCH(line->cigars[line->n_cigar], GENEPRED_CIGAR_MATCH_TYPE);
                    break;
                    
                case GENEPRED_CIGAR_INSERT_BASE:
                    // Insertion in the transcript sequence.
                    BRANCH(line->cigars[line->n_cigar], GENEPRED_CIGAR_INSERT_TYPE);
                    break;
                      
                case GENEPRED_CIGAR_DELETE_BASE:
                    BRANCH(line->cigars[line->n_cigar], GENEPRED_CIGAR_DELETE_TYPE);
                    break;
                    
                default:
                    error_print("Unknown CIGAR type %c", line->realn[i]);
                    return 1;
            }
#undef BRANCH
            // not allocated memory for cigars, get out of realign...
            if ( m == 0 )
                break;
            line->n_cigar++;
            i++;
        }                
    }

    // Calculate the length of function regions, utr5_length is the length of UTR5, and utr3_length
    // is the length of UTR3, without consider of the strand of transcript.
    // int read_length = 0;
    int utr5_length = 0;
    int utr3_length = 0;
    // int cds_length = 0;
    int loc = 0;
    int i;
    int exon_start, exon_end, exon_length;
    // Check the strand.
    int is_strand = line->strand == '+';
    // Check the transcript type, for noncoding RNA cdsstart == cdsend.
    int is_coding = line->cdsstart < line->cdsend ? 0 : 1;
    
    for ( i = 0; i < 2; i++ ) {
        line->loc[i] = (int*)calloc(line->exon_count, sizeof(int));
    }
            
    // First loop. Purpose of this loop is trying to calculate the forward and backward length.
    // Meanwhile, the related location of the transcripts block edges will also be calculated.
    for ( i = 0; i < line->exon_count; ++i ) {

        exon_start = read_start(line->exons, i);
        exon_end = read_end(line->exons, i);
        // Because we just convert 0 based start to 1 based, so here we need plus 1 to calculate exon_length.
        exon_length = exon_end - exon_start +1;

        // Set related location of the transcript block edges. Assume all the transcripts are plus strand in
        // the first plus, then translocate the blocks if transcript is minus.
        line->loc[BLOCK_START][i] = ++loc;
        loc = loc + exon_length - 1;
        line->loc[BLOCK_END][i] = loc;

        // Add the exon length to dna reference length. The reference length should be the sum of all exons.
        line->reference_length += exon_length;

        // If noncoding transcript skip next steps.
        if ( is_coding == 0 )
            continue;

        // Count forward length.
        if ( line->exons[BLOCK_END][i] <= line->cdsstart ) {
            utr5_length += exon_length;
        } else if ( line->cdsstart > line->exons[BLOCK_START][i] ) {
            // First cds, exon consist of UTR and cds.
            utr5_length += line->cdsstart - line->exons[BLOCK_START][i];
        }
        if ( line->cdsend <= line->exons[BLOCK_START][i]) {
            utr3_length += exon_length;
        } else if ( line->cdsend < line->exons[BLOCK_END][i]) {
            utr3_length += line->exons[BLOCK_END][i] - line->cdsend;
        }
    }

    // read_length is the coding reference length for coding transcript or length of noncoding transcript
    // if ( is_coding == 0 ) {
    // read_length = line->reference_length;
    // }  else {
    //    read_length = line->reference_length - utr5_length - utr3_length;
    // Init forward and backward length. For minus strand, reverse the backward and forward length.
    if ( is_strand ) {
        line->utr5_length = utr5_length;
        // line->utr3_length = utr3_length;
        line->cds_length = line->reference_length - utr3_length;
    } else {
        line->utr5_length = utr3_length;
        // line->utr3_length = utr5_length;
        line->cds_length = line->reference_length - utr5_length;
        for ( i = 0; i < line->exon_count; ++i ) {
            line->loc[BLOCK_START][i] = line->reference_length - line->loc[BLOCK_START][i] + 1;
            line->loc[BLOCK_END][i] = line->reference_length - line->loc[BLOCK_END][i] + 1;
        }
    }
    // }

    // For minus strand, reverse the locations.
    // if ( is_strand == 0 ) {
    // for ( i = 0; i < line->exon_count; ++i ) {
    // line->loc[BLOCK_START][i] = line->reference_length - line->loc[BLOCK_START][i] + 1;
    // line->loc[BLOCK_END][i] = line->reference_length - line->loc[BLOCK_END][i] + 1;
    // }
    // }
    
    // realign genome locations
    if ( line->n_cigar > 0 ) {
        int match = 0;
        int del = 0;
        int ins = 0;
        int offset = 0;
        
        int i = 0, j = 0;
        int lock_utr5_realign = is_coding ? 0 : 1;
        int lock_cds_realign = is_coding ? 0 : 1;

        for ( ;; ) {
            // Offset need to be added to utr5 length only if first match block smaller than it.
            if ( lock_utr5_realign == 0 && line->utr5_length <= match) {
                line->utr5_length += offset;
                lock_utr5_realign = 1;
            }
            
            // Adjust the cds length.
            if ( lock_cds_realign == 0 && line->cds_length >= match) {
                line->cds_length += offset;
                lock_cds_realign = 1;
            }

            if ( i == line->exon_count*2 || j == line->n_cigar )
                break;
            
            if ( line->cigars[j] & GENEPRED_CIGAR_MATCH_TYPE ) {
                match += line->cigars[j] >> GENEPRED_CIGAR_PACKED_FIELD;
            }
            // deletions should be consider as match when count locs
            else if ( line->cigars[j] & GENEPRED_CIGAR_DELETE_TYPE ) {
                del = line->cigars[j] >> GENEPRED_CIGAR_PACKED_FIELD;
                match += del;
                offset -= del;
            }
            // If insertion is tail-As, skip count offset.
            else if ( line->cigars[j] & GENEPRED_CIGAR_INSERT_TYPE ) {
                ins = line->cigars[j] >> GENEPRED_CIGAR_PACKED_FIELD;
                if ( j != line->n_cigar -1) 
                    offset += ins;
            }
            for (; i < line->exon_count*2;) {
                int *loc = is_strand ? &line->loc[i%2][i/2] : &line->loc[i%2 ? 0 : 1][line->exon_count-i /2-1];
                if (*loc > match)
                    break;
                
                *loc += offset;
                i++;
                // debug_print("%d, %d, ", i,*loc, match);
            }
            j++;
        }
    }

    return 0;
}

int read_line(struct genepred_spec *spec, kstring_t *string)
{
    // int dret, ret;
    for ( ;; ) {
        if ( hts_getline(spec->fp, KS_SEP_LINE, string) < 0 )
            return 1;

        if ( string->l == 0 )
            continue;
        if ( string->s[0] == '#' || string->s[0] == '/')
            continue;
        
        break;
    }
    return 0;
}
// return 0 on success checked.
int check_gene_trans(struct genepred_spec *spec,  struct genepred_line *line)
{
    if ( spec->genes ) {
        if ( check_list(spec->genes, line->name2) == 0 ) {
            return 0;
        } else {        
            if ( spec->trans ) {
                // Only check transcript name, no version check.
                if ( check_list(spec->trans, line->name1) == 0 ) {
                    return 0;
                }
                return 1;
            }
        }
    } else {
        if ( spec->trans ) {
            // Only check transcript name, no version check.
            if ( check_list(spec->trans, line->name1) == 0 ) {
                return 0;
            }
            return 1;
        } else {
            // if no gene list or transcript list specified, skip check.
            return 0;
        }
    }
    return 1;
}

int genepred_read_line(struct genepred_spec *spec, struct genepred_line *line)
{
    kstring_t string = KSTRING_INIT;
    
    for ( ;; ) {
        if ( read_line(spec, &string) == 1)
            break;
        if ( parse_line(&string, line))
	    continue;
        if ( check_gene_trans(spec, line) == 1 )
            continue;
        if ( string.m ) 
            free(string.s), string.l = 0, string.m = 0;
        return 0;
    }
    if ( string.m ) 
        free(string.s), string.l = 0, string.m = 0;

    return 1;
}
struct genepred_line *genepred_line_copy(struct genepred_line *line)
{
    if ( line == NULL )
        error("Try to copy a NULL line. Are you serious ?");

    int i;
    struct genepred_line *nl = genepred_line_create();
    //memcpy(nl, line, sizeof(struct genepred_line));
    nl->next = line->next;
    nl->chrom = strdup(line->chrom);
    nl->txstart = line->txstart;
    nl->txend = line->txend;
    nl->strand = line->strand;
    nl->name1 = strdup(line->name1);
    nl->name2 = strdup(line->name2);
    nl->cdsstart = line->cdsstart;
    nl->cdsend = line->cdsend;
    nl->utr5_length = line->utr5_length;
    //nl->utr3_length = line->utr3_length;
    nl->cds_length = line->cds_length;
    nl->reference_length = line->reference_length;
    nl->exon_count = line->exon_count;
    nl->loc_parsed = line->loc_parsed;
    for ( i = 0; i < 2; ++i ) {
        nl->exons[i] = (int*)calloc(line->exon_count, sizeof(int));
        memcpy(nl->exons[i], line->exons[i], sizeof(int) *line->exon_count);
    }
    if ( nl->loc_parsed ) {
        for ( i = 0; i < 2; ++i ) {
            nl->loc[i] = (int*)calloc(line->exon_count, sizeof(int));
            memcpy(nl->loc[i], line->loc[i], sizeof(int) *line->exon_count);
        }
    }
    if ( line->realn )
        nl->realn = strdup(line->realn);

    if ( line->n_cigar ) {
        nl->n_cigar = line->n_cigar;
        nl->cigars = (int*)malloc(line->n_cigar*sizeof(int));
        memcpy(nl->cigars, line->cigars, line->n_cigar*sizeof(int));
    }
    
    return nl;
}
void genepred2line(struct genepred_line *line, kstring_t *str)
{
    int i;
    str->l = 0;
    ksprintf(str, "%s\t%s\t%c\t%d\t%d\t%d\t%d\t%d\t",
            line->name1, line->chrom, line->strand, line->txstart-1, line->txend, line->cdsstart-1, line->cdsend, line->exon_count
        );
    for ( i = 0; i < line->exon_count; ++i )
        ksprintf(str, "%d,", line->exons[BLOCK_START][i]-1);

    kputc('\t', str);

    for ( i = 0; i < line->exon_count; ++i )
        ksprintf(str, "%d,", line->exons[BLOCK_END][i]);
    kputc('\t', str);
    kputs(line->name2, str);

}

char *generate_dbref_header()
{
    kstring_t string = KSTRING_INIT;
    kputs("#Chrom\tStart\tEnd\tStrand\tGene\tTranscript\tExon\tStart(p.)\tEnd(p.)\tStart(c.)\tEnd(c.)", &string);
    return string.s;
}
void generate_dbref_database(struct genepred_line *line)
{
    int i;
    kstring_t temp[2] = { KSTRING_INIT, KSTRING_INIT };
    // int types[2];
    for ( i = 0; i < line->exon_count; ++i ) {
        temp[0].l = temp[1].l = 0;        
        int exon_id;
        exon_id = line->strand == '+' ? i + 1 : line->exon_count -i;
        if ( line->loc[BLOCK_START][i] < line->utr5_length ) {
            kputc('-', &temp[0]);
            kputw(line->utr5_length - line->loc[BLOCK_START][i] + 1, &temp[0]);
        } else if ( line->loc[BLOCK_START][i] > line->cds_length ) {
            kputc('*', &temp[0]);
            kputw(line->loc[BLOCK_START][i] - line->cds_length, &temp[0]);
        } else {
            kputw(line->loc[BLOCK_START][i] - line->utr5_length, &temp[0]);
        }

        if ( line->loc[BLOCK_END][i] < line->utr5_length ) {
            kputc('-', &temp[1]);
            kputw(line->utr5_length - line->loc[BLOCK_END][i] + 1, &temp[1]);
        } else if ( line->loc[BLOCK_END][i] > line->cds_length ) {
            kputc('*', &temp[1]);
            kputw(line->loc[BLOCK_END][i] - line->cds_length, &temp[1]);
        } else {
            kputw(line->loc[BLOCK_END][i] - line->utr5_length, &temp[1]);
        }            
	// format: CHROM,START,END,STRAND, GENE, TRANSCRIPT, EXON, START_LOC, END_LOC
        fprintf(stdout, "%s\t%d\t%d\t%c\t%s\t%s\tEX%d\t%d\t%d\t%s\t%s\n",
                line->chrom,  // chromosome
                read_start(line->exons, i)-1, // start, 0 based
                read_end(line->exons,i),
                line->strand,
                line->name2,
                line->name1,
                exon_id,
                read_start(line->loc, i),
                read_end(line->loc, i),
                temp[0].s,
                temp[1].s);
        
    }
    free(temp[0].s);
    free(temp[1].s);
}
struct genepred_line *genepred_retrieve_gene(struct genepred_spec *spec, const char *name)
{
    kstring_t string = KSTRING_INIT;
    struct genepred_line *head = NULL;
    struct genepred_line *temp = NULL;
    struct genepred_line node;
    memset(&node, 0, sizeof(node));
    // Rewind file.
    if ( file_seek(spec->fp, 0, SEEK_SET) < 0 )
        return NULL;

    // For gene could transcript to several transcripts, so need to read all database through to retrieve all transcripts.
    for ( ;; ) {
        if ( read_line(spec, &string) )
            break;
        if ( parse_line(&string, &node) )
	    continue;
        if ( strcasecmp(node.name2, name) == 0 ) {
            struct genepred_line *temp1 = genepred_line_copy(&node);
            if ( head == NULL )
                head = temp1;
            if ( temp ) {
                temp->next = temp1;
            }
            temp = temp1;
        }
    }
    clear_genepred_line(&node);
    if ( string.m )
        free(string.s);
    if ( head == NULL )
        return NULL;
    return head;
}
struct genepred_line *genepred_retrieve_trans(struct genepred_spec *spec, const char *name)
{
    kstring_t string = KSTRING_INIT;
    struct genepred_line *head = NULL;
    struct genepred_line *temp = NULL;
    struct genepred_line node;
    memset(&node, 0, sizeof(node));
    char *ss = (char*)name;
    int check_version = 0;
    int length;
    length = strlen(name);
    for ( ; ss != NULL && *ss; ss++ ) {
        if ( *ss == '.' ) {
            check_version = 1;
            length = ss-name;
            break;
        }
    }

    if ( file_seek(spec->fp, 0, SEEK_SET) < 0 )
        return NULL;

    // For transcript, may align to different contigs or alternative locuses, so also need to read all database throght.    
    for ( ;; ) {
        if ( read_line(spec, &string) )
            break;
        if ( parse_line(&string, &node) )
	    continue;

        if ( check_version == 0 ) {
            char *ss;
            //int i;
            ss = node.name1;
            if ( strncmp(name, node.name1, length) != 0 )
                continue;
            /* for ( i = 0; ss; ++ss, ++i) { */
            /*     if ( *ss == '.') */
            /*         break; */
            /*     if (name[i] != *ss) */
            /*         break; */
            /* } */
            /* if ( *ss != '.' || name[i] != '\0' ) */
            /*     continue; */
            
            // if (strncasecmp(node.name1, name, i ) == 0 ) {
            struct genepred_line *temp1 = genepred_line_copy(&node);
            if ( head == NULL )
                head = temp1;
            if ( temp ) {
                temp->next = temp1;
            }
            temp = temp1;
        } else if ( strcasecmp(node.name1, name) == 0 ) {            
            struct genepred_line *temp1 = genepred_line_copy(&node);
            if ( head == NULL )
                head = temp1;
            if ( temp ) {
                temp->next = temp1;
            }
            temp = temp1;
        }
    }
    clear_genepred_line(&node);
    if ( string.m )
        free(string.s);
    return head;
}

struct genepred_line *genepred_retrieve_region(struct genepred_spec *spec, char *name, int start, int end)
{
    int id;
    id = tbx_name2id(spec->idx, name);
    if ( id == -1 )
        return NULL;

    hts_itr_t *itr = tbx_itr_queryi(spec->idx, id, start, end);
    kstring_t string = KSTRING_INIT;
    struct genepred_line *head = NULL;
    struct genepred_line *temp = NULL;
    while ( tbx_itr_next(spec->fp, spec->idx, itr, &string) >= 0 ) {
        struct genepred_line *line = genepred_line_create();
        if ( parse_line(&string, line) )
	    continue;

        if ( check_gene_trans(spec, line) == 1 )
            continue;
        
        if ( head == NULL ) {
            head = line;
        }
        if ( temp )
            temp->next = line;
        temp = line;
        string.l = 0;
    }    
    free(string.s);
    tbx_itr_destroy(itr);
    return head;
}

#ifdef GENEPRED_TEST_MAIN

struct args {
    const char *format;
    int noheader;
    const char *fast;
    struct genepred_spec *spec;
} args = {
    .format = NULL,
    .fast = NULL,
    .noheader = 0,
    .spec = NULL,
};
void destroy_args()
{
    genepred_spec_destroy(args.spec);
}
int usage()
{
    fprintf(stderr, "retrievebed\n");
    fprintf(stderr, "    -nm transcripts.txt\n");
    fprintf(stderr, "    -gene genes.txt\n");
    fprintf(stderr, "    -fast < one gene or transcript name>\n");
    fprintf(stderr, "    -format [ genepred | refgene | refflat ]\n");
    fprintf(stderr, "    -noheader\n");
    fprintf(stderr, "   [genepred.tsv.gz]\n");
    return 1;
}
int parse_args(int argc, char **argv)
{
    if ( argc == 1)
        return usage();
    --argc, ++argv;
    int i;
    const char *data_fname = 0;
    const char *transcripts = 0;
    const char *genes = 0;
    for ( i = 0; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if ( strcmp(a, "-nm") == 0 && transcripts == 0 )
            var = &transcripts;
        else if ( strcmp(a, "-gene") == 0 && genes == 0 )
            var = &genes;
        else if ( strcmp(a, "-format") == 0 && args.format == 0)
            var = &args.format;
        else if ( strcmp(a, "-fast") == 0 && args.fast == 0)
            var = &args.fast;
        if (var != 0 ) {
            if ( i == argc ) {
                fprintf(stderr,"Missing an argument after %s", a);
                return 1;
            }
            *var = argv[i++];
            continue;
        }
        if ( strcmp(a, "-noheader") == 0 ) {
            args.noheader = 1;
            continue;
        }
        if ( data_fname == 0) {
            data_fname = a;
            continue;
        }
        fprintf(stderr,"Unknown argument : %s", a);
        return 1;
    }
    
    if ( genes == 0 && transcripts == 0 && args.fast == 0)
        return usage();

    args.spec = genepred_spec_init();
    if ( args.format == 0 )
        args.format = "genepred";

    if ( strcmp(args.format, "genepred") == 0 ) {
        set_format_genepred();
    } else if ( strcmp(args.format, "refgene") == 0 ) {
        set_format_refgene();
    } else if ( strcmp(args.format, "refflat") == 0 ) {
        set_format_refflat();
    } else {
        fprintf(stderr, "Unknown format, %s", args.format);
        return 1;
    }
    // if -fast specified, ignore gene or transcript list
    if ( args.fast == 0) {
        genepred_load_genes(args.spec, genes);
        genepred_load_trans(args.spec, transcripts);
    }

    genepred_load_data(args.spec, data_fname);
    return 0;
}

int retrieve_bed()
{
    if ( args.noheader == 0 ) {
        char *header = generate_dbref_header();
        puts(header);
        free(header);
    }
    if ( args.fast ) {
        struct genepred_line *node = NULL;
        node = genepred_retrieve_gene(args.spec, args.fast);
        if ( node == NULL )
            node = genepred_retrieve_trans(args.spec, args.fast);
        struct genepred_line *temp = node;
        for ( ; temp; temp = temp->next )  {
            if ( parse_line_locs(temp) )
                return 1;
            generate_dbref_database(temp);
        }
        list_lite_del(&node, genepred_line_destroy);
    } else {
        struct genepred_line node;
        memset(&node, 0, sizeof(struct genepred_line));
        for ( ;; ) {
            if ( genepred_read_line(args.spec, &node) == 0 ) {
                parse_line_locs(&node);
                generate_dbref_database(&node);
            } else {
                break;
            }
        }
        clear_genepred_line(&node);
    }
    return 0;
}
int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    if ( retrieve_bed() )
        return 1;
    destroy_args();
    
    return 0;
}

#endif

/*
    Copyright (C) 2016,2017  BGI Research

    Author: Shi Quan <shiquan@genomics.cn>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE. 
*/

// Annotate bed format to bedDetail format. Description of bedDetail refer to
// http://rohsdb.cmb.usc.edu/GBshape/FAQ/FAQformat.html#format1.7

#include "utils.h"
#include "genepred.h"
#include "number.h"
#include "sort_list.h"
#include "htslib/hts.h"
#include "htslib/tbx.h"
#include "htslib/kstring.h"

const char *annobed_version = "1.0";
// History:
// v1.0   - init version, simple annotate simple 3 columns bed file with genepred and id databases

struct args {

    // input file, force to be 3 column BED format
    const char *input_fname;

    //  genepredPlus format database to annotate Gene and cds position
    const char *data_fname;

    // ID database, four columns seperated by tab: chrom,start(0 based),end,id(RS id or genome loci id)
    const char *name_fname;

    // track name for annotated bedDetail format
    const char *track_name;

    // description string for annotated bedDetail format
    const char *description;

    // db mark for annotated bedDetail format, should be consistant with UCSC genome version
    const char *db_version;

    // if set , output standard bedDetail format with track name
    int         UCSC_bedDetails_flag;
    
    htsFile    *fp_data;
    htsFile    *fp_name;
    tbx_t      *idx_data;
    tbx_t      *idx_name;
    
} args = {
    .input_fname = NULL,    
    .data_fname  = NULL,
    .name_fname  = NULL,
    .track_name  = "track name",
    .description = "description",
    .db_version  = "unknown",
    .fp_data     = NULL,
    .fp_name     = NULL,
    .idx_data    = NULL,
    .idx_name    = NULL,
};

int usage()
{
    fprintf(stderr, "annobed [options] [input.bed]\n"
            "-data  <FILE>   genepredPlus database\n"
            "-name  <FILE>   database for annotate name column\n"
            "-ucsc           output UCSC bedDatail format flag\n"
        );    
    return 1;
}
struct bedDetail {
    char *chrom;
    int   start;
    int   end;
    char *name;
    char *id;
    char *description;
};

struct IDdata {
    char *chrom;
    int   start;
    int   end;
    char *extra;
};

void clean_IDdata(struct IDdata *IDdata)
{
    // clean buffer
    if ( IDdata->chrom )
        free(IDdata->chrom);

    if ( IDdata->extra )
        free(IDdata->extra);
    memset(IDdata, 0, sizeof(struct IDdata));    
}
static int parse_IDdata(char *str, int l, struct IDdata *IDdata)
{
    int i;
    int col = 0;
    int start = 0;

    for ( i = 0; i < l; ++i) {

#define BRANCH(_key, _func, _len) do {                                  \
            if ( str[i] == '\t' || str[i] == '\n' || str[i] == '\0' || i == l -1) {\
                _key = _func(str+start, _len);\
                col++;\
                start = i + 1;\
            }\
        } while(0)\
            
        if ( col == 0 ) {
            BRANCH(IDdata->chrom, strndup, i);
        }
        else if ( col == 1 ) {
            BRANCH(IDdata->start, str2int_l, i - start);
        }
        else if ( col == 2 ) {
            BRANCH(IDdata->end, str2int_l, i - start);
        }
        else if ( col == 3 ) {
            BRANCH(IDdata->extra, strndup, i - start);
        }
        else {
            break;
        }
#undef BRANCH
    }
    
    return 0;
}

struct transLite {
    // position in gene coordinate
    int pos;
    // location for function region (UTR,CDS)
    int loc;
    int offset;    
    char strand;
    int exon_id;
};

static int find_the_block(struct genepred_line *line, int *blk_start, int *blk_end, int pos)
{
    *blk_start = 0;
    *blk_end = 0;

    int i;
    for ( i = 0; i < line->exon_count; ++i ) {
        int start = line->exons[BLOCK_START][i];
        int end = line->exons[BLOCK_END][i];
        if ( pos < start) {
            *blk_end = i;
            break;
        }
        else {
            *blk_start = i;
            if ( pos <= end ) {
                *blk_end = i;
                break;
            }
        }
    }
    // Always return 0 ? if out of range return 1??
    return 0;
}
static int genepred_find_location(struct genepred_line *line, int *pos, int *offset, int start, int *exon_id)
{
    parse_line_locs(line);
    
    int block1 = 0;
    int block2 = 0;
    if ( find_the_block(line, &block1, &block2, start ) )
        return 1;

    if ( block1 == block2 ) {

        if ( line->strand == '+' ) {
            *pos = line->loc[BLOCK_START][block1] + start - line->exons[BLOCK_START][block1];
        }
        else {
            *pos = line->loc[BLOCK_END][block1] + ( line->exons[BLOCK_END][block1] - start);
        }
        *offset = 0;
    }
    else {
        int upstream =  start - line->exons[BLOCK_END][block1];
        int downstream = line->exons[BLOCK_START][block2]-start;
        if ( upstream > downstream ) {
            *pos = line->loc[BLOCK_START][block2];
            *offset = line->strand == '+' ? -downstream : downstream;
        }
        else {
            *pos = line->loc[BLOCK_END][block1];
            *offset = line->strand == '+' ? upstream : -upstream;
        }
    }
    *exon_id = block1;

    // adjust position if located in exon and there are realignments in the transcript sequence
    if ( *offset != 0 )
        return 0;

    if ( line->n_cigar == 0 )
        return 0;

    if ( line->n_cigar == 1 && (line->cigars[0] & GENEPRED_CIGAR_MATCH_TYPE) )
        return 0;
    
    int adjust = 0;
    int i;
    int match = 0;
    int ins, del;
    for ( i = 0; i < line->n_cigar; ++i ) {            
        if ( line->cigars[i] & GENEPRED_CIGAR_MATCH_TYPE )  {
            match += line->cigars[i] >> GENEPRED_CIGAR_PACKED_FIELD;
        }
        else if ( line->cigars[i] & GENEPRED_CIGAR_DELETE_TYPE ) {

            del = line->cigars[i] >> GENEPRED_CIGAR_PACKED_FIELD;
 
            // there is no need to check match and *pos, becase for plus strand match will be always smaller than *pos in  this function,
            // check if this deletion in the target block
            if ( line->loc[line->strand == '+' ? BLOCK_START : BLOCK_END][i] <= match && *pos > match) {
                adjust -= del;            
                // if this variant located in the deletion, just put pos to the edge of this gap
                if ( *pos < match +  del ) {
                    *pos += 1;
                    return 0;
                }
            }
            match += del;
        }
        else if ( line->cigars[i] & GENEPRED_CIGAR_INSERT_TYPE ) {
            ins = line->cigars[i] >> GENEPRED_CIGAR_PACKED_FIELD;
            if ( line->loc[line->strand == '+' ? BLOCK_START : BLOCK_END][block1] <= match  && *pos > match )
                adjust += ins;
        }
        if ( line->strand == '+' && match >= *pos )     
            break;
        else if ( line->strand == '-' && match < *pos )
            break;
    }
    *pos += adjust;

    return 0;
}

static int convert_transLite(int pos, struct genepred_line *gp, struct transLite *trans) {

    memset(trans, 0, sizeof(struct transLite));

    if ( genepred_find_location(gp, &trans->pos, &trans->offset, pos, &trans->exon_id) )
        return 1;
    
    return 0;
}
// pos is 1 based
static struct genepred_line *genepred_retrieve_target(char *name, int pos)
{
    int id;
    id = tbx_name2id(args.idx_data, name);
    if ( id == -1 )
        return NULL;
    
    hts_itr_t *itr = tbx_itr_queryi(args.idx_data, id, pos-1, pos);
    if ( itr == NULL )
        return NULL;
    
    kstring_t string = {0, 0, 0};
    struct genepred_line *head = NULL;
    struct genepred_line *temp = NULL;
    
    while ( tbx_itr_next(args.fp_data, args.idx_data, itr, &string) >= 0 ) {
        struct genepred_line *line = genepred_line_create();
        if ( parse_line(&string, line) )
	    continue;
        
        if ( head == NULL )
            head = line;
        
        if ( temp )
            temp->next = line;

        temp = line;
        string.l = 0;
    }  
    free(string.s);
    tbx_itr_destroy(itr);
    return head;
}

struct bedHandler {
    FILE                 *fp_input;
    FILE                 *fp_output;
    hts_itr_t            *itr;
    int                   n_lines;
    struct bedDetail      bed;
    struct IDdata         IDdata;
    struct genepred_line *start;
    struct genepred_line *end;
} handler = {
    .fp_input  = NULL,
    .fp_output = NULL,
    .itr       = NULL,
    .n_lines   = 1,
    .bed       = { NULL, 0, 0, NULL, NULL, NULL },
    .IDdata    = { NULL, 0, 0, NULL },
};

int init_handler ()
{
    // default is stdin
    handler.fp_input = strcmp(args.input_fname, "-") == 0 ? stdin : fopen(args.input_fname, "r");
    if ( handler.fp_input == NULL )
        error("%s : %s", args.input_fname, strerror(errno));

    handler.fp_output = stdout;

    // todo: check the first line if track information specified
    return 0;
}

void clean_bed(struct bedDetail *bed)
{
    if ( bed->chrom )
        free(bed->chrom);

    if ( bed->name )
        free(bed->name);

    if ( bed->id )
        free(bed->id);

    if ( bed->description )
        free(bed->description);
                    
    memset(bed, 0, sizeof(struct bedDetail));    
}

void clean_handler()
{
    clean_bed(&handler.bed);
    clean_IDdata(&handler.IDdata);
    
    if ( handler.itr )
        tbx_itr_destroy(handler.itr);

    handler.itr = NULL;
    
    if ( handler.start )
        list_lite_del(&handler.start, genepred_line_destroy);

    handler.start = NULL;
    
    if ( handler.end )
        list_lite_del(&handler.end, genepred_line_destroy);
    handler.end = NULL;
    
}
int fill_handler()
{
    if ( feof(handler.fp_input ))
        return 1;

    clean_handler();
    
    kstring_t str = { 0, 0, 0};
    int       col = 0;
    
    struct bedDetail *bed = &handler.bed;
    
    for ( ;; ) {
        // end of line
        if ( feof(handler.fp_input) )
            return 1;
        
        char c = fgetc(handler.fp_input);
        
#define BRANCH(_key, _func, _len) do {                                       \
            if ( c == '\t' || c == '\n' || feof(handler.fp_input) ) {   \
                if ( _len == 0 ) break;\
                _key = _func(str.s, _len);                                   \
                str.l = 0;\
                col++;\
            }\
            else {\
                kputc(c, &str);\
            }\
        } while(0)\
            
        if ( col == 0 ) {
            BRANCH(bed->chrom, strndup, str.l);
        }
        else if ( col == 1 ) {
            BRANCH(bed->start, str2int_l, str.l);
        }
        else if ( col == 2 ) {
            BRANCH(bed->end, str2int_l, str.l);
        }
        else { 
            for ( ; c != '\n' && !feof(handler.fp_input); )
                c = fgetc(handler.fp_input);
            str.l = 0;                
            break;
        }

#undef BRANCH
    }
    // if truncated return -1
    if ( bed->end == 0 ) {
        bed->end = bed->start;
        bed->start--;
        if ( bed->start < 0 ) {
            warnings("Truncated line. %d ", handler.n_lines);
            return -1;
        }
    }

    handler.n_lines ++;
    // clean memory
    if ( str.m ) free(str.s);
    
    return 0;
}

void output_handler()
{
    kstring_t str  = { 0, 0, 0 };
    kstring_t name = { 0, 0, 0 };
    kstring_t des  = { 0, 0, 0 };

    struct bedDetail *bed = &handler.bed;
    struct transLite transLite;
    
    if ( handler.itr != NULL ) {
        for ( ;; ) {
            if ( tbx_itr_next(args.fp_name, args.idx_name, handler.itr, &str) < 0)
                break;
            
            parse_IDdata(str.s, str.l, &handler.IDdata);
            
            if ( bed->start >= handler.IDdata.start && bed->end <= handler.IDdata.end ) {
                if ( name.l ) kputc(',', &name);
                kputs(handler.IDdata.extra, &name);
            }
        }
    }

    if ( handler.start == NULL ) {
        kputs("Start=Intergenic;", &des);
    }
    else {
        kputs("Start=", &des);
        struct genepred_line *line = handler.start;
        for ( ;; ) {
            if ( line == NULL ) break;
            convert_transLite(bed->start, line, &transLite);
            if ( transLite.offset > 0 )
                ksprintf(&des, "%s(%s):n.%d+%d,EX%d|", line->name1, line->name2, transLite.pos, transLite.offset, transLite.exon_id+1);
            else if ( transLite.offset < 0 )
                ksprintf(&des, "%s(%s):n.%d%d,EX%d|", line->name1, line->name2, transLite.pos, transLite.offset, transLite.exon_id+1);
            else
                ksprintf(&des, "%s(%s):n.%d,EX%d|", line->name1, line->name2, transLite.pos, transLite.exon_id+1);
            line = line->next;
        }
        kputc(';', &des);
    }
    
    if ( handler.end == NULL ) {
        kputs("End=Intergenic;", &des);
    }
    else {
        kputs("End=", &des);
        struct genepred_line *line = handler.end;
        for ( ;; ) {
            if ( line == NULL ) break;
            convert_transLite(bed->end, line, &transLite);
            if ( transLite.offset > 0 )
                ksprintf(&des, "%s(%s):n.%d+%d,EX%d|", line->name1, line->name2, transLite.pos, transLite.offset, transLite.exon_id+1);
            else if ( transLite.offset < 0 )
                ksprintf(&des, "%s(%s):n.%d%d,EX%d|", line->name1, line->name2, transLite.pos, transLite.offset, transLite.exon_id+1);
            else
                ksprintf(&des, "%s(%s):n.%d,EX%d|", line->name1, line->name2, transLite.pos, transLite.exon_id+1);

            line = line->next;
        }
        kputc(';', &des);        
    }

    if ( name.l )
        handler.bed.name = name.s;

    if ( des.l )
        handler.bed.description = des.s;

    if ( str.m)
        free(str.s);
    
    fprintf(handler.fp_output, "%s\t%d\t%d\t%s\t%s\t%s\n",handler.bed.chrom, handler.bed.start, handler.bed.end,
            handler.bed.name == NULL ? "." : handler.bed.name,
            handler.bed.id   == NULL ? "." : handler.bed.id,
            handler.bed.description == NULL ? "." : handler.bed.description);
}

// Input format:
//  chrom \t start \t end 
int parse_args(int argc, char **argv)
{
    if ( argc == 1 )
        return usage();

    int i;
    for ( i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;

        if ( strcmp(a, "-data") == 0 )
            var = & args.data_fname;
        else if ( strcmp(a, "-name") == 0 )
            var = & args.name_fname;

        if ( var != 0 ) {
            if ( argc == i )
                error("Missing argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if ( args.input_fname == 0 ) {
            args.input_fname = a;
            continue;
        }

        error("Unknown argument, %s.", a);
    }

    if ( args.input_fname == NULL && !isatty(fileno((FILE*)stdin)) )
        args.input_fname = "-";

    if ( args.input_fname == NULL )
        error("No input bed file specified");

    if ( args.data_fname == NULL )
        error("genepredPlus database must be specified with -data.");

    if ( args.name_fname == NULL )
        error("Four-columns name database must specified with -name.");


    args.fp_data = hts_open(args.data_fname, "r");
    if ( args.fp_data == NULL )
        error("%s : %s.", args.data_fname, strerror(errno));

    args.fp_name = hts_open(args.name_fname, "r");
    if ( args.fp_name == NULL )
        error("%s : %s.", args.name_fname, strerror(errno));

    args.idx_name = tbx_index_load(args.name_fname);
    if ( args.idx_name == NULL )
        error("Failed to load index file of %s.", args.name_fname);
    
    args.idx_data = tbx_index_load(args.data_fname);
    if ( args.idx_data == NULL )
        error("Failed to load index file of %s.", args.data_fname);
    
    return 0;
}

void output_header()
{
    if ( args.UCSC_bedDetails_flag == 1 ) {
        fprintf(handler.fp_output, "track name=%s type=bedDetail description=\"%s\" db=%s", args.track_name, args.description, args.db_version);
    }
    else {
        fprintf(handler.fp_output, "#filetype=bedDetail\n##track_name=%s\n##description=%s\n##db=%s\n", args.track_name, args.description, args.db_version);        
        fprintf(handler.fp_output, "#chrom\tstart\tend\tname\tid\tdescription\n");
    }
    
}

int anno_bed()
{
    int n_lines = 0;
    int ret;
    init_handler();

    struct bedDetail *bed = &handler.bed;
    
    set_format_genepredPlus();

           
    for ( ;; ) {
        
        ret = fill_handler();

        // Truncated line.
        if ( ret == -1 )
            continue;
        // End of file.
        else if ( ret == 1 )
            break;
        
        // id, set only if region length is smaller than 100bp
        if ( bed->end - bed->start < 100 ) {
            int tid = tbx_name2id(args.idx_name, bed->chrom);
            //temp.l = 0;
            
            if ( tid > -1 )
                handler.itr = tbx_itr_queryi(args.idx_name, tid, bed->start, bed->end);
        }
        // convert location information to description column
        handler.start = genepred_retrieve_target(bed->chrom, bed->start);

        if ( bed->end > bed->start +1 ) 
            handler.end = genepred_retrieve_target(bed->chrom, bed->end);
        
        // id
        if ( bed->id == NULL ) {
            kstring_t str1 = {0, 0, 0};
            kputw(n_lines, &str1);
            bed->id = strdup(str1.s);
            free(str1.s);
        }

        output_handler();

        // increase line count
        n_lines++;
    }

    // empty bed file
    if ( n_lines == 0 )
        return 1;
    
    return 0;
}

int clean_memory()
{
    clean_bed(&handler.bed);
    fclose(handler.fp_input);
    fclose(handler.fp_output);
    tbx_destroy(args.idx_data);
    tbx_destroy(args.idx_name);
    hts_close(args.fp_data);
    hts_close(args.fp_name);
    return 0;
}

// Annotation results:

int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;
    
    if ( anno_bed() )
        return 1;

    clean_memory();

    return 0;
}

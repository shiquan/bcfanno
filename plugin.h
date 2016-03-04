#ifndef VCFANNO_PLUGIN_HEADER
#define VCFANNO_PLUGIN_HEADER

typedef struct {
    char **cols;
    int ncols, mcols;
    char **als;
    kstring_t line;
    int rid, start, end;
} annot_line_t;


    


#endif

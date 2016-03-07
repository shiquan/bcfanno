#ifndef VCFANNO_PLUGIN_HEADER
#define VCFANNO_PLUGIN_HEADER

/* try to convert any data from SQL into annot_line struct and annotated by core_annotate()
 */
typedef struct {
    char **cols;
    int ncols, mcols;
    char **als;
    kstring_t line;
    int rid, start, end;
} annot_line_t;


    


#endif

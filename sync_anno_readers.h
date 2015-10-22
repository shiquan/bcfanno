/*   anno.h -- annotations file type and data structure APIs

     (C) shiquan@genomics.cn

     @ MIT
 */
#ifndef _VCF_ANNO_H
#define _VCF_ANNO_H

#define REPLACE_MISSING  0 // replace only missing values
#define REPLACE_ALL      1 // replace both missing and existing values
#define REPLACE_EXISTING 2 // replace only if tgt is not missing

struct _args_t;

typedef struct _anno_col_t
{
    int icol;
    int replace; // REPLACE_*
    int number; // one of BCF_VL_* types
    char *hdr_key;
    int (*setter)(struct _args_t *, bcf1_t *, struct _anno_col_t *, void *);
    
}
anno_col_t;

typedef struct
{
}
anno_hdr_t;

typedef struct
{
    char **cols;
    int n_col, m_col;
    char **als;
    int n_als, m_als;
    kstring_t line;
    int rid, start, end;
}
anno_line_t;

typedef struct
{
    anno_hdr_t *header;
    const char *fname;
    anno_line_t **buffer;
    int nbuffer, mbuffer;
}
anno_sr_t;

typedef enum
{
    open_failed, connect_error, api_usage_error, header_error
}
anno_sr_error;

typedef enum
{
    anno_is_vcf, anno_is_dynlib, anno_is_tbx,
    anno_unsupport
}
anno_file_type;

typedef struct
{
    int *has_line;
    anno_sr_error errnum;

    anno_sr_t *readers;
    int nreaders;
    int streaming; // reading mode, SQL or streaming
}
anno_srs_t;

typedef struct _anno_api
{
    anno_file_type anno_type;
    int isrc;
    char *file;
    char *dynlib_path;
    anno_col_t *cols;
    int ncols;
    
}
anno_api_t;

#ifdef __cplusplus
extern "C" {
#endif
    /** Init and destroy anno_srs_t struct **/
    anno_srs_t * anno_sr_init(void);
    void anno_sr_destroy(anno_srs_t *readers);

    char *anno_sr_strerror(int errnum);

    /**
     * anno_sr_add_reader() - open new reader api
     * @readers: holder of the open readers
     * @fname: the dataset file
     *
     * Returns 1 if the call successed, or 0 on error.
     */
    int anno_sr_add_reader(anno_srs_t *readers, const char **fname);
    void anno_sr_remove_reader(anno_srs_t *readers, int ir);

    /**
     * anno_sr_next_line() - the iterator
     * @readers: holder of the open readers
     *
     * Returns the number of readers which have the current line
     * (anno_sr_t.buffer[0]) set at this position. Use the anno_sr_has_line macro to
     * determine which of the readers are set.
     */
    int anno_sr_next_line(anno_srs_t *readers);
#define anno_sr_has_line(readers, i) (readers)->has_line[i]
#define anno_sr_get_line(readers, i) ((readers)->has_line[i] ? ((readers)->readers[i].buffer[0]) : NULL)
#define anno_sr_region_done(readers, i) (!(readers)->has_line[i] && !(readers)->readers[i].nbuffer ? 1 : 0)
#define anno_sr_get_header(readers, i) (readers)->readers[i].header
#define anno_sr_get_reader(readers, i) &((readers)->readers[i])

    /**
     * anno_sr_seek() - set all readers to selected position
     * @seq: sequence name; NULL to seek to start
     * @pos: 0-based coordinate
     */
    int anno_sr_seek(anno_srs_t *readers, const char *seq, int pos);

    
#ifdef __cplusplus
}
#endif

#endif

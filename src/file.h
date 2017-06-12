#ifndef FILE_HEADER
#define FILE_HEADER

#include <stdio.h>
#include <stdlib.h>

// The stream_lite structure is defined for fast access htsFile.
/* struct stream_lite { */
/*     int begin, end; */
/*     int is_eof:2, bufsize:30; */
/*     uint64_t seek_pos; */
/*     BGZF* f; */
/*     unsigned char *buf; */
/* }; */

/* extern int file_seek(htsFile *fp, long offset, int where); */

extern size_t fread_noeof(void *ptr, size_t size, size_t nmemb, FILE *stream);

extern long fread_safe(FILE *fp, long size, void *a);

#endif

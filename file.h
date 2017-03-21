#ifndef FILE_HEADER
#define FILE_HEADER

#include <stdio.h>
#include <stdlib.h>
#include "htslib/hts.h"
#include "htslib/bgzf.h"

// The stream_lite structure is defined for fast access htsFile.
struct stream_lite {
    int begin, end;
    int is_eof:2, bufsize:30;
    uint64_t seek_pos;
    BGZF* f;
    unsigned char *buf;
};

extern int file_seek(htsFile *fp, long offset, int where);

#endif

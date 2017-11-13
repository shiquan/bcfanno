/*  
    Copyright (C) 2016,2017  BGI Research

    Author: Shi Quan (shiquan@genomics.cn)

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

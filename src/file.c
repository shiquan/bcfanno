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

#include "utils.h"
#include "file.h"
#include "htslib/kseq.h"

size_t fread_noeof(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
    size_t ret = fread(ptr, size, nmemb, stream);
    if (ret != nmemb )
        error("%s", ferror(stream) ? strerror(errno) : "Unexpected end of file");
    return ret;
}

// Mac / Darwim has a bug when reading data longer than 2GB. This function fixes this issue
// by reading data in small chunks
long fread_safe(FILE *fp, long size, void *a)
{
    const int bufsize = 0x1000000; // 16M block
    long offset = 0;
    while ( size ) {
        int x = bufsize < size ? bufsize : size;
        if ( ( x = fread_noeof(a + offset, 1, x, fp) ) == 0 )
            break;
        size -= x;
        offset += x;
    }
    return offset;
}

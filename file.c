#include "utils.h"
#include "file.h"
#include "htslib/kseq.h"

/* int file_seek(htsFile *fp, long offset, int where) */
/* { */
/*     if ( fp->is_bin ) { */
/*         return bgzf_seek(fp->fp.bgzf, offset, where); */
/*     } else { */
/*         ks_rewind((struct stream_lite*)fp->fp.voidp); */
/*         ((struct stream_lite*)fp->fp.voidp)->seek_pos = offset; */
/*         return bgzf_seek(((struct stream_lite*)fp->fp.voidp)->f, offset, where); */
/*     } */
/* } */

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

#include "utils.h"
#include "file.h"
#include "htslib/kseq.h"

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

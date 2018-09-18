// sequence.c - handle ACGT sequences

#include "utils.h"
#include "motif_encode.h"

// 1 on failure, 0 on success
int is_atcg(char *s)
{
    int i,l;
    l = strlen(s);
    if ( l == 0 ) return 1;
    for ( i = 0; i < l; ++i )
        if ( _enc[s[i]] == 0 ) return 1;
    return 0;
}

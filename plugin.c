#include "plugin.h"

void *dlopen_plugin(const char *fname)
{
    void *handle = NULL;
    handle = dlopen(fname, RTLD_NOW);
    if ( !handle ) {
	error("%s:\n\tdlopen .. %s\n", fname, dlerror());
    }
    return handle;
}



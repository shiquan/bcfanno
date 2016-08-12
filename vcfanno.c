#include "utils.h"
#include "anno.h"
#include "version.h"

#ifndef KSTRING_INIT
#define KSTRING_INIT {0, 0, 0}
#endif

const char *config_file = NULL;
const char *input_fname = NULL;
const char *output_fname = NULL;
int silence_mode = 0;
int file_type = 0;

static void prase_argv(int argc, char **argv)
{
    char *avs = NULL;
    kstring_t str = KSTRING_INIT;
    int i;
    

    if (silence_mode == 0) {
	LOG_print("Compile platform: %s, complie time: %s", COMPLIE_OS, COMPLIE_TIME);
	LOG_print("Args: %s", avs);
    }
}

int main(int argc, char **argv)
{
    prase_argv(argc, argv);
    
}

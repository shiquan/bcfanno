#include "utils.h"
#include "anno.h"
#include "version.h"

#ifndef KSTRING_INIT
#define KSTRING_INIT {0, 0, 0}
#endif

struct args {
    const char *config_file;
    const char *input_fname;
    const char *output_fname;
    
    int silence_mode;
    int file_type;
    int test_mode;
};

struct args args = {
    .config_file = 0,
    .input_fname = 0,
    .output_fname = 0,
    .silence_mode = 0,
    .file_type = 0,
    .test_mode = 0,
};

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

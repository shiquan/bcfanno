#ifndef OPTIONS_HEADER
#define OPTIONS_HEADER
#include "khash.h"

// option type
enum option_type {
    option_boolean,
    option_string,
    option_int,
    option_float,
    option_long,
    option_double,
};

// Sepcification of options. Consist of an array of a single option are passed to
// options_init() to validate options.
struct options_spec {
    char *name;
    unsigned int flags;
};

struct options_compact {
};
int options_init(int argc, char *argv[], char *avs);
int option_exist(char *opt);
char *get_option(int id);

#endif

PROG=       vcfanno
DEBUG_PROG= vcfanno_debug

all: $(PROG)

debug: $(DEBUG_PROG)

HTSDIR = htslib-1.3
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a

CC       = gcc
CFLAGS   = -Wall -O3 -DHTS3
DEBUG_CFLAGS   = -g -Wall -O0 -DHTS3
DFLAGS   =
INCLUDES = -I. -I$(HTSDIR)/


all:$(PROG)

# See htslib/Makefile
PACKAGE_VERSION = 0.01
ifneq "$(wildcard .git)" ""
PACKAGE_VERSION := $(shell git describe --always --dirty)
DOC_VERSION :=  $(shell git describe --always)+
DOC_DATE := $(shell date +'%Y-%m-%d %R %Z')
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))
endif
version.h:
	echo '#define VCFANNO_VERSION "$(PACKAGE_VERSION)"' > $@


.SUFFIXES:.c .o
.PHONY:all clean clean-all clean-plugins distclean install lib tags test testclean force plugins docs

force:

.c.o:
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@


# Plugin rules
#PLUGINC = $(foreach dir, plugins, $(wildcard $(dir)/*.c))
#PLUGINS = $(PLUGINC:.c=.so)
#PLUGINM = $(PLUGINC:.c=.mk)

#%.so: %.c version.h version.c $(HTSDIR)/libhts.so
#	$(CC) $(CFLAGS) $(INCLUDES) -fPIC -shared -o $@ version.c $< -L$(HTSDIR) -lhts

#-include $(PLUGINM)

#plugins: $(PLUGINS)

#hgvs_generate: $(HTSLIB) 
#	$(CC) -D_HGVS_MAIN $(DEBUG_CFLAGS) $(INCLUDES) -pthread -lz -o $@ sequence.c genepred.c hgvs_generate.c $(HTSLIB)

vcfadd: $(HTSLIB) 
	$(CC) -D_VCF_ANNOS_MAIN $(DEBUG_CFLAGS) $(INCLUDES) -pthread -lz -o $@ vcf_annos.c json_config.c config.c kson.c vcmp.c $(HTSLIB)

bedadd: $(HTSLIB) 
	$(CC) -D_BED_ANNOS_MAIN $(DEBUG_CFLAGS) $(INCLUDES) -pthread -lz -o $@ anno_bed.c json_config.c config.c kson.c $(HTSLIB)

vcf2tsv: $(HTSLIB) version.h 
	$(CC) $(CFLAGS) $(INCLUDES) -pthread -lz -o $@ misc/vcf2tsv.c $(HTSLIB)

tsv2vcf: $(HTSLIB) version.h
	$(CC) $(CFLAGS) $(INCLUDES) -pthread -lz -o $@ misc/tsv2vcf.c $(HTSLIB)

vcf_rename_tags: $(HTSLIB) version.h 
	$(CC) $(CFLAGS) $(INCLUDES) -pthread -lz -o $@ misc/vcf_rename_tags.c $(HTSLIB)

vcfanno: $(HTSLIB) version.h 
	$(CC) $(CFLAGS) $(INCLUDES) -lz -pthread -o $@ vcfanno_main.c vcf_annos.c anno_bed.c sequence.c genepred.c hgvs.c hgvs_vcf.c number.c vcmp.c json_config.c config.c kson.c file.c sort_list.c $(HTSLIB)

vcfanno_debug: $(HTSLIB) version.h
	$(CC) -DDEBUG_MODE $(DEBUG_CFLAGS) $(INCLUDES) -lz -pthread -o $@ vcfanno_main.c vcf_annos.c anno_bed.c sequence.c genepred.c hgvs.c hgvs_vcf.c number.c vcmp.c json_config.c config.c kson.c file.c sort_list.c $(HTSLIB)

test: $(HTSLIB) version.h vcfadd bedadd

clean: testclean
	-rm -f gmon.out *.o *~ $(PROG) version.h 
	-rm -rf *.dSYM plugins/*.dSYM test/*.dSYM
	-rm -f anno_vcf bedadd vcfadd vcfanno anno_bed hgvs_generate
	-rm -f config vcfanno_debug vcf2tsv tsv2vcf vcf_rename_tags

testclean:
	-rm -f test/*.o test/*~ $(TEST_PROG)

distclean: clean
	-rm -f TAGS

clean-all: clean clean-htslib

tags:
	ctags -f TAGS *.[ch] plugins/*.[ch]

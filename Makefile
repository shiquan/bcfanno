PROG=       vcfanno
##TEST_PROG=  test/test-rbuf

all: $(PROG) $(TEST_PROG)

# Adjust $(HTSDIR) to point to your top-level htslib directory
HTSDIR = htslib-1.3
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
#BGZIP  = $(HTSDIR)/bgzip
#TABIX  = $(HTSDIR)/tabix

CC       = gcc
CFLAGS   = -g -Wall -Wc++-compat -O2
DFLAGS   =
#OBJS     = main.o vcfindex.o tabix.o \
#           vcfstats.o vcfisec.o vcfmerge.o vcfquery.o vcffilter.o filter.o vcfsom.o \
#           vcfnorm.o vcfgtcheck.o vcfview.o vcfannotate.o vcfroh.o vcfconcat.o \
#           vcfcall.o mcall.o vcmp.o gvcf.o reheader.o convert.o vcfconvert.o tsv2vcf.o \
#           vcfcnv.o HMM.o vcfplugin.o consensus.o ploidy.o version.o \
#           ccall.o em.o prob1.o kmin.o # the original samtools calling
INCLUDES = -I. -I$(HTSDIR)

# The polysomy command is not compiled by default because it brings dependency
# on libgsl. The command can be compiled wth `make USE_GPL=1`. See the INSTALL
# and LICENSE documents to understand license implications.

all:$(PROG) plugins

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

vcfanno: $(HTSLIB) version.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ anno_core.c vcmp.c config.c anno_setter.c kson.c vcfannotate.c $(HTSLIB) -lz


#docs: doc/bcftools.1 doc/bcftools.html

#install: $(PROG) doc/bcftools.1
#	$(INSTALL_DIR) $(DESTDIR)$(bindir) $(DESTDIR)$(man1dir)
#	$(INSTALL_PROGRAM) $(PROG) plot-vcfstats vcfutils.pl $(DESTDIR)$(bindir)
#	$(INSTALL_DATA) doc/bcftools.1 $(DESTDIR)$(man1dir)

clean: testclean
	-rm -f gmon.out *.o *~ $(PROG) version.h 
	-rm -rf *.dSYM plugins/*.dSYM test/*.dSYM


testclean:
	-rm -f test/*.o test/*~ $(TEST_PROG)

distclean: clean
	-rm -f TAGS

clean-all: clean clean-htslib

tags:
	ctags -f TAGS *.[ch] plugins/*.[ch]

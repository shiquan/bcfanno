#ifndef C_HGVS_VCF_HEADER
#define C_HGVS_VCF_HEADER

#include "utils.h"
#include "htslib/vcf.h"


extern int init_hgvs_anno(const char *data, const char *fasta, bcf_hdr_t *hdr);
extern int close_hgvs_anno();
extern int setter_hgvs_vcf(bcf_hdr_t *hdr, bcf1_t *line);






#endif 

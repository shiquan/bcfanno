// rename tags for vcf/bcf file, only support INFO fields
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include "utils.h"
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#define KSTRING_INIT {0, 0, 0}

static const char *list_fname = 0;

// only one flag could be set
static int convert_ucsc_name_flag = 0;
static int convert_ncbi_name_flag = 0;

struct args  {
    htsFile *fp_input;
    htsFile *fp_output;
    bcf_hdr_t *hdr_in;
    bcf_hdr_t *hdr_out;
    kstring_t command;
};

struct args args = {
    .fp_input = NULL,
    .fp_output = NULL,
    .hdr_in = NULL,
    .hdr_out = NULL,
    .command = KSTRING_INIT,
};

void args_destroy() {
    hts_close(args.fp_input);
    hts_close(args.fp_output);
    bcf_hdr_destroy(args.hdr_in);
    bcf_hdr_destroy(args.hdr_out);
    if ( args.command.m )
	free(args.command.s);
}
// tags list with two columns, old name -> new name

int usage() {
    fprintf(stderr, "vcf_renames_tags  input_file\n"
            "Options:\n"
            "  -list tags.txt     This file include two columns, they are old and new tags.\n"
            "  -ncbi              Convert UCSC chromosome name to NCBI type.\n"
            "  -ucsc              Convert NCBI chromosome name to UCSC type.\n"
            "  -O <type>          Output type.[bzu]\n"
            "  -o output_file     Output file.\n"
        );
    return 1;
}
static const char *hts_bcf_wmode(int file_type)
{
    if ( file_type == FT_BCF )
	return "wbu";    // uncompressed BCF
    if ( file_type & FT_BCF )
	return "wb";      // compressed BCF
    if ( file_type & FT_GZ )
	return "wz";       // compressed VCF
    return "w";                                 // uncompressed VCF
}
struct chr_pair {
    const char *ncbi;
    const char *ucsc;
};

static const struct chr_pair chr_pairs[] = {    
    { "1", "chr1", },
    { "2", "chr2", },
    { "3", "chr3", },
    { "4", "chr4", },
    { "5", "chr5", },
    { "6", "chr6", },
    { "7", "chr7", },
    { "8", "chr8", },
    { "9", "chr9", },
    { "10", "chr10", },
    { "11", "chr11", },
    { "12", "chr12", },
    { "13", "chr13", },
    { "14", "chr14", },
    { "15", "chr15", },
    { "16", "chr16", },
    { "17", "chr17", },
    { "18", "chr18", },
    { "19", "chr19", },
    { "20", "chr20", },
    { "21", "chr21", },
    { "22", "chr22", },
    { "X", "chrX", },
    { "Y", "chrY", },

    // mito chromosome
    { "MT", "chrM", },
    { "M", "chrM", },

    // alternative locus, unplaced contigs
    // hg19
    { "GL000191.1", "chr1_gl000191_random", },
    { "GL000192.1", "chr1_gl000192_random", },
    { "GL000193.1", "chr4_gl000193_random", },
    { "GL000194.1", "chr4_gl000194_random", },
    { "GL000195.1", "chr7_gl000195_random", },
    { "GL000196.1", "chr8_gl000196_random", },
    { "GL000197.1", "chr8_gl000197_random", },
    { "GL000198.1", "chr9_gl000198_random", },
    { "GL000199.1", "chr9_gl000199_random", },
    { "GL000200.1", "chr9_gl000200_random", },
    { "GL000201.1", "chr9_gl000201_random", },
    { "GL000202.1", "chr11_gl000202_random", },
    { "GL000203.1", "chr17_gl000203_random", },
    { "GL000204.1", "chr17_gl000204_random", },
    { "GL000205.1", "chr17_gl000205_random", },
    { "GL000206.1", "chr17_gl000206_random", },
    { "GL000207.1", "chr18_gl000207_random", },
    { "GL000208.1", "chr19_gl000208_random", },
    { "GL000209.1", "chr19_gl000209_random", },
    { "GL000210.1", "chr21_gl000210_random", },
    { "GL000211.1", "chrUn_gl000211", },
    { "GL000212.1", "chrUn_gl000212", },
    { "GL000213.1", "chrUn_gl000213", },
    { "GL000214.1", "chrUn_gl000214", },
    { "GL000215.1", "chrUn_gl000215", },
    { "GL000216.1", "chrUn_gl000216", },
    { "GL000217.1", "chrUn_gl000217", },
    { "GL000218.1", "chrUn_gl000218", },
    { "GL000219.1", "chrUn_gl000219", },
    { "GL000220.1", "chrUn_gl000220", },
    { "GL000221.1", "chrUn_gl000221", },
    { "GL000222.1", "chrUn_gl000222", },
    { "GL000223.1", "chrUn_gl000223", },
    { "GL000224.1", "chrUn_gl000224", },
    { "GL000225.1", "chrUn_gl000225", },
    { "GL000226.1", "chrUn_gl000226", },
    { "GL000227.1", "chrUn_gl000227", },
    { "GL000228.1", "chrUn_gl000228", },
    { "GL000229.1", "chrUn_gl000229", },
    { "GL000230.1", "chrUn_gl000230", },
    { "GL000231.1", "chrUn_gl000231", },
    { "GL000232.1", "chrUn_gl000232", },
    { "GL000233.1", "chrUn_gl000233", },
    { "GL000234.1", "chrUn_gl000234", },
    { "GL000235.1", "chrUn_gl000235", },
    { "GL000236.1", "chrUn_gl000236", },
    { "GL000237.1", "chrUn_gl000237", },
    { "GL000238.1", "chrUn_gl000238", },
    { "GL000239.1", "chrUn_gl000239", },
    { "GL000240.1", "chrUn_gl000240", },
    { "GL000241.1", "chrUn_gl000241", },
    { "GL000242.1", "chrUn_gl000242", },
    { "GL000243.1", "chrUn_gl000243", },
    { "GL000244.1", "chrUn_gl000244", },
    { "GL000245.1", "chrUn_gl000245", },
    { "GL000246.1", "chrUn_gl000246", },
    { "GL000247.1", "chrUn_gl000247", },
    { "GL000248.1", "chrUn_gl000248", },
    { "GL000249.1", "chrUn_gl000249", },
    { "GL000250.1", "chr6_apd_hap1", },
    { "GL000251.1", "chr6_cox_hap2", },
    { "GL000252.1", "chr6_dbb_hap3", },
    { "GL000253.1", "chr6_mann_hap4", },
    { "GL000254.1", "chr6_mcf_hap5", },
    { "GL000255.1", "chr6_qbl_hap6", },
    { "GL000256.1", "chr6_ssto_hap7", },
    { "GL000257.1", "chr4_ctg9_hap1", },
    { "GL000258.1", "chr17_ctg5_hap1", },

    // hg38
    { "KI270706.1", "chr1_KI270706v1_random", },
    { "KI270707.1", "chr1_KI270707v1_random", },
    { "KI270708.1", "chr1_KI270708v1_random", },
    { "KI270709.1", "chr1_KI270709v1_random", },
    { "KI270710.1", "chr1_KI270710v1_random", },
    { "KI270711.1", "chr1_KI270711v1_random", },
    { "KI270712.1", "chr1_KI270712v1_random", },
    { "KI270713.1", "chr1_KI270713v1_random", },
    { "KI270714.1", "chr1_KI270714v1_random", },
    { "KI270715.1", "chr2_KI270715v1_random", },
    { "KI270716.1", "chr2_KI270716v1_random", },
    { "GL000221.1", "chr3_GL000221v1_random", },
    { "GL000008.2", "chr4_GL000008v2_random", },
    { "GL000208.1", "chr5_GL000208v1_random", },
    { "KI270717.1", "chr9_KI270717v1_random", },
    { "KI270718.1", "chr9_KI270718v1_random", },
    { "KI270719.1", "chr9_KI270719v1_random", },
    { "KI270720.1", "chr9_KI270720v1_random", },
    { "KI270721.1", "chr11_KI270721v1_random", },
    { "GL000009.2", "chr14_GL000009v2_random", },
    { "GL000194.1", "chr14_GL000194v1_random", },
    { "GL000225.1", "chr14_GL000225v1_random", },
    { "KI270722.1", "chr14_KI270722v1_random", },
    { "KI270723.1", "chr14_KI270723v1_random", },
    { "KI270724.1", "chr14_KI270724v1_random", },
    { "KI270725.1", "chr14_KI270725v1_random", },
    { "KI270726.1", "chr14_KI270726v1_random", },
    { "KI270727.1", "chr15_KI270727v1_random", },
    { "KI270728.1", "chr16_KI270728v1_random", },
    { "GL000205.2", "chr17_GL000205v2_random", },
    { "KI270729.1", "chr17_KI270729v1_random", },
    { "KI270730.1", "chr17_KI270730v1_random", },
    { "KI270731.1", "chr22_KI270731v1_random", },
    { "KI270732.1", "chr22_KI270732v1_random", },
    { "KI270733.1", "chr22_KI270733v1_random", },
    { "KI270734.1", "chr22_KI270734v1_random", },
    { "KI270735.1", "chr22_KI270735v1_random", },
    { "KI270736.1", "chr22_KI270736v1_random", },
    { "KI270737.1", "chr22_KI270737v1_random", },
    { "KI270738.1", "chr22_KI270738v1_random", },
    { "KI270739.1", "chr22_KI270739v1_random", },
    { "KI270740.1", "chrY_KI270740v1_random", },
    { "GL000195.1", "chrUn_GL000195v1", },
    { "GL000213.1", "chrUn_GL000213v1", },
    { "GL000214.1", "chrUn_GL000214v1", },
    { "GL000216.2", "chrUn_GL000216v2", },
    { "GL000218.1", "chrUn_GL000218v1", },
    { "GL000219.1", "chrUn_GL000219v1", },
    { "GL000220.1", "chrUn_GL000220v1", },
    { "GL000224.1", "chrUn_GL000224v1", },
    { "GL000226.1", "chrUn_GL000226v1", },
    { "KI270302.1", "chrUn_KI270302v1", },
    { "KI270303.1", "chrUn_KI270303v1", },
    { "KI270304.1", "chrUn_KI270304v1", },
    { "KI270305.1", "chrUn_KI270305v1", },
    { "KI270310.1", "chrUn_KI270310v1", },
    { "KI270311.1", "chrUn_KI270311v1", },
    { "KI270312.1", "chrUn_KI270312v1", },
    { "KI270315.1", "chrUn_KI270315v1", },
    { "KI270316.1", "chrUn_KI270316v1", },
    { "KI270317.1", "chrUn_KI270317v1", },
    { "KI270320.1", "chrUn_KI270320v1", },
    { "KI270322.1", "chrUn_KI270322v1", },
    { "KI270329.1", "chrUn_KI270329v1", },
    { "KI270330.1", "chrUn_KI270330v1", },
    { "KI270333.1", "chrUn_KI270333v1", },
    { "KI270334.1", "chrUn_KI270334v1", },
    { "KI270335.1", "chrUn_KI270335v1", },
    { "KI270336.1", "chrUn_KI270336v1", },
    { "KI270337.1", "chrUn_KI270337v1", },
    { "KI270338.1", "chrUn_KI270338v1", },
    { "KI270340.1", "chrUn_KI270340v1", },
    { "KI270362.1", "chrUn_KI270362v1", },
    { "KI270363.1", "chrUn_KI270363v1", },
    { "KI270364.1", "chrUn_KI270364v1", },
    { "KI270366.1", "chrUn_KI270366v1", },
    { "KI270371.1", "chrUn_KI270371v1", },
    { "KI270372.1", "chrUn_KI270372v1", },
    { "KI270373.1", "chrUn_KI270373v1", },
    { "KI270374.1", "chrUn_KI270374v1", },
    { "KI270375.1", "chrUn_KI270375v1", },
    { "KI270376.1", "chrUn_KI270376v1", },
    { "KI270378.1", "chrUn_KI270378v1", },
    { "KI270379.1", "chrUn_KI270379v1", },
    { "KI270381.1", "chrUn_KI270381v1", },
    { "KI270382.1", "chrUn_KI270382v1", },
    { "KI270383.1", "chrUn_KI270383v1", },
    { "KI270384.1", "chrUn_KI270384v1", },
    { "KI270385.1", "chrUn_KI270385v1", },
    { "KI270386.1", "chrUn_KI270386v1", },
    { "KI270387.1", "chrUn_KI270387v1", },
    { "KI270388.1", "chrUn_KI270388v1", },
    { "KI270389.1", "chrUn_KI270389v1", },
    { "KI270390.1", "chrUn_KI270390v1", },
    { "KI270391.1", "chrUn_KI270391v1", },
    { "KI270392.1", "chrUn_KI270392v1", },
    { "KI270393.1", "chrUn_KI270393v1", },
    { "KI270394.1", "chrUn_KI270394v1", },
    { "KI270395.1", "chrUn_KI270395v1", },
    { "KI270396.1", "chrUn_KI270396v1", },
    { "KI270411.1", "chrUn_KI270411v1", },
    { "KI270412.1", "chrUn_KI270412v1", },
    { "KI270414.1", "chrUn_KI270414v1", },
    { "KI270417.1", "chrUn_KI270417v1", },
    { "KI270418.1", "chrUn_KI270418v1", },
    { "KI270419.1", "chrUn_KI270419v1", },
    { "KI270420.1", "chrUn_KI270420v1", },
    { "KI270422.1", "chrUn_KI270422v1", },
    { "KI270423.1", "chrUn_KI270423v1", },
    { "KI270424.1", "chrUn_KI270424v1", },
    { "KI270425.1", "chrUn_KI270425v1", },
    { "KI270429.1", "chrUn_KI270429v1", },
    { "KI270435.1", "chrUn_KI270435v1", },
    { "KI270438.1", "chrUn_KI270438v1", },
    { "KI270442.1", "chrUn_KI270442v1", },
    { "KI270448.1", "chrUn_KI270448v1", },
    { "KI270465.1", "chrUn_KI270465v1", },
    { "KI270466.1", "chrUn_KI270466v1", },
    { "KI270467.1", "chrUn_KI270467v1", },
    { "KI270468.1", "chrUn_KI270468v1", },
    { "KI270507.1", "chrUn_KI270507v1", },
    { "KI270508.1", "chrUn_KI270508v1", },
    { "KI270509.1", "chrUn_KI270509v1", },
    { "KI270510.1", "chrUn_KI270510v1", },
    { "KI270511.1", "chrUn_KI270511v1", },
    { "KI270512.1", "chrUn_KI270512v1", },
    { "KI270515.1", "chrUn_KI270515v1", },
    { "KI270516.1", "chrUn_KI270516v1", },
    { "KI270517.1", "chrUn_KI270517v1", },
    { "KI270518.1", "chrUn_KI270518v1", },
    { "KI270519.1", "chrUn_KI270519v1", },
    { "KI270521.1", "chrUn_KI270521v1", },
    { "KI270522.1", "chrUn_KI270522v1", },
    { "KI270528.1", "chrUn_KI270528v1", },
    { "KI270529.1", "chrUn_KI270529v1", },
    { "KI270530.1", "chrUn_KI270530v1", },
    { "KI270538.1", "chrUn_KI270538v1", },
    { "KI270539.1", "chrUn_KI270539v1", },
    { "KI270544.1", "chrUn_KI270544v1", },
    { "KI270548.1", "chrUn_KI270548v1", },
    { "KI270579.1", "chrUn_KI270579v1", },
    { "KI270580.1", "chrUn_KI270580v1", },
    { "KI270581.1", "chrUn_KI270581v1", },
    { "KI270582.1", "chrUn_KI270582v1", },
    { "KI270583.1", "chrUn_KI270583v1", },
    { "KI270584.1", "chrUn_KI270584v1", },
    { "KI270587.1", "chrUn_KI270587v1", },
    { "KI270588.1", "chrUn_KI270588v1", },
    { "KI270589.1", "chrUn_KI270589v1", },
    { "KI270590.1", "chrUn_KI270590v1", },
    { "KI270591.1", "chrUn_KI270591v1", },
    { "KI270593.1", "chrUn_KI270593v1", },
    { "KI270741.1", "chrUn_KI270741v1", },
    { "KI270742.1", "chrUn_KI270742v1", },
    { "KI270743.1", "chrUn_KI270743v1", },
    { "KI270744.1", "chrUn_KI270744v1", },
    { "KI270745.1", "chrUn_KI270745v1", },
    { "KI270746.1", "chrUn_KI270746v1", },
    { "KI270747.1", "chrUn_KI270747v1", },
    { "KI270748.1", "chrUn_KI270748v1", },
    { "KI270749.1", "chrUn_KI270749v1", },
    { "KI270750.1", "chrUn_KI270750v1", },
    { "KI270751.1", "chrUn_KI270751v1", },
    { "KI270752.1", "chrUn_KI270752v1", },
    { "KI270753.1", "chrUn_KI270753v1", },
    { "KI270754.1", "chrUn_KI270754v1", },
    { "KI270755.1", "chrUn_KI270755v1", },
    { "KI270756.1", "chrUn_KI270756v1", },
    { "KI270757.1", "chrUn_KI270757v1", },
    { "GL383518.1", "chr1_GL383518v1_alt", },
    { "GL383519.1", "chr1_GL383519v1_alt", },
    { "GL383520.2", "chr1_GL383520v2_alt", },
    { "KI270759.1", "chr1_KI270759v1_alt", },
    { "KI270760.1", "chr1_KI270760v1_alt", },
    { "KI270761.1", "chr1_KI270761v1_alt", },
    { "KI270762.1", "chr1_KI270762v1_alt", },
    { "KI270763.1", "chr1_KI270763v1_alt", },
    { "KI270764.1", "chr1_KI270764v1_alt", },
    { "KI270765.1", "chr1_KI270765v1_alt", },
    { "KI270766.1", "chr1_KI270766v1_alt", },
    { "GL383521.1", "chr2_GL383521v1_alt", },
    { "GL383522.1", "chr2_GL383522v1_alt", },
    { "GL582966.2", "chr2_GL582966v2_alt", },
    { "KI270767.1", "chr2_KI270767v1_alt", },
    { "KI270768.1", "chr2_KI270768v1_alt", },
    { "KI270769.1", "chr2_KI270769v1_alt", },
    { "KI270770.1", "chr2_KI270770v1_alt", },
    { "KI270771.1", "chr2_KI270771v1_alt", },
    { "KI270772.1", "chr2_KI270772v1_alt", },
    { "KI270773.1", "chr2_KI270773v1_alt", },
    { "KI270774.1", "chr2_KI270774v1_alt", },
    { "KI270775.1", "chr2_KI270775v1_alt", },
    { "KI270776.1", "chr2_KI270776v1_alt", },
    { "GL383526.1", "chr3_GL383526v1_alt", },
    { "JH636055.2", "chr3_JH636055v2_alt", },
    { "KI270777.1", "chr3_KI270777v1_alt", },
    { "KI270778.1", "chr3_KI270778v1_alt", },
    { "KI270779.1", "chr3_KI270779v1_alt", },
    { "KI270780.1", "chr3_KI270780v1_alt", },
    { "KI270781.1", "chr3_KI270781v1_alt", },
    { "KI270782.1", "chr3_KI270782v1_alt", },
    { "KI270783.1", "chr3_KI270783v1_alt", },
    { "KI270784.1", "chr3_KI270784v1_alt", },
    { "GL000257.2", "chr4_GL000257v2_alt", },
    { "GL383527.1", "chr4_GL383527v1_alt", },
    { "GL383528.1", "chr4_GL383528v1_alt", },
    { "KI270785.1", "chr4_KI270785v1_alt", },
    { "KI270786.1", "chr4_KI270786v1_alt", },
    { "KI270787.1", "chr4_KI270787v1_alt", },
    { "KI270788.1", "chr4_KI270788v1_alt", },
    { "KI270789.1", "chr4_KI270789v1_alt", },
    { "KI270790.1", "chr4_KI270790v1_alt", },
    { "GL339449.2", "chr5_GL339449v2_alt", },
    { "GL383530.1", "chr5_GL383530v1_alt", },
    { "GL383531.1", "chr5_GL383531v1_alt", },
    { "GL383532.1", "chr5_GL383532v1_alt", },
    { "GL949742.1", "chr5_GL949742v1_alt", },
    { "KI270791.1", "chr5_KI270791v1_alt", },
    { "KI270792.1", "chr5_KI270792v1_alt", },
    { "KI270793.1", "chr5_KI270793v1_alt", },
    { "KI270794.1", "chr5_KI270794v1_alt", },
    { "KI270795.1", "chr5_KI270795v1_alt", },
    { "KI270796.1", "chr5_KI270796v1_alt", },
    { "GL000250.2", "chr6_GL000250v2_alt", },
    { "GL383533.1", "chr6_GL383533v1_alt", },
    { "KB021644.2", "chr6_KB021644v2_alt", },
    { "KI270797.1", "chr6_KI270797v1_alt", },
    { "KI270798.1", "chr6_KI270798v1_alt", },
    { "KI270799.1", "chr6_KI270799v1_alt", },
    { "KI270800.1", "chr6_KI270800v1_alt", },
    { "KI270801.1", "chr6_KI270801v1_alt", },
    { "KI270802.1", "chr6_KI270802v1_alt", },
    { "GL383534.2", "chr7_GL383534v2_alt", },
    { "KI270803.1", "chr7_KI270803v1_alt", },
    { "KI270804.1", "chr7_KI270804v1_alt", },
    { "KI270805.1", "chr7_KI270805v1_alt", },
    { "KI270806.1", "chr7_KI270806v1_alt", },
    { "KI270807.1", "chr7_KI270807v1_alt", },
    { "KI270808.1", "chr7_KI270808v1_alt", },
    { "KI270809.1", "chr7_KI270809v1_alt", },
    { "KI270810.1", "chr8_KI270810v1_alt", },
    { "KI270811.1", "chr8_KI270811v1_alt", },
    { "KI270812.1", "chr8_KI270812v1_alt", },
    { "KI270813.1", "chr8_KI270813v1_alt", },
    { "KI270814.1", "chr8_KI270814v1_alt", },
    { "KI270815.1", "chr8_KI270815v1_alt", },
    { "KI270816.1", "chr8_KI270816v1_alt", },
    { "KI270817.1", "chr8_KI270817v1_alt", },
    { "KI270818.1", "chr8_KI270818v1_alt", },
    { "KI270819.1", "chr8_KI270819v1_alt", },
    { "KI270820.1", "chr8_KI270820v1_alt", },
    { "KI270821.1", "chr8_KI270821v1_alt", },
    { "KI270822.1", "chr8_KI270822v1_alt", },
    { "GL383539.1", "chr9_GL383539v1_alt", },
    { "GL383540.1", "chr9_GL383540v1_alt", },
    { "GL383541.1", "chr9_GL383541v1_alt", },
    { "GL383542.1", "chr9_GL383542v1_alt", },
    { "KI270823.1", "chr9_KI270823v1_alt", },
    { "GL383545.1", "chr10_GL383545v1_alt", },
    { "GL383546.1", "chr10_GL383546v1_alt", },
    { "KI270824.1", "chr10_KI270824v1_alt", },
    { "KI270825.1", "chr10_KI270825v1_alt", },
    { "GL383547.1", "chr11_GL383547v1_alt", },
    { "JH159136.1", "chr11_JH159136v1_alt", },
    { "JH159137.1", "chr11_JH159137v1_alt", },
    { "KI270826.1", "chr11_KI270826v1_alt", },
    { "KI270827.1", "chr11_KI270827v1_alt", },
    { "KI270829.1", "chr11_KI270829v1_alt", },
    { "KI270830.1", "chr11_KI270830v1_alt", },
    { "KI270831.1", "chr11_KI270831v1_alt", },
    { "KI270832.1", "chr11_KI270832v1_alt", },
    { "GL383549.1", "chr12_GL383549v1_alt", },
    { "GL383550.2", "chr12_GL383550v2_alt", },
    { "GL383551.1", "chr12_GL383551v1_alt", },
    { "GL383552.1", "chr12_GL383552v1_alt", },
    { "GL383553.2", "chr12_GL383553v2_alt", },
    { "GL877875.1", "chr12_GL877875v1_alt", },
    { "GL877876.1", "chr12_GL877876v1_alt", },
    { "KI270833.1", "chr12_KI270833v1_alt", },
    { "KI270834.1", "chr12_KI270834v1_alt", },
    { "KI270835.1", "chr12_KI270835v1_alt", },
    { "KI270836.1", "chr12_KI270836v1_alt", },
    { "KI270837.1", "chr12_KI270837v1_alt", },
    { "KI270838.1", "chr13_KI270838v1_alt", },
    { "KI270839.1", "chr13_KI270839v1_alt", },
    { "KI270840.1", "chr13_KI270840v1_alt", },
    { "KI270841.1", "chr13_KI270841v1_alt", },
    { "KI270842.1", "chr13_KI270842v1_alt", },
    { "KI270843.1", "chr13_KI270843v1_alt", },
    { "KI270844.1", "chr14_KI270844v1_alt", },
    { "KI270845.1", "chr14_KI270845v1_alt", },
    { "KI270846.1", "chr14_KI270846v1_alt", },
    { "KI270847.1", "chr14_KI270847v1_alt", },
    { "GL383554.1", "chr15_GL383554v1_alt", },
    { "GL383555.2", "chr15_GL383555v2_alt", },
    { "KI270848.1", "chr15_KI270848v1_alt", },
    { "KI270849.1", "chr15_KI270849v1_alt", },
    { "KI270850.1", "chr15_KI270850v1_alt", },
    { "KI270851.1", "chr15_KI270851v1_alt", },
    { "KI270852.1", "chr15_KI270852v1_alt", },
    { "GL383556.1", "chr16_GL383556v1_alt", },
    { "GL383557.1", "chr16_GL383557v1_alt", },
    { "KI270853.1", "chr16_KI270853v1_alt", },
    { "KI270854.1", "chr16_KI270854v1_alt", },
    { "KI270855.1", "chr16_KI270855v1_alt", },
    { "KI270856.1", "chr16_KI270856v1_alt", },
    { "GL000258.2", "chr17_GL000258v2_alt", },
    { "GL383563.3", "chr17_GL383563v3_alt", },
    { "GL383564.2", "chr17_GL383564v2_alt", },
    { "GL383565.1", "chr17_GL383565v1_alt", },
    { "GL383566.1", "chr17_GL383566v1_alt", },
    { "JH159146.1", "chr17_JH159146v1_alt", },
    { "JH159147.1", "chr17_JH159147v1_alt", },
    { "KI270857.1", "chr17_KI270857v1_alt", },
    { "KI270858.1", "chr17_KI270858v1_alt", },
    { "KI270859.1", "chr17_KI270859v1_alt", },
    { "KI270860.1", "chr17_KI270860v1_alt", },
    { "KI270861.1", "chr17_KI270861v1_alt", },
    { "KI270862.1", "chr17_KI270862v1_alt", },
    { "GL383567.1", "chr18_GL383567v1_alt", },
    { "GL383568.1", "chr18_GL383568v1_alt", },
    { "GL383569.1", "chr18_GL383569v1_alt", },
    { "GL383570.1", "chr18_GL383570v1_alt", },
    { "GL383571.1", "chr18_GL383571v1_alt", },
    { "GL383572.1", "chr18_GL383572v1_alt", },
    { "KI270863.1", "chr18_KI270863v1_alt", },
    { "KI270864.1", "chr18_KI270864v1_alt", },
    { "GL383573.1", "chr19_GL383573v1_alt", },
    { "GL383574.1", "chr19_GL383574v1_alt", },
    { "GL383575.2", "chr19_GL383575v2_alt", },
    { "GL383576.1", "chr19_GL383576v1_alt", },
    { "GL949746.1", "chr19_GL949746v1_alt", },
    { "KI270865.1", "chr19_KI270865v1_alt", },
    { "KI270866.1", "chr19_KI270866v1_alt", },
    { "KI270867.1", "chr19_KI270867v1_alt", },
    { "KI270868.1", "chr19_KI270868v1_alt", },
    { "GL383577.2", "chr20_GL383577v2_alt", },
    { "KI270869.1", "chr20_KI270869v1_alt", },
    { "KI270870.1", "chr20_KI270870v1_alt", },
    { "KI270871.1", "chr20_KI270871v1_alt", },
    { "GL383578.2", "chr21_GL383578v2_alt", },
    { "GL383579.2", "chr21_GL383579v2_alt", },
    { "GL383580.2", "chr21_GL383580v2_alt", },
    { "GL383581.2", "chr21_GL383581v2_alt", },
    { "KI270872.1", "chr21_KI270872v1_alt", },
    { "KI270873.1", "chr21_KI270873v1_alt", },
    { "KI270874.1", "chr21_KI270874v1_alt", },
    { "GL383582.2", "chr22_GL383582v2_alt", },
    { "GL383583.2", "chr22_GL383583v2_alt", },
    { "KI270875.1", "chr22_KI270875v1_alt", },
    { "KI270876.1", "chr22_KI270876v1_alt", },
    { "KI270877.1", "chr22_KI270877v1_alt", },
    { "KI270878.1", "chr22_KI270878v1_alt", },
    { "KI270879.1", "chr22_KI270879v1_alt", },
    { "KI270880.1", "chrX_KI270880v1_alt", },
    { "KI270881.1", "chrX_KI270881v1_alt", },
    { "KI270892.1", "chr1_KI270892v1_alt", },
    { "KI270893.1", "chr2_KI270893v1_alt", },
    { "KI270894.1", "chr2_KI270894v1_alt", },
    { "KI270895.1", "chr3_KI270895v1_alt", },
    { "KI270896.1", "chr4_KI270896v1_alt", },
    { "KI270897.1", "chr5_KI270897v1_alt", },
    { "KI270898.1", "chr5_KI270898v1_alt", },
    { "GL000251.2", "chr6_GL000251v2_alt", },
    { "KI270899.1", "chr7_KI270899v1_alt", },
    { "KI270900.1", "chr8_KI270900v1_alt", },
    { "KI270901.1", "chr8_KI270901v1_alt", },
    { "KI270902.1", "chr11_KI270902v1_alt", },
    { "KI270903.1", "chr11_KI270903v1_alt", },
    { "KI270904.1", "chr12_KI270904v1_alt", },
    { "KI270905.1", "chr15_KI270905v1_alt", },
    { "KI270906.1", "chr15_KI270906v1_alt", },
    { "JH159148.1", "chr17_JH159148v1_alt", },
    { "KI270907.1", "chr17_KI270907v1_alt", },
    { "KI270908.1", "chr17_KI270908v1_alt", },
    { "KI270909.1", "chr17_KI270909v1_alt", },
    { "KI270910.1", "chr17_KI270910v1_alt", },
    { "KI270911.1", "chr18_KI270911v1_alt", },
    { "KI270912.1", "chr18_KI270912v1_alt", },
    { "GL949747.2", "chr19_GL949747v2_alt", },
    { "KB663609.1", "chr22_KB663609v1_alt", },
    { "KI270913.1", "chrX_KI270913v1_alt", },
    { "KI270924.1", "chr3_KI270924v1_alt", },
    { "KI270925.1", "chr4_KI270925v1_alt", },
    { "GL000252.2", "chr6_GL000252v2_alt", },
    { "KI270926.1", "chr8_KI270926v1_alt", },
    { "KI270927.1", "chr11_KI270927v1_alt", },
    { "GL949748.2", "chr19_GL949748v2_alt", },
    { "KI270928.1", "chr22_KI270928v1_alt", },
    { "KI270934.1", "chr3_KI270934v1_alt", },
    { "GL000253.2", "chr6_GL000253v2_alt", },
    { "GL949749.2", "chr19_GL949749v2_alt", },
    { "KI270935.1", "chr3_KI270935v1_alt", },
    { "GL000254.2", "chr6_GL000254v2_alt", },
    { "GL949750.2", "chr19_GL949750v2_alt", },
    { "KI270936.1", "chr3_KI270936v1_alt", },
    { "GL000255.2", "chr6_GL000255v2_alt", },
    { "GL949751.2", "chr19_GL949751v2_alt", },
    { "KI270937.1", "chr3_KI270937v1_alt", },
    { "GL000256.2", "chr6_GL000256v2_alt", },
    { "GL949752.1", "chr19_GL949752v1_alt", },
    { "KI270758.1", "chr6_KI270758v1_alt", },
    { "GL949753.2", "chr19_GL949753v2_alt", },
    { "KI270938.1", "chr19_KI270938v1_alt", },
    { "KI270882.1", "chr19_KI270882v1_alt", },
    { "KI270883.1", "chr19_KI270883v1_alt", },
    { "KI270884.1", "chr19_KI270884v1_alt", },
    { "KI270885.1", "chr19_KI270885v1_alt", },
    { "KI270886.1", "chr19_KI270886v1_alt", },
    { "KI270887.1", "chr19_KI270887v1_alt", },
    { "KI270888.1", "chr19_KI270888v1_alt", },
    { "KI270889.1", "chr19_KI270889v1_alt", },
    { "KI270890.1", "chr19_KI270890v1_alt", },
    { "KI270891.1", "chr19_KI270891v1_alt", },
    { "KI270914.1", "chr19_KI270914v1_alt", },
    { "KI270915.1", "chr19_KI270915v1_alt", },
    { "KI270916.1", "chr19_KI270916v1_alt", },
    { "KI270917.1", "chr19_KI270917v1_alt", },
    { "KI270918.1", "chr19_KI270918v1_alt", },
    { "KI270919.1", "chr19_KI270919v1_alt", },
    { "KI270920.1", "chr19_KI270920v1_alt", },
    { "KI270921.1", "chr19_KI270921v1_alt", },
    { "KI270922.1", "chr19_KI270922v1_alt", },
    { "KI270923.1", "chr19_KI270923v1_alt", },
    { "KI270929.1", "chr19_KI270929v1_alt", },
    { "KI270930.1", "chr19_KI270930v1_alt", },
    { "KI270931.1", "chr19_KI270931v1_alt", },
    { "KI270932.1", "chr19_KI270932v1_alt", },
    { "KI270933.1", "chr19_KI270933v1_alt", },
    { "GL000209.2", "chr19_GL000209v2_alt", },
};
    
int parse_args(int argc, char **argv)
{
    if (argc == 1)
	error("Failed to parse arguments. Use -h for more information.");
    int i;
    for ( i = 0; i < argc; ++i )
	kputs(argv[i], &args.command);    
    --argc, ++argv;
    const char *input_fname = 0;
    const char *output_fname = 0;
    const char *out_type_string = 0;

    for ( i = 0; i < argc; ) {
	const char *a = argv[i++];
	if ( strcmp(a, "-h") == 0 )
	    return usage();
	const char **arg_var = 0;
	if ( strcmp(a, "-list") == 0 )
	    arg_var = &list_fname;
	else if ( strcmp(a, "-O") == 0 && out_type_string == 0 )
	    arg_var = &out_type_string;
	else if ( strcmp(a, "-o") == 0 && output_fname == 0 )
	    arg_var = &output_fname;

	if ( arg_var != 0 ) {
	    if ( i == argc )
		error("Missing arg after %s.", a);
	    *arg_var = argv[i++];
	    continue;
	}

        if ( strcmp(a, "-ncbi") == 0 ) {
            if ( convert_ucsc_name_flag == 1 )
                error( "-ncbi and -ucsc conflict.");
            convert_ncbi_name_flag = 1;
            continue;            
        }
        else if ( strcmp(a, "-ucsc") == 0 ) {
            if ( convert_ncbi_name_flag == 1 )
                error( "-ncbi and -ucsc conflict.");                        
            convert_ucsc_name_flag = 1;
            continue;
        }
        
	if ( input_fname == 0 ) {
	    input_fname = a;
	    continue;
	}
	error("Unknown arg %s.", a);
    }

    if ( input_fname == 0 && (!isatty(fileno(stdin))) )
	input_fname = "-";
    if ( input_fname == 0 )
	error("No input file!");
    if ( list_fname == 0 && convert_ucsc_name_flag == 0 && convert_ncbi_name_flag == 0)
	error("No list file or -ucsc or -ncbi specified.");

    args.fp_input = hts_open(input_fname, "r");
    if (args.fp_input == NULL)
	error("Failed to open %s : %s.", input_fname, strerror(errno));

    htsFormat type = *hts_get_format(args.fp_input);
    if (type.format  != vcf && type.format != bcf)
	error("Unsupported input format! %s", input_fname);

    int out_type = FT_VCF;
    if (out_type_string != 0) {
	switch (out_type_string[0]) {
	    case 'b':
		out_type = FT_BCF_GZ; break;
	    case 'u':
		out_type = FT_BCF; break;
	    case 'z':
		out_type = FT_VCF_GZ; break;
	    case 'v':
		out_type = FT_VCF; break;
	    default :
		error("The output type %s not recognised\n", out_type_string);
	};
    }
    args.hdr_in = bcf_hdr_read(args.fp_input);
    args.fp_output = output_fname == 0 ? hts_open("-", hts_bcf_wmode(out_type)) : hts_open(output_fname, hts_bcf_wmode(out_type));
    args.hdr_out = bcf_hdr_dup(args.hdr_in);
    bcf_hdr_append(args.hdr_out, args.command.s);
    bcf_hdr_sync(args.hdr_out);
    return 0;
}

int tags_convert()
{
    if (list_fname == NULL)
        goto ucsc_ncbi_flag;
    
    int n, i;
    char **names = hts_readlist(list_fname, 1, &n);
    if (n == 0)
	return 1;

    kstring_t string = KSTRING_INIT;
    int nfields = 0;
    
    for ( i = 0; i < n; ++i ) {
	// skip comment lines
	if (names[i] == NULL || names[i][0] == '#')
	    continue;
	// init string
	string.l = 0;
	kputs(names[i], &string);
	// split string by tab
	int *splits = ksplit(&string, '\t', &nfields);
	if (nfields != 2) {
	    warnings("Error format; only accept two columns per line. %s", names[i]);
	    continue;
	}
	char *ss = string.s + splits[0];
	char *se = string.s + splits[1];
	if ( strlen(ss) > 4 && strncmp(ss, "CTG/", 4) == 0 ) {
	    ss += 4;
	    if ( strlen(se) > 4 && strncmp(se, "CTG/", 4) == 0 )
		ss += 4;
	    int rid = bcf_hdr_name2id(args.hdr_out, ss);
	    bcf_hrec_t *rec = bcf_hdr_get_hrec(args.hdr_out, BCF_HL_CTG, "ID", ss, NULL);
	    if ( !rec )
		continue;
	    int j = bcf_hrec_find_key(rec, "ID");
	    assert(j >= 0);
	    free(rec->vals[j]);
	    rec->vals[j] = strdup(se);
	    args.hdr_out->id[BCF_DT_CTG][rid].key = rec->vals[j];
	    free(splits);
	    continue;
	}
	    
	if ( strlen(ss) > 5 && strncmp(ss, "INFO/", 5) == 0 )
	    ss += 5;
	if ( strlen(se) > 5 && strncmp(se, "INFO/", 5) == 0 )
	    se += 5;
	int rid = bcf_hdr_name2id(args.hdr_out, ss);
	bcf_hrec_t *rec = bcf_hdr_get_hrec(args.hdr_out, BCF_HL_INFO, "ID", ss, NULL);
	if ( !rec )
	    continue;
	int j = bcf_hrec_find_key(rec, "ID");
	assert(j >= 0);
	free(rec->vals[j]);
	rec->vals[j] = strdup(se);
	args.hdr_out->id[BCF_DT_ID][rid].key = rec->vals[j];	
	free(splits);
    }
    for ( i = 0; i < n; ++i )
	free(names[i]);
    free(names);
    if ( string.m )
	free(string.s);

  ucsc_ncbi_flag:
    // convert NCBI/UCSC chromosome names
    if ( convert_ncbi_name_flag == 1 || convert_ucsc_name_flag == 1) {
        int i, n;
        n = sizeof(chr_pairs) / sizeof (chr_pairs[0]);
        for ( i = 0; i < n; ++i ) {
            char *name, *repl;
            if ( convert_ncbi_name_flag == 1 ) {
                name = chr_pairs[i].ucsc;
                repl = chr_pairs[i].ncbi;
            }
            else {
                name = chr_pairs[i].ncbi;
                repl = chr_pairs[i].ucsc;
            }
            
            //int rid = bcf_hdr_id2int(args.hdr_out, BCF_DT_ID, name);
            int rid = bcf_hdr_name2id(args.hdr_out, name);
            if ( rid == -1 )
                continue;
            bcf_hrec_t *rec = bcf_hdr_get_hrec(args.hdr_out, BCF_HL_CTG, "ID", name, NULL);
            if ( !rec )
                continue;
            int j = bcf_hrec_find_key(rec, "ID");
            assert(j >= 0);
            free(rec->vals[j]);
            rec->vals[j] = strdup(repl);
            args.hdr_out->id[BCF_DT_CTG][rid].key = rec->vals[j];            
        }
    }
    
    return 0;
}
int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) == 1)
	return 1;
    tags_convert();
    bcf_hdr_write(args.fp_output, args.hdr_out);
    bcf1_t *line = bcf_init();
    while ( bcf_read(args.fp_input, args.hdr_in, line) == 0 ) {
	bcf_write1(args.fp_output, args.hdr_out, line);
    }
    bcf_destroy(line);
    args_destroy();
    return 0;
}

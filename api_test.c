#include <my_global.h>
#include <mysql.h>

int main(int argc, char **argv)
{
    MYSQL *connect;
    connect = mysql_init(NULL);
    if ( connect == NULL ) { 
        fprintf(stderr, "%s\n", mysql_error(connect));
        return 1;
    }
    my_bool secure_auth = 0;
    mysql_options(connect, MYSQL_SECURE_AUTH, &secure_auth);
    if ( mysql_real_connect(connect, "10.1.10.58", "xuhuixin", "ss@Ed2Th", "hgmd_pro", 0, NULL, 0) == NULL ) {
        fprintf(stderr, "%s\n", mysql_error(connect));
        mysql_close(connect);
        return 1;
    }

    if ( mysql_query(connect, "SELECT hgmd_hg38_vcf.chrom, hgmd_hg38_vcf.pos, hgmd_hg38_vcf.ref, hgmd_hg38_vcf.alt, allmut.* FROM allmut, hgmd_hg38_vcf where allmut.acc_num = hgmd_hg38_vcf.id") ) {
        fprintf(stderr, "%s\n", mysql_error(connect));
        return 1;
    }

    MYSQL_RES *result = mysql_store_result(connect);

    if ( result == NULL ) {
        // finish_with_error(connect);
        fprintf(stderr, "%s\n", mysql_error(connect));
        return 1;
    }

    int num_fields = mysql_num_fields(result);

    MYSQL_ROW row;
    int i;
    while ((row = mysql_fetch_row(result))) {
        for ( i=0; i<num_fields; ++i ) {
            printf("%s\t", row[i] ? row[i] : "NULL\t");
        }
        printf("\n");
    }

    mysql_free_result(result);
    mysql_close(connect);
    return 0;
}

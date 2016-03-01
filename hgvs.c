#include "hgvs.h"
#include <htslib/faidx.h>

struct refgene * retrieve_refgene_from_local(const char *fname, int retrieve_rule)
{
}

void release_refgene_list(struct refgene *root)
{
    if (root == NULL) return;
    assert(root->is_root == 1);
    struct refgene ** pp = &root;
    struct refgene *head = root;
    while (head)
    {
	if (head->gene_name != NULL) {
	    safe_free(head->gene_name);
	}
	if (head->trans_name != NULL) {
	    safe_free(head->trans_name);

	if (head->exon_count != 0) {
	    safe_free(head->exon_starts);
	    safe_free(head->exon_stops);
	}
	*p = head->next;
	safe_free(head);
	head = *p;
    }
}

#include "hgvs.h"
#include <htslib/faidx.h>

void extract_refgene(struct refgene_entry *entry, int type)
{
    
}

struct refgene_entry *refgene_entry_praser (struct refgene_entry *entry, char *string, int rule)
{
    assert(entry->buffer == NULL);
    entry->buffer = string;
    extract_refgene(entry, rule);
    return entry;
}

/*
  input: 
    start: begin position
    end:  end position
    
  return
    *buffers: the refgene entries buffer
    *nbuffers: the entry number
 */
struct refgene_entry * retrieve_refgene_from_local(const char *fname, int retrieve_rule, int start, int end, int *nbuffer)
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

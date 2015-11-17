#include "anno_setter.h"
#include "vcmp.h"

int anno_setter_id(anno_setters_t *handler, bcf1_t *line, anno_col_t *col, void *data)
{
    anno_line_t *tab = (anno_line_t*)data;
    if ( tab->cols[col->icol] && tab->cols[col->icol][0]=='.' && !tab->cols[col->icol][1])
	return 0; // donot replace with '.'
    if ( col->replace != REPLACE_MISSING)
	return bcf_update_id(handler->hdr_out, line, tab->cols[col->icol]);

    if ( !line->d.id || (line->d.id[0]=='.' && !line->d.id[1]) )
	return bcf_update_id()
}

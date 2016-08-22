/* I/O of multiple sequence alignments in PHYLIP format
 */
#ifndef eslMSAFILE_PHYLIP_INCLUDED
#define eslMSAFILE_PHYLIP_INCLUDED

#include "esl_msa.h"
#include "esl_msafile.h"

extern int esl_msafile_phylip_SetInmap     (ESLX_MSAFILE *afp);
extern int esl_msafile_phylip_GuessAlphabet(ESLX_MSAFILE *afp, int *ret_type);
extern int esl_msafile_phylip_Read         (ESLX_MSAFILE *afp, ESL_MSA **ret_msa);
extern int esl_msafile_phylip_Write        (FILE *fp, const ESL_MSA *msa, int format, ESLX_MSAFILE_FMTDATA *opt_fmtd);

extern int esl_msafile_phylip_CheckFileFormat(ESL_BUFFER *bf, int *ret_format, int *ret_namewidth);

#endif /* eslMSAFILE_PHYLIP_INCLUDED */
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version i1.1.1; July 2014
 * Copyright (C) 2014 HHMI Janelia Farm Research Campus
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *
 * SVN $Id: esl_msafile_phylip.h 708 2011-07-20 12:49:10Z eddys $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/easel/branches/infernal/1.1/esl_msafile_phylip.h $
 *****************************************************************/

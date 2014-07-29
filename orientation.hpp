#ifndef _ORIENTATION_H_
#define _ORIENTATION_H_

#include "sam.h"

enum Orientation { FR, RF, FF, RR, UNMAPPED };

bool pair_is_mapped(bam1_t* read1, bam1_t* read2);

// Functions to determine orientation of a pair
bool pair_is_FR(bam1_t* read1, bam1_t* read2);
bool pair_is_RF(bam1_t* read1, bam1_t* read2);
bool pair_is_RR(bam1_t* read1, bam1_t* read2);
bool pair_is_FF(bam1_t* read1, bam1_t* read2);
Orientation get_orientation(bam1_t *read1, bam1_t *read2);

#endif /* _ORIENTATION_H_ */

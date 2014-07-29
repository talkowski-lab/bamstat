#ifndef _BAMSTAT_H_
#define _BAMSTAT_H_

#include <cstdlib>
#include <cstdarg>
#include <iostream>
#include <string>
#include <unistd.h>
#include <getopt.h>
#include "sam.h"
#include "alignment_stats.hpp"
#include "running_mean.hpp"
#include "running_median.hpp"
#include "variant_type.hpp"
#include "samstring.hpp"
#include "orientation.hpp"

#define MAX_INSERT 5000

int usage();
void sample_statistics(samfile_t* samfile, int* FR_interval, int* RF_interval,
    int threshold, int max_insert);
VariantType call_variant(bam1_t* read1, bam1_t* read2, AlignmentStats* stats,
    int max_insert);

bool is_proper_sample(bam1_t* read1, bam1_t* read2, int max_insert);
bool is_proper_RF(bam1_t* read1, bam1_t* read2, int max_insert);

#endif /* _BAMSTAT_H_ */

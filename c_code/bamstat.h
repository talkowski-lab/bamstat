#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>

#include "sam.h"
#include "malloc.h"
#include "alignment_stats.h"
#include "running_stat.h"

int sample_statistics(samfile_t*);
int is_proper_sample(bam1_t*, bam1_t*);
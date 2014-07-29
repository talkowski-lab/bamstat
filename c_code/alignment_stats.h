#include <stdio.h>

/* Definition of AlignmentStats structure
   Used to track statistics of reads in a bam file */
typedef struct _AlignmentStats {

    int mapped;
    int unmapped;
    int total;
    int pos_mapq;

    int proper_RR_count;
    int proper_RF_count;
    int proper_FR_count;
    int proper_FF_count;

    int notproper_RR_count;
    int notproper_RF_count;
    int notproper_FR_count;
    int notproper_FF_count;

    int notproper_neither_mapped;
    int notproper_one_mapped;
    int notproper_both_mapped;

    int proper_same_chr;
    int proper_diff_chr;

    int notproper_same_chr;
    int notproper_diff_chr;

    int putative_RR_deletions;
    int putative_RF_deletions;
    int putative_FR_deletions;
    int putative_FF_deletions;

    int putative_RR_insertions;
    int putative_RF_insertions;
    int putative_FR_insertions;
    int putative_FF_insertions;

} __AlignmentStats;

typedef struct _AlignmentStats AlignmentStats;

/* Function prototypes */
void init_stats(AlignmentStats*);
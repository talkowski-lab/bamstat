#include "alignment_stats.h"

/* Initialize all counters to zero */
void init_stats(AlignmentStats *stats) {

    stats->mapped = 0;
    stats->unmapped = 0;
    stats->total = 0;
    stats->pos_mapq = 0;

    stats->proper_RR_count = 0;
    stats->proper_RF_count = 0;
    stats->proper_FR_count = 0;
    stats->proper_FF_count = 0;

    stats->notproper_RR_count = 0;
    stats->notproper_RF_count = 0;
    stats->notproper_FR_count = 0;
    stats->notproper_FF_count = 0;

    stats->notproper_neither_mapped = 0;
    stats->notproper_one_mapped = 0;
    stats->notproper_both_mapped = 0;

    stats->proper_same_chr = 0;
    stats->proper_diff_chr = 0;

    stats->notproper_same_chr = 0;
    stats->notproper_diff_chr = 0;

    stats->putative_RR_deletions = 0;
    stats->putative_RF_deletions = 0;
    stats->putative_FR_deletions = 0;
    stats->putative_FF_deletions = 0;

    stats->putative_RR_insertions = 0;
    stats->putative_RF_insertions = 0;
    stats->putative_FR_insertions = 0;
    stats->putative_FF_insertions = 0;
}

/* Output final results to file */
void write_stats(AlignmentStats *stats, FILE *file) {
    /* FILE *file = fopen("stats.file", 'w'); */

    fprintf(file, "Mapped Reads\t%d\n", stats->mapped);
    fprintf(file, "Total Reads\t%d\n", stats->total);
    fprintf(file, "\n");

    fprintf(file, "Counts of (Mapped) Pairs\n");
    fprintf(file, "----\n");
    fprintf(file, "Orientation\tProper Pairs\tImproper Pairs\n");
    fprintf(file, "RR\t\t%d\t\t%d\n",
        stats->proper_RR_count, stats->notproper_RR_count);
    fprintf(file, "RF\t\t%d\t\t%d\n",
        stats->proper_RF_count, stats->notproper_RF_count);
    fprintf(file, "FR\t\t%d\t\t%d\n",
        stats->proper_FR_count, stats->notproper_FR_count);
    fprintf(file, "FF\t\t%d\t\t%d\n",
        stats->proper_FF_count, stats->notproper_FF_count);
    fprintf(file, "\n");

    fprintf(file, "Properly Paired / Same Chr\t%d\n",
        stats->proper_same_chr);
    fprintf(file, "Properly Paired / Diff Chr\t%d\n",
        stats->proper_diff_chr);
    fprintf(file, "\n");

    fprintf(file, "Not Paired / Same Chr\t\t%d\n",
        stats->notproper_same_chr);
    fprintf(file, "Not Paired / Diff Chr\t\t%d\n",
        stats->notproper_diff_chr);
    fprintf(file, "\n");

    fprintf(file, "Not Paired / Neither Mapped\t%d\n",
        stats->notproper_neither_mapped);
    fprintf(file, "Not Paired / One End Mapped\t%d\n",
        stats->notproper_one_mapped);
    fprintf(file, "Not Paired / Both Ends Mapped\t%d\n",
        stats->notproper_both_mapped);
    fprintf(file, "\n");

    fprintf(file, "Reads Mapped\t\t\t%d\n", stats->mapped);
    fprintf(file, "Reads Unmapped\t\t\t%d\n", stats->unmapped);
    fprintf(file, "Reads Total\t\t\t%d\n", stats->total);
    fprintf(file, "Reads with Positive MAPQ\t%d\n", stats->pos_mapq);
    fprintf(file, "\n");

}
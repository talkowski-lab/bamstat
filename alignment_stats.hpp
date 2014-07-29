#ifndef _ALIGNMENT_STATS_H_
#define _ALIGNMENT_STATS_H_

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include "variant_type.hpp"
#include "orientation.hpp"
#include "running_mean.hpp"
#include "running_median.hpp"

class AlignmentStats {
    int mapped_reads;
    int unmapped_reads;
    int total_reads;

    int proper_pairs;
    int mapped_pairs;

    int pos_mapq;

    int notproper_RR_count;
    int notproper_RF_count;
    int notproper_FR_count;
    int notproper_FF_count;

    int notproper_neither_mapped;
    int notproper_one_mapped;
    int notproper_both_mapped;

    int deletions;
    int insertions;
    int inversions;
    int translocations;
    int short_RF;

    RunningMean FR_insert_mean;
    RunningMean RF_insert_mean;
    RunningMedian *FR_insert_median;
    RunningMedian *RF_insert_median;

public:
    AlignmentStats();
    ~AlignmentStats();
    void print(std::string filename);
    void update(VariantType variant, int isize, Orientation o);

};

#endif /* _ALIGNMENT_STATS_H_ */

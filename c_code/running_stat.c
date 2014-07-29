#include "running_stat.h"

void reset_running(RunningStat *stat) {
    stat->num_data = 0;
    stat->old_mean = 0.0;
    stat->mean = 0.0;
    stat->old_s = 0.0;
    stat->s = 0.0;
    stat->variance = 0.0;
    stat->std_dev = 0.0;
    stat->min = 0;
    stat->max = 0;
}

void push_stat(RunningStat *stat, int value) {
    stat->num_data++;

    if (stat->num_data == 1) {
        stat->old_mean = value;
        stat->mean = value;
        stat->min = value;
        stat->max = value;
    /*    } else if (stat->num_data > 1) { */
    } else {
        stat->mean = stat->old_mean +
            ((value - stat->old_mean) / stat->num_data);
        stat->s = stat->old_s +
            ((value - stat->old_mean) * (value - stat->mean));

        stat->variance = stat->s / (stat->num_data - 1);
        stat->std_dev = sqrt(stat->variance);

        stat->old_mean = stat->mean;
        stat->old_s = stat->s;

        if (value > stat->max) {
            stat->max = value;
        }
        if (value < stat->min) {
            stat->min = value;
        }
    }
    /*
    } else {
        fprintf(stderr, "Error: ")
    }
    */
}

void print_running(RunningStat *stat) {
    printf("Mean: %f\n", stat->mean);
    printf("Standard Deviation: %f\n", stat->std_dev);
}
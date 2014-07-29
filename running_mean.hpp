#ifndef _RUNNING_MEAN_H_
#define _RUNNING_MEAN_H_

#include <math.h>

class RunningMean {
    int num_data;
    float old_mean;
    float old_s;
    float s;
    float variance;
    float min;
    float max;

public:
    float mean;
    float std_dev;

    RunningMean();
    void push(int value);
    void clear();
    void write_stats();
    // void RunningMean::push(int);
    // void RunningMean::write_stats();
};

#endif /* _RUNNING_MEAN_H_ */

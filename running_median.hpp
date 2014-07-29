#ifndef _RUNNING_MEDIAN_H_
#define _RUNNING_MEDIAN_H_

#include <cstdlib>

class RunningMedian {
    int min;
    int max;
    int num_counts;
    int total_values;
    int *counts;

public:
    RunningMedian(int min, int max);
    ~RunningMedian();
    void push(int value);
    int get_median();
    int mad();
};

#endif /* _RUNNING_MEDIAN_H_ */

#include <math.h>
#include <stdio.h>

typedef struct _RunningStat {

    int num_data;
    float old_mean;
    float mean;
    float old_s;
    float s;
    float variance;
    float std_dev;
    float min;
    float max;

} __RunningStat;

typedef struct _RunningStat RunningStat;

void reset_running(RunningStat*);
void push_stat(RunningStat*, int);
void print_running(RunningStat*);
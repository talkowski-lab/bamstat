#include "running_mean.hpp"

RunningMean::RunningMean() {
    this->num_data = 0;
    this->old_mean = 0.0;
    this->old_s = 0.0;
    this->s = 0.0;
    this->variance = 0.0;
    this->min = 0;
    this->max = 0;
    this->mean = 0.0;
    this->std_dev = 0.0;
}

void RunningMean::push(int value) {
    this->num_data++;

    if (this->num_data == 1) {
        this->old_mean = value;
        this->mean = value;
        this->min = value;
        this->max = value;
    } else {
        this->mean = this->old_mean +
            ((value - this->old_mean) / this->num_data);
        this->s = this->old_s +
            ((value - this->old_mean) * (value - this->mean));

        this->variance = this->s / (this->num_data - 1);
        this->std_dev = sqrt(this->variance);

        this->old_mean = this->mean;
        this->old_s = this->s;

        if (value > this->max){
            this->max = value;
        }
        if (value < this->min) {
            this->min = value;
        }
    }
}

void RunningMean::clear() {
    this->num_data = 0;
    this->old_mean = 0.0;
    this->mean = 0.0;
    this->old_s = 0.0;
    this->s = 0.0;
    this->variance = 0.0;
    this->std_dev = 0.0;
    this->min = 0;
    this->max = 0;
}
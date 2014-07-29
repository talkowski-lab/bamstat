#include "running_median.hpp"

RunningMedian::RunningMedian(int min, int max) {
    this->min = min;
    this->max = max;
    this->num_counts = max - min + 3;
    this->total_values = 0;
    this->counts = new int[this->num_counts];

    for (int i = 0; i < this->num_counts; i++) {
        this->counts[i] = 0;
    }
}

RunningMedian::~RunningMedian() {
    delete counts;
}

void RunningMedian::push(int value) {
    this->total_values++;

    if (value < this->min) {
        this->counts[0]++;
    } else if (value > this->max) {
        this->counts[this->num_counts - 1]++;
    } else {
        this->counts[value - this->min + 1]++;
    }
}

int RunningMedian::get_median() {
    int current_count = 0;
    int half_count = (this->total_values / 2) + 1;
    int last_bin = 0;
    int i = 0;
    bool odd_count = ((current_count % 2) != 0);

    if (this->total_values == 0) {
        return 0;
    }

    while (current_count < half_count) {
        current_count += this->counts[i];
        if (this->counts[i] > 0) {
            last_bin = i;
        }
        i++;
    }

    if (current_count == half_count){
        return i - 2 + this->min;
    } else {
        return last_bin - 2 + this->min;
    }
}

// Calculate median absolute deviation
int RunningMedian::mad() {
    RunningMedian abs_deviations(0, this->max);
    int median = this->get_median();
    int value;

    for (int i = 0; i < num_counts; i++) {
        value = i + min - 1;
        for (int j = 0; j < counts[i]; j++) {
            abs_deviations.push(abs(median - value));
        }
    }

    return abs_deviations.get_median();
}
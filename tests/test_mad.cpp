#include <cstdlib>
#include <cstdarg>
#include <iostream>

#include "running_median.hpp"

void print_array(int array[], int length) {
    std::cout << "{";
    int i;
    for (i = 0; i < length - 1; i++) {
        std::cout << array[i] << ",  ";
    }
    std::cout << array[i] << "}" << std::endl;
}

int main(int argc, char** argv) {
    RunningMedian *run = new RunningMedian(10, 50);
    RunningMedian *deviations = new RunningMedian(0,50);
    int med;

    int data[40] = {48, 35, 43, 15, 10, 19, 22, 36, 10, 42, 33, 31, 10, 16, 13, 42, 33, 24, 16, 12, 42, 26, 14, 16, 28, 47, 44, 23, 28, 26, 43, 44, 27, 27, 37, 14, 10, 12, 46, 46};
    for (int i = 0; i < 40; i++) {
        run->push(data[i]);
    }

    std::cout << run->mad() << std::endl;

    // std::cout << "Median: " << run->get_median() << std::endl;
    // std::cout << "MAD:    " << run->mad() << std::endl;
    // std::cout << "Median: " << run->get_median() << std::endl;
    // std::cout << "MAD:    " << run->mad() << std::endl;
    // std::cout << "Median: " << run->get_median() << std::endl;
    // std::cout << "MAD:    " << run->mad() << std::endl;
    // print_array(data, 40);
    // print_array(run->counts, 43);
    // std::cout << "     ";
    // for (int i = 10; i <= 50; i++) {
    //     std::cout << i << "  ";
    // }
    // std::cout << std::endl;
    // std::cout << "total values: " << run->total_values << std::endl;
    // std::cout << "Median 2: " << median(run) << std::endl;
    // std::cout << "total values: " << run->total_values << std::endl;
    // std::cout << "Median 2: " << median(run) << std::endl;
    // std::cout << "Median: " << run->get_median() << std::endl;
    // std::cout << "Median 2: " << median(run) << std::endl;

    // print_array(run->counts, 43);
    // std::cout << "     ";
    // for (int i = 10; i <= 50; i++) {
    //     std::cout << i << "  ";
    // }
    // std::cout << std::endl;
    // std::cout << "Median: " << run->get_median() << std::endl;
    // std::cout << "Median 2: " << median(run) << std::endl;

    // med = run->get_median();
    // // std::cout << med << std::endl;

    // for (int i = 0; i < 40; i++) {
    //     deviations->push(abs(med - data[i]));
    // }

    // // int med = median.get_median();
    // std::cout << "Median: " << median->get_median() << std::endl;
    // std::cout << "MAD:    " << deviations->get_median() << std::endl;
}
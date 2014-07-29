#ifndef _READ_PAIR_H_
#define _READ_PAIR_H_

#include <iostream>
#include <string>
#include <sstream>
#include "sam.h"

class ReadPair {
    // bam1_t *readA;
    // bam1_t *readB;

    // std::string readID;

    int tidA;
    std::string chrA;
    int posA;
    int scoreA;
    int lenA;
    int flagA;
    std::string seqA;
    std::string qualA;

    int tidB;
    std::string chrB;
    int posB;
    int scoreB;
    int lenB;
    int flagB;
    std::string seqB;
    std::string qualB;

public:
    std::string readID;
    int isize;

    ReadPair(bam1_t* readA, bam1_t* readB, samfile_t* samfile);
    ~ReadPair();

    void set_pair(bam1_t *read1, bam1_t *read2, samfile_t *samfile);
    int get_isize();
    bool is_FF();
    bool is_RR();
    bool is_FR();
    bool is_RF();
    bool is_proper(int max_insert);
    void print();
};

std::string get_seq(bam1_t* read);
std::string get_qual(bam1_t* read);

#endif /* _READ_PAIR_H_ */

#include "samstring.hpp"

// Returns query sequence as string
std::string get_seq(bam1_t *read){
    uint8_t *s = bam1_seq(read);
    std::stringstream seq;

    for (int i = 0; i < read->core.l_qseq; i++) {
        seq << bam_nt16_rev_table[bam1_seqi(s, i)];
    }

    return seq.str();
}

// Returns string of quality scores for each base in query sequence
std::string get_qual(bam1_t *read){
    uint8_t *q = bam1_qual(read);
    std::stringstream qual;

    if (q[0] == 0xff) {
        qual << "*";
    } else for (int i = 0; i < read->core.l_qseq; i++) {
        qual << (char) (q[i] + 33);
    }

    return qual.str();
}

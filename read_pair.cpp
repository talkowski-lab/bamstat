#include "read_pair.hpp"

// Construct read pair object from two aligned reads
ReadPair::ReadPair(bam1_t *read1, bam1_t *read2, samfile_t *samfile) {
    bam1_t *temp; // = bam_init1();

    // Order reads appropriately
    // Assign upstream read to read1 if both reads map to same chromosome
    // Otherwise assign read with lower chromosome number to read1
    if (((read1->core.tid == read2->core.tid) &&
         (read1->core.pos > read2->core.pos)) ||
        (read1->core.tid > read2->core.tid)) {
        temp = read1;
        read1 = read2;
        read2 = temp;
    }

    // Parse all necessary data from read pairs
    this->readID = bam1_qname(read1);
    this->isize = abs(read1->core.isize);

    this->tidA = read1->core.tid;
    if (read1->core.tid != -1) {
        this->chrA = samfile->header->target_name[read1->core.tid];
    } else {
        this->chrA = "chrM";
    }
    this->posA = read1->core.pos;
    this->scoreA = read1->core.qual;
    this->lenA = read1->core.l_qseq;
    this->flagA = read1->core.flag;
    this->seqA = get_seq(read1);
    this->qualA = get_qual(read1);

    this->tidB = read2->core.tid;
    if (read2->core.tid != -1) {
        this->chrB = samfile->header->target_name[read2->core.tid];
    } else {
        this->chrB = "chrM";
    }
    this->posB = read2->core.pos;
    this->scoreB = read2->core.qual;
    this->lenB = read2->core.l_qseq;
    this->flagB = read2->core.flag;
    this->seqB = get_seq(read2);
    this->qualB = get_qual(read2);

    // bam_destroy1(temp);
}

ReadPair::~ReadPair() {
    //delete readID;
}

void ReadPair::set_pair(bam1_t *read1, bam1_t *read2, samfile_t *samfile) {
    bam1_t *temp;

    // Order reads appropriately
    // Assign upstream read to read1 if both reads map to same chromosome
    // Otherwise assign read with lower chromosome number to read1
    if (((read1->core.tid == read2->core.tid) &&
         (read1->core.pos > read2->core.pos)) ||
        (read1->core.tid > read2->core.tid)) {
        temp = read1;
        read1 = read2;
        read2 = temp;
    }

    // Parse all necessary data from read pairs
    this->readID = bam1_qname(read1);
    this->isize = abs(read1->core.isize);

    this->tidA = read1->core.tid;
    this->chrA = samfile->header->target_name[read1->core.tid];
    this->posA = read1->core.pos;
    this->scoreA = read1->core.qual;
    this->lenA = read1->core.l_qseq;
    this->flagA = read1->core.flag;
    this->seqA = get_seq(read1);
    this->qualA = get_qual(read1);

    this->tidB = read2->core.tid;
    this->chrB = samfile->header->target_name[read2->core.tid];
    this->posB = read2->core.pos;
    this->scoreB = read2->core.qual;
    this->lenB = read2->core.l_qseq;
    this->flagB = read2->core.flag;
    this->seqB = get_seq(read2);
    this->qualB = get_qual(read2);
}

void ReadPair::print() {
    std::cout << this->readID << "\t";
    std::cout << this->isize << "\t";

    std::cout << this->tidA << "\t";
    std::cout << this->chrA << "\t";
    std::cout << this->posA << "\t";
    std::cout << this->scoreA << "\t";
    std::cout << this->lenA << "\t";
    std::cout << this->flagA << "\t";
    std::cout << this->seqA << "\t";
    std::cout << this->qualA << "\t";

    std::cout << this->tidB << "\t";
    std::cout << this->chrB << "\t";
    std::cout << this->posB << "\t";
    std::cout << this->scoreB << "\t";
    std::cout << this->lenB << "\t";
    std::cout << this->flagB << "\t";
    std::cout << this->seqB << "\t";
    std::cout << this->qualB << std::endl;
}

int ReadPair::get_isize() {
    return this->isize;
}

// Methods to determine pair orientation
bool ReadPair::is_FF() {
    return !(this->flagA & BAM_FREVERSE) && !(this->flagB & BAM_FREVERSE);
}

bool ReadPair::is_RR() {
    return (this->flagA & BAM_FREVERSE) && (this->flagB & BAM_FREVERSE);
}

bool ReadPair::is_FR() {
    return !(this->flagA & BAM_FREVERSE) && (this->flagB & BAM_FREVERSE);
}

bool ReadPair::is_RF() {
    return (this->flagA & BAM_FREVERSE) && !(this->flagB & BAM_FREVERSE);
}

// Determines whether pair represents a proper sample
// Proper sample is defined to be any pair fulfilling the following conditions:
//      Insert size greater than 0 and less than a given maximum
//      Both reads mapped
//      Both reads mapped to same chromosome
bool ReadPair::is_proper(int max_insert) {
    return ((this->isize > 0) && (this->isize < max_insert) &&
        !(this->flagA & BAM_FUNMAP) &&
        !(this->flagB & BAM_FUNMAP) &&
        (this->tidA == this->tidB));
}

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

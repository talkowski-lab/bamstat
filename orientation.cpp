#include "orientation.hpp"

// Return whether pair is of FF orientation
bool pair_is_FF(bam1_t *read1, bam1_t *read2) {
    return !(read1->core.flag & BAM_FREVERSE) &&
        !(read2->core.flag & BAM_FREVERSE);
}

// Return whether pair is of RR orientation
bool pair_is_RR(bam1_t *read1, bam1_t *read2) {
    return (read1->core.flag & BAM_FREVERSE) &&
        (read2->core.flag & BAM_FREVERSE);
}

// Return whether pair is of FR orientation
// If pair maps to different chromosomes, call read on lower
// chromosome as upstream read
bool pair_is_FR(bam1_t *read1, bam1_t *read2) {
    if (read1->core.tid == read2->core.tid) {
        if (read1->core.pos < read2->core.pos) {
            return !(read1->core.flag & BAM_FREVERSE) &&
                (read2->core.flag & BAM_FREVERSE);
        } else {
            return !(read2->core.flag & BAM_FREVERSE) &&
                (read1->core.flag & BAM_FREVERSE);
        }
    } else if (read1->core.tid < read2->core.tid) {
        return !(read1->core.flag & BAM_FREVERSE) &&
            (read2->core.flag & BAM_FREVERSE);
    } else {
        return !(read2->core.flag & BAM_FREVERSE) &&
            (read1->core.flag & BAM_FREVERSE);
    }
}


// Return whether pair is of RF orientation
// If pair maps to different chromosomes, call read on lower
// chromosome as upstream read
bool pair_is_RF(bam1_t *read1, bam1_t *read2) {
    if (read1->core.tid == read2->core.tid) {
        if (read1->core.pos < read2->core.pos) {
            return (read1->core.flag & BAM_FREVERSE) &&
                !(read2->core.flag & BAM_FREVERSE);
        } else {
            return (read2->core.flag & BAM_FREVERSE) &&
                !(read1->core.flag & BAM_FREVERSE);
        }
    } else if (read1->core.tid < read2->core.tid) {
        return (read1->core.flag & BAM_FREVERSE) &&
            !(read2->core.flag & BAM_FREVERSE);
    } else {
        return (read2->core.flag & BAM_FREVERSE) &&
            !(read1->core.flag & BAM_FREVERSE);
    }
}


// Return orientation of a pair
Orientation get_orientation(bam1_t *read1, bam1_t *read2) {
    if (pair_is_mapped(read1, read2)) {
        if (pair_is_FR(read1, read2))
            return FR;
        if (pair_is_RF(read1, read2))
            return RF;
        if (pair_is_FF(read1, read2))
            return FF;
        if (pair_is_RR(read1, read2))
            return RR;
    } else {
        return UNMAPPED;
    }
}

// Return whether both reads in pair are mapped
bool pair_is_mapped(bam1_t *read1, bam1_t *read2) {
    return !(read1->core.flag & BAM_FUNMAP) &&
        !(read2->core.flag & BAM_FUNMAP);
}
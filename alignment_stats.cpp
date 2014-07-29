#include "alignment_stats.hpp"

AlignmentStats::AlignmentStats() {
    mapped_reads = 0;
    unmapped_reads = 0;
    total_reads = 0;

    proper_pairs = 0;
    mapped_pairs = 0;

    pos_mapq = 0;

    notproper_RR_count = 0;
    notproper_RF_count = 0;
    notproper_FR_count = 0;
    notproper_FF_count = 0;

    notproper_neither_mapped = 0;
    notproper_one_mapped = 0;
    notproper_both_mapped = 0;

    deletions = 0;
    insertions = 0;
    inversions = 0;
    translocations = 0;
    short_RF = 0;

    FR_insert_median = new RunningMedian(500, 20000);
    RF_insert_median = new RunningMedian(1, 3000);
}

AlignmentStats::~AlignmentStats() {
    delete FR_insert_median;
    delete RF_insert_median;
}

void AlignmentStats::update(VariantType variant, int isize, Orientation o) {

    total_reads += 2;

    switch (variant) {
        case NORMAL:
            proper_pairs++;
            mapped_pairs++;
            mapped_reads += 2;
            FR_insert_median->push(isize);
            FR_insert_mean.push(isize);
            break;

        case DELETION:
            mapped_reads += 2;
            mapped_pairs++;
            notproper_both_mapped++;
            deletions++;
            break;

        case INSERTION:
            mapped_reads += 2;
            mapped_pairs++;
            notproper_both_mapped++;
            insertions++;
            break;

        case INVERSION:
            mapped_reads += 2;
            mapped_pairs++;
            notproper_both_mapped++;
            inversions++;
            break;

        case TRANSLOCATION:
            mapped_reads += 2;
            mapped_pairs++;
            notproper_both_mapped++;
            translocations++;
            break;

        case SHORT:
            mapped_reads += 2;
            mapped_pairs++;
            notproper_both_mapped++;
            short_RF++;
            RF_insert_median->push(isize);
            RF_insert_mean.push(isize);
            break;

        case ONE_UNMAPPED:
            notproper_one_mapped++;
            mapped_reads++;
            unmapped_reads++;
            break;

        case BOTH_UNMAPPED:
            notproper_neither_mapped++;
            unmapped_reads += 2;
            break;

        default:
            break;
    }

    if (variant != NORMAL) {
        switch(o) {
            case FR:
                notproper_FR_count++;
                break;

            case RF:
                notproper_RF_count++;
                break;

            case FF:
                notproper_FF_count++;
                break;

            case RR:
                notproper_RR_count++;
                break;

            // Unmapped, skip
            default:
                break;
        }
    }
}

void AlignmentStats::print(std::string filename) {

    std::ofstream statsfile(filename.c_str());

    statsfile << "Total Reads\t\t" << total_reads << std::endl;
    statsfile << "Mapped Reads\t\t" << mapped_reads << std::endl;
    statsfile << "Unmapped Reads\t\t" << unmapped_reads << std::endl;

    statsfile << "Mapped Read Pairs\t" << mapped_pairs << std::endl;
    statsfile << "Proper pairs\t\t" << proper_pairs << std::endl;
    statsfile << std::endl;

    statsfile << "Counts of Mapped Improper Pairs" << std::endl;
    statsfile << "----" << std::endl;
    statsfile << "Orientation" << std::endl;
    statsfile << "RR\t\t" << notproper_RR_count << std::endl;
    statsfile << "RF\t\t" << notproper_RF_count << std::endl;
    statsfile << "FR\t\t" << notproper_FR_count << std::endl;
    statsfile << "FF\t\t" << notproper_FF_count << std::endl;
    statsfile << std::endl;

    statsfile << "Improper pairs, neither end mapped\t" <<
        notproper_neither_mapped << std::endl;
    statsfile << "Improper pairs, one end mapped\t\t" <<
        notproper_one_mapped << std::endl;
    statsfile << "Improper pairs, both ends mapped\t" <<
        notproper_both_mapped << std::endl;
    statsfile << std::endl;

    statsfile << "Potential structural variants" << std::endl;
    statsfile << "Putative deletion pairs\t\t\t" << deletions << std::endl;
    statsfile << "Putative insertion pairs\t\t" << insertions << std::endl;
    statsfile << "Putative inversion pairs\t\t" << inversions << std::endl;
    statsfile << "Putative translocation pairs\t\t" << translocations
        << std::endl;
    statsfile << std::endl;

    statsfile << "Actual FR mean insert size:\t\t" << FR_insert_mean.mean
        << std::endl;
    statsfile << "Actual FR insert standard deviation:\t"
        << FR_insert_mean.std_dev << std::endl;
    statsfile << "Actual FR median insert size:\t\t"
        << FR_insert_median->get_median() << std::endl;
    statsfile << "Actual FR median absolute deviation:\t"
        << FR_insert_median->mad() << "\n" << std::endl;
    statsfile << "Actual RF mean insert size:\t\t" << RF_insert_mean.mean
        << std::endl;
    statsfile << "Actual RF insert standard deviation:\t"
        << RF_insert_mean.std_dev << std::endl;
    statsfile << "Actual RF median insert size:\t\t"
        << RF_insert_median->get_median() << std::endl;
    statsfile << "Actual RF median absolute deviation:\t"
        << RF_insert_median->mad() << std::endl;

    statsfile.close();
}

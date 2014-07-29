#include "bamstat.hpp"

/* Print usage information */
int usage() {
    std::cout << std::endl;

    std::cout << "Program: bamstat (SAM read statistics and sorting)"
        << std::endl;
    std::cout << "Version: 0.2.0\n" << std::endl;

    std::cout << "Usage: bamstat -i <filename> -o <directory> [options]\n"
        << std::endl;

    std::cout << "Required arguments:" << std::endl;
    std::cout << "    -i --input-file <filename>\t\tRead in from <filename>"
        << std::endl;
    // std::cout << "    -o --output-dir <directory>\t\tStore output in "
    //     "<directory>" << std::endl;
    std::cout << std::endl;

    std::cout << "Optional arguments, with defaults in [ ]:" << std::endl;
    std::cout << "    -h --help\t\t\t\tShow this message" << std::endl;
    std::cout << "    --insert-deviations=<threshold>\t"
        "Define normal insert size to be within\n"
        "\t\t\t\t\t<threshold> absolute deviations of the median [5]" << std::endl;
    std::cout << "    --max-insert=<insert-size>\t\t"
        "Define maximum insert size for calculation\n"
        "\t\t\t\t\tof sample mean to be <insert-size> [5000]" << std::endl;
    std::cout << "    --bam-input\t\t\t\tRead from .bam format instead "
        "of .sam" << std::endl;
    std::cout << "    --save-unmapped\t\t\tOutput unmapped read pairs\n"
        << std::endl;

    return 0;
}

// Write pair of reads to plaintext file
// Assumes reads have already been sorted by orientation
void write_pair(std::ofstream &outfile, bam1_t *read1, bam1_t *read2,
    samfile_t *samfile) {

    // if (read1->core.tid != -1 && read2->core.tid != -1) {
    //     outfile << bam1_qname(read1) << "\t"
    //         << abs(read1->core.isize) << "\t"
    //         << samfile->header->target_name[read1->core.tid] << "\t"
    //         << samfile->header->target_name[read2->core.tid] << "\t"
    //         << read1->core.pos << "\t"
    //         << read2->core.pos << "\t"
    //         << read1->core.qual << "\t"
    //         << read2->core.qual << "\t"
    //         << read1->core.l_qseq << "\t"
    //         << read2->core.l_qseq << "\t"
    //         << read1->core.flag << "\t"
    //         << read2->core.flag << "\t"
    //         << get_seq(read1) << "\t"
    //         << get_seq(read2) << "\t"
    //         << get_qual(read1) << "\t"
    //         << get_qual(read2) << std::endl;
    // }

    // old formatting to preserve compatibility with readPairCluster
    // rpcluster will take above formatting
    if (read1->core.tid != -1 && read2->core.tid != -1) {
        outfile << bam1_qname(read1) << "\t"
            << samfile->header->target_name[read1->core.tid] << "\t"
            << read1->core.pos << "\t"
            << read1->core.flag << "\t"
            << samfile->header->target_name[read2->core.tid] << "\t"
            << read2->core.pos << "\t"
            << read2->core.flag << "\t"
            << read1->core.l_qseq << "\t"
            << read1->core.qual << "\t"
            << read2->core.l_qseq << "\t"
            << read2->core.qual << "\t"
            << get_seq(read1) << "\t"
            << get_qual(read1) << "\t"
            << get_seq(read2) << "\t"
            << get_qual(read2) << std::endl;
    }
}


// Return whether a read pair represents a proper sample
// Proper sample is defined as having both pairs on the same chromosome,
// both mapped, and having FR orientation
bool is_proper_sample(bam1_t *read1, bam1_t *read2, int max_insert=MAX_INSERT) {
    int size1 = abs(read1->core.isize);
    int size2 = abs(read2->core.isize);
    return ((size1 > 0) && (size1 < max_insert) &&
        (size2 > 0) && (size2 < max_insert) &&
        pair_is_mapped(read1, read2) &&
    // return pair_is_mapped(read1, read2) &&
        pair_is_FR(read1, read2) &&
        (read1->core.tid == read2->core.tid));
}

bool is_proper_RF(bam1_t *read1, bam1_t *read2, int max_insert=MAX_INSERT) {
    int size1 = abs(read1->core.isize);
    int size2 = abs(read2->core.isize);
    return ((size1 > 0) && (size1 < max_insert) &&
        (size2 > 0) && (size2 < max_insert) &&
        pair_is_mapped(read1, read2) &&
        pair_is_RF(read1, read2) &&
        (read1->core.tid == read2->core.tid));
}


// Samples the first ten million reads in the input bamfile
// in order to generate read statistics used for calling variants
void sample_statistics(samfile_t *samfile, int *max_insert,
    int threshold, int sample_max_insert=MAX_INSERT) {

    int sample_size = 1000000;
    int sample_count = 0;
    int proper_sample_count = 0;

    bam1_t *read1 = bam_init1();
    bam1_t *read2 = bam_init1();

    RunningMean longread_sample_mean;
    RunningMean shortread_sample_mean;

    RunningMedian longread_sample_median(500, 20000);
    RunningMedian shortread_sample_median(1, sample_max_insert);

    while (sample_count < sample_size) {

        if (samread(samfile, read1) > 0) {
            samread(samfile, read2);
            sample_count += 2;
        } else {
            break;
        }

        if (is_proper_sample(read1, read2, sample_max_insert)) {
            proper_sample_count += 2;
            longread_sample_mean.push(abs(read1->core.isize));
            longread_sample_median.push(abs(read1->core.isize));
        } else if (is_proper_RF(read1, read2, sample_max_insert)) {
            shortread_sample_mean.push(abs(read1->core.isize));
            shortread_sample_median.push(abs(read1->core.isize));
        }
    }

    *max_insert = (threshold * longread_sample_median.mad()) +
        longread_sample_median.get_median();

    std::cout << "Sampled " << sample_count << "; " << proper_sample_count <<
        " in proper pairs" << std::endl;
    std::cout << std::endl;

    std::cout << "Sample proper FR mean insert size: "
        << longread_sample_mean.mean << std::endl;
    std::cout << "Sample proper FR standard deviation: "
        << longread_sample_mean.std_dev << std::endl;
    std::cout << std::endl;

    std::cout << "Sample proper FR median insert size: "
        << longread_sample_median.get_median() << std::endl;
    std::cout << "Sample proper FR median abs deviation: "
        << longread_sample_median.mad() << std::endl;
    std::cout << "Acceptable max FR insert size: "
        << *max_insert << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "Sample RF mean insert size: "
        << shortread_sample_mean.mean << std::endl;
    std::cout << "Sample RF standard deviation: "
        << shortread_sample_mean.std_dev << std::endl;
    std::cout << "Sample RF median insert size: "
        << shortread_sample_median.get_median() << std::endl;
    std::cout << "Sample proper RF median abs deviation: "
        << shortread_sample_median.mad() << std::endl;

    std::cout << std::endl;

    bam_destroy1(read1);
    bam_destroy1(read2);
}

// Calls which type of variant (if any) a pair of reads represents
// and updates samfile stats
VariantType call_variant(bam1_t *read1, bam1_t *read2, int max_insert) {

    int insert = abs(read1->core.isize);

    // Call a proper sample.
    // Pair is FR, mapped, both reads map to same chromosome,
    // and insert size is less than max.
    //if (is_proper_sample(read1, read2, max_insert)) {
    //    return NORMAL;

    if (read1->core.flag & BAM_FPROPER_PAIR) {
        return NORMAL;

    // Otherwise, call an improper sample
    } else {
        // Mapped pairs
        if (pair_is_mapped(read1, read2)) {
            // Reads map to different chromosomes
            // Call translocation
            if (read1->core.tid != read2->core.tid) {
                return TRANSLOCATION;

            // Both reads map to same chromosome
            } else {
                // Each read has same orientation
                // Represents inversion
                if ((read1->core.flag & BAM_FREVERSE) ==
                    (read2->core.flag & BAM_FREVERSE)) {
                    return INVERSION;

                // Insert size is greater than
                // median + (INSERT_DEVIATIONS * mad)
                // Uses same median and mad criteria as deletion calling
                // Since
                } else if (pair_is_RF(read1, read2)) {
                    // Insertions will have similar insert size
                    // characteristics as deletions, but inverted orientation
                    if (insert > max_insert) {
                        return INSERTION;

                    // Insert size less than that of a deletion
                    // Pair fits criteria for short read
                    } else {
                        return SHORT;
                    }

                // Large insert size exceeding our threshold
                // represents a deletion
                // (defined as n*mad + median, n specified by user)
                } else if (pair_is_FR(read1, read2)) {
                    return DELETION;
                }
            }

        // Unmapped pairs
        } else if ((read1->core.flag & BAM_FUNMAP) !=
            (read2->core.flag & BAM_FUNMAP)) {
            return ONE_UNMAPPED;
        } else {
            return BOTH_UNMAPPED;
        }
    }
}

int main(int argc, char* argv[]) {

    /* Command line arguments and options */
    int c, length;
    int OPT_bam_input = 0;
    int OPT_write_unmapped = 0;
    int INSERT_DEVIATIONS = 5;
    char *bam_filename = NULL;
    char *output_directory = NULL;
    int sample_max_insert = MAX_INSERT;

    int max_insert;

    /* Data Structures */
    RunningMean longread_mean;
    RunningMean shortread_mean;

    AlignmentStats *stats = new AlignmentStats;
    // RunningMedian longread_median;
    // RunningMedian shortread_median;

    VariantType type;
    Orientation orientation;

    /* Input samfile */
    samfile_t *samfile;

    /* Samfile outputs */
    samfile_t *deletions;
    samfile_t *insertions;
    samfile_t *inversions;
    samfile_t *translocations;
    samfile_t *unmapped;
    samfile_t *both_unmapped;

    /* Plaintext outputs */
    std::ofstream deletion_pairs;
    std::ofstream insertion_pairs;
    std::ofstream inversion_pairs;
    std::ofstream transloc_pairs;
    std::ofstream unmapped_pairs;
    std::ofstream both_unmapped_pairs;

    /* SAM read structures */
    bam1_t *read1 = bam_init1();
    bam1_t *read2 = bam_init1();
    bam1_t *temp;

    /* Define command line arguments */
    const struct option longopts[] = {
        { "help",              no_argument,       NULL, 'h' },
        { "input-file",        required_argument, NULL, 'i' },
        // { "output-dir",        required_argument, NULL, 'o' },
        { "bam-input",         no_argument,       NULL, 'b' },
        { "save-unmapped",     no_argument,       NULL, 'u' },
        { "insert-deviations", required_argument, NULL, 'd' },
        { "max-insert",        required_argument, NULL, 'm' },
        { 0, 0, 0, 0 }
    };

    /* If no arguments are passed print usage information */
    if (argc == 1) {
        usage();
        exit(0);
    }

    /* Get command line arguments */
    while ((c = getopt_long(argc, argv, "hi:o:bud:m:", longopts, NULL)) != -1) {
        switch (c) {
            /* Print help */
            case 'h':
                usage();
                exit(0);
                break;
            /* Get name of input file */
            case 'i':
                bam_filename = optarg;
                break;

            /* Get name of output directory */
            // case 'o':
            //     output_directory = optarg;
            //     break;

            /* Read .bam format instead of .sam */
            case 'b':
                OPT_bam_input = 1;
                break;

            /* Save unmapped read pairs for re-mapping */
            case 'u':
                OPT_write_unmapped = 1;
                break;

            /* Number of standard deviations to use as threshold
               for acceptable insert size */
            case 'd':
                INSERT_DEVIATIONS = atoi(optarg);
                break;

            /* Maximum acceptable insert size for a sampled pair */
            case 'm':
                sample_max_insert = atoi(optarg);
                break;

            case 0:
                break;

            /* Exit if presented with an invalid flag */
            default:
                std::cerr << "bamstat: Use --help for more information.";
                std::cerr << std::endl;
                exit(-1);
        }
    }

    std::cout << "Std deviations: " << INSERT_DEVIATIONS << std::endl;
    std::cout << "Bam input?: " << OPT_bam_input << std::endl;
    std::cout << "Write unmapped?: " << OPT_write_unmapped << std::endl;

    /* Require both an input file and an output directory */
    if (bam_filename == NULL) {
	bam_filename = "-";
        // std::cerr << "bamstat: No input file specified." << std::endl;
        // std::cerr << "bamstat: Use --help for more information." << std::endl;
        // exit(-1);
    }

    // if (output_directory == NULL) {
    //     std::cerr << "bamstat: No output directory specified." << std::endl;
    //     std::cerr << "bamstat: Use --help for more information." << std::endl;
    //     exit(-1);
    // }

    // if (chdir(output_directory) != 0) {
    //     std::cerr << "bamstat: Output directory does not exist." << std::endl;
    //     exit(-1);
    // }

    /* Open input file as appropriate type */
    if (OPT_bam_input) {
        samfile = samopen(bam_filename, "rb", NULL);
    } else {
        samfile = samopen(bam_filename, "r", NULL);
    }

    // Sample first million reads from samfile in order to calculate
    // insert size that is likely to represent a deletion
    sample_statistics(samfile, &max_insert, INSERT_DEVIATIONS,
        sample_max_insert);
    samclose(samfile);

    // Reset samfile to beginning
    if (OPT_bam_input) {
        samfile = samopen(bam_filename, "rb", NULL);
    } else {
        samfile = samopen(bam_filename, "r", NULL);
    }

    // Open output samfiles
    deletions = samopen("deletions.bam", "wb", samfile->header);
    insertions = samopen("insertions.bam", "wb", samfile->header);
    inversions = samopen("inversions.bam", "wb", samfile->header);
    translocations = samopen("translocations.bam", "wb", samfile->header);

    // Open output plaintext files
    deletion_pairs.open("deletion_pairs.txt");
    insertion_pairs.open("insertion_pairs.txt");
    inversion_pairs.open("inversion_pairs.txt");
    transloc_pairs.open("translocation_pairs.txt");

    if (OPT_write_unmapped) {
        unmapped = samopen("unmapped.bam", "wb", samfile->header);
        unmapped_pairs.open("unmapped_pairs.txt");
        both_unmapped = samopen("both_unmapped.bam", "wb", samfile->header);
        both_unmapped_pairs.open("both_unmapped_pairs.txt");
    }

    // Call each pair in samfile
    while (samread(samfile, read1) > 0) {
        samread(samfile, read2);

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

        // Call variant type and update AlignmentStats accordingly
        type = call_variant(read1, read2, max_insert);
        orientation = get_orientation(read1, read2);

        stats->update(type, abs(read1->core.isize), orientation);

        // Write each pair to respective file
        switch (type) {
            case DELETION:
                samwrite(deletions, read1);
                samwrite(deletions, read2);
                write_pair(deletion_pairs, read1, read2, samfile);
                break;

            case INSERTION:
                samwrite(insertions, read1);
                samwrite(insertions, read2);
                write_pair(insertion_pairs, read1, read2, samfile);
                break;

            case INVERSION:
                samwrite(inversions, read1);
                samwrite(inversions, read2);
                write_pair(inversion_pairs, read1, read2, samfile);
                break;

            case TRANSLOCATION:
                samwrite(translocations, read1);
                samwrite(translocations, read2);
                write_pair(transloc_pairs, read1, read2, samfile);
                break;

            case ONE_UNMAPPED:
                if (OPT_write_unmapped) {
                    samwrite(unmapped, read1);
                    samwrite(unmapped, read2);
                    write_pair(unmapped_pairs, read1, read2, samfile);
                }
                break;

            case BOTH_UNMAPPED:
                if (OPT_write_unmapped) {
                    samwrite(both_unmapped, read1);
                    samwrite(both_unmapped, read2);
                    write_pair(both_unmapped_pairs, read1, read2, samfile);
                }
                break;

            // NORMAL, SHORT_RF, BOTH_UNMAPPED
            default:
                break;
        }
    }

    // Write out statistic results to file
    stats->print("stats.file");

    /* Clean up */
    samclose(samfile);
    samclose(deletions);
    samclose(insertions);
    samclose(inversions);
    samclose(translocations);

    deletion_pairs.close();
    insertion_pairs.close();
    inversion_pairs.close();
    transloc_pairs.close();

    if (OPT_write_unmapped) {
        samclose(unmapped);
        unmapped_pairs.close();
    }

    delete stats;

    bam_destroy1(read1);
    bam_destroy1(read2);

    return 0;
}

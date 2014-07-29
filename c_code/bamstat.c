#include "bamstat.h"

/* Print usage information */
int usage() {
    printf("\n");
    printf("Program: bamstat (Read statistics and sorting)\n");
    printf("Version: 0.2.0\n\n");
    printf("Usage: bamstat -i <filename> -o <directory> [options]\n\n");
    printf("Options:\n");
    printf("    -h --help\t\t\tShow this message\n");
    printf("    -i --input-file <filename>\tRead in from <filename>\n");
    printf("    -o --output-dir <directory>\tStore output in <directory>\n");
    printf("    -s --insert-deviations <threshold>\tDefine normal insert size "
        "to be within\n\t\t\t\t<threshold> standard deviations of the mean\n");
    printf("    -b --bam-input\t\tRead from .bam format instead "
        "of .sam\n");
    printf("    -u --save-unmapped\t\tOutput unmapped read pairs\n\n");

    return 0;
}

int is_proper_sample(bam1_t *read1, bam1_t *read2) {
    int max_insert = 20000;
    return ((abs(read1->core.isize) > 0) &&
            (abs(read1->core.isize) < max_insert) &&
            (abs(read2->core.isize) > 0) &&
            (abs(read2->core.isize) < max_insert) &&
            !(read1->core.flag & BAM_FUNMAP) &&
            !(read1->core.flag & BAM_FUNMAP) &&
            (read1->core.tid == read2->core.tid));

}

int sample_statistics(samfile_t *samfile){
    int sample_size = 1000000;
    int sample_count = 0;
    int max_per_tile = 20;
    int current_tile_count = 0;
    int proper_sample_count = 0;

    char *last_tile = NULL;
    char *current_tile;

    bam1_t *read1 = bam_init1();
    bam1_t *read2 = bam_init1();

    RunningStat *longread_sample = malloc(sizeof(RunningStat));
    MALLOC_CHECK(longread_sample);
    BZERO(longread_sample, sizeof(RunningStat));
    reset_running(longread_sample);

    RunningStat *shortread_sample = malloc(sizeof(RunningStat));
    MALLOC_CHECK(shortread_sample);
    BZERO(shortread_sample, sizeof(RunningStat));
    reset_running(shortread_sample);

   // while (sample_count < sample_size) {
        if (samread(samfile, read1) > 0) {
            samread(samfile, read2);
        } else {
            // break;
        }

        current_tile = strtok(bam1_qname(read1),":");
        current_tile = strtok(NULL, ":");
        current_tile = strtok(NULL, ":");
        current_tile = strtok(NULL, ":");

        if (last_tile == NULL) {
            last_tile = current_tile;
        }

        if (strcmp(current_tile, last_tile) == 0) {
            printf("same tile\n");
            printf("%s\t%s\n", current_tile, last_tile);
            if (current_tile_count < max_per_tile) {
                sample_count += 2;
                current_tile_count += 2;
            }
        } else {
            printf("different tile\n");
            printf("%s\n", current_tile);
            printf("%s\n", last_tile);
            printf("%d\t%d\n", current_tile_count, sample_count);
            last_tile = current_tile;
            current_tile_count = 0;
        }
     /*   for (int i = 0; i < 3; i++) {
            last_tile = strtok(bam1_qname(read1))
        }*/

        //sample_count
 //   }
    bam_destroy1(read1);
    bam_destroy1(read2);
    SAFE_FREE(longread_sample);
    SAFE_FREE(shortread_sample);

}

int main(int argc, char* argv[]) {

    /* Command line arguments and options */
    int c, length;
    int OPT_bam_input = 0;
    int OPT_write_unmapped = 0;
    int INSERT_DEVIATIONS = 3;
    char *bam_filename = NULL;
    char *output_directory = NULL;

    /* Statistic data structures */
    AlignmentStats *stats = malloc(sizeof(AlignmentStats));
    MALLOC_CHECK(stats);
    BZERO(stats, sizeof(AlignmentStats));
    init_stats(stats);

    RunningStat *longread_sizes = malloc(sizeof(RunningStat));
    MALLOC_CHECK(longread_sizes);
    BZERO(longread_sizes, sizeof(RunningStat));
    reset_running(longread_sizes);

    RunningStat *shortread_sizes = malloc(sizeof(RunningStat));
    MALLOC_CHECK(shortread_sizes);
    BZERO(shortread_sizes, sizeof(RunningStat));
    reset_running(shortread_sizes);

    /* Input samfile */
    samfile_t *samfile;

    /* Samfile outputs */
    samfile_t *deletions;
    samfile_t *insertions;
    samfile_t *inversions;
    samfile_t *translocations;

    /* Plain text file outputs */
    FILE *deletion_pairs;
    FILE *insertion_pairs;
    FILE *inversion_pairs;
    FILE *translocation_pairs;

    /* SAM read structures */
    bam1_t *read1 = bam_init1();
    bam1_t *read2 = bam_init1();

    /* Define command line arguments */
    static struct option longopts[] = {
        { "help",              no_argument,       NULL, 'h' },
        { "input-file",        required_argument, NULL, 'i' },
        { "output-dir",        required_argument, NULL, 'o' },
        { "bam-input",         no_argument,       NULL, 'b' },
        { "save-unmapped",     no_argument,       NULL, 'u' },
        { "insert-deviations", required_argument, NULL, 's' },
        { 0, 0, 0, 0 }
    };

    /* If no arguments are passed print usage information */
    if (argc == 1) {
        usage();
        exit(0);
    }

    /* Get command line arguments */
    while ((c = getopt_long(argc, argv, "i:hbuo:s:", longopts, NULL)) != -1) {
        switch (c) {
            /* Print help */
            case 'h':
                usage();
                exit(0);
                break;
            /* Get name of input file */
            case 'i':
                length = strlen(optarg) + 1;
                bam_filename = malloc(length * sizeof(char));
                MALLOC_CHECK(bam_filename);
                BZERO(bam_filename, (length * sizeof(char)));
                strcpy(bam_filename, optarg);
                break;

            /* Get name of output directory */
            case 'o':
                length = strlen(optarg) + 1;
                output_directory = malloc(length * sizeof(char));
                MALLOC_CHECK(output_directory);
                BZERO(output_directory, (length * sizeof(char)));
                strcpy(output_directory, optarg);
                break;

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
            case 's':
                INSERT_DEVIATIONS = atoi(optarg);
                break;

            case 0:
                break;

            /* Exit if presented with an invalid flag */
            default:
                fprintf(stderr, "bamstat: Use --help for more information.\n");
                exit(-1);
        }
    }

    /* Require both an input file and an output directory */
    if (bam_filename == NULL) {
        fprintf(stderr, "bamstat: No input file specified.\n");
        fprintf(stderr, "bamstat: Use --help for more information.\n");
        exit(-1);
    }

    if (output_directory == NULL) {
        fprintf(stderr, "bamstat: No output directory specified.\n");
        fprintf(stderr, "bamstat: Use --help for more information.\n");
        exit(-1);
    }

    /* Open input file as appropriate type */
    if (OPT_bam_input) {
        samfile = samopen(bam_filename, "rb", NULL);
    } else {
        samfile = samopen(bam_filename, "r", NULL);
    }

    sample_statistics(samfile);

    /* Open output files using input header as template */
    deletions = samopen("deletions.bam", "wb", samfile->header);
    insertions = samopen("insertions.bam", "wb", samfile->header);
    inversions = samopen("inversions.bam", "wb", samfile->header);
    translocations = samopen("translocations.bam", "wb", samfile->header);


    /* Clean up */
    /* Free data structures */
    SAFE_FREE(output_directory);
    SAFE_FREE(bam_filename);
    SAFE_FREE(stats);
    SAFE_FREE(longread_sizes);
    SAFE_FREE(shortread_sizes);
    bam_destroy1(read1);
    bam_destroy1(read2);

    /* Close files */
    samclose(samfile);
    samclose(deletions);
    samclose(insertions);
    samclose(inversions);
    samclose(translocations);

    return 0;
}
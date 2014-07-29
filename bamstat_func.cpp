int bamstat(samfile_t *samfile, ) {

    /* Samfile outputs */
    samfile_t *deletions = samopen("deletions.bam", "wb", samfile->header);
    samfile_t *insertions = samopen("insertions.bam", "wb", samfile->header);
    samfile_t *inversions = samopen("inversions.bam", "wb", samfile->header);
    samfile_t *translocations = samopen("translocations.bam", "wb", samfile->header);
    samfile_t *unmapped = samopen("unmapped.bam", "wb", samfile->header);

    /* Plaintext outputs */
    std::ofstream deletion_pairs("deletions.txt");
    std::ofstream insertion_pairs("insertions.txt");
    std::ofstream inversion_pairs("inversions.txt");
    std::ofstream transloc_pairs("translocations.txt");
    std::ofstream unmapped_reads("unmapped.txt");

    /* SAM read structures */
    bam1_t *read1 = bam_init1();
    bam1_t *read2 = bam_init1();
    bam1_t *temp;

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

        type = call_variant(read1, read2, stats, FR_max_insert);

        // Write each pair to respective file
        switch (type) {
            case INVERSION:
                samwrite(inversions, read1);
                samwrite(inversions, read2);
                write_pair(inversion_pairs, read1, read2, samfile);
                break;
            case DELETION:
                samwrite(deletions, read1);
                samwrite(deletions, read2);
                write_pair(deletion_pairs, read1, read2, samfile);
                break;
            case TRANSLOCATION:
                samwrite(translocations, read1);
                samwrite(translocations, read2);
                write_pair(transloc_pairs, read1, read2, samfile);
                break;
            case INSERTION:
                samwrite(insertions, read1);
                samwrite(insertions, read2);
                write_pair()
            case NORMAL:
                break;
            default:
                break;
        }
    }
}
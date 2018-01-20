#include <stdlib.h>
#include <stdio.h>
#include <plinkio/plinkio.h>

int
main(int argc, char *argv[])
{
    if (argc < 2)
    {
        printf("Usage:\n");
        printf("    plinkio <bfile>\n");
        printf("where <bfile> is path to plink bed/bim/fam file wihtout extension\n");
        return EXIT_FAILURE;
    }

    struct pio_file_t plink_file;
    snp_t *snp_buffer;
    int sample_id;
    int locus_id;

    if( pio_open( &plink_file, argv[ 1 ] ) != PIO_OK )
    {
        printf( "Error: Could not open %s\n", argv[ 1 ] );
        return EXIT_FAILURE;
    }

    if( !pio_one_locus_per_row( &plink_file ) )
    {
        printf( "This script requires that snps are rows and samples columns.\n" );
        return EXIT_FAILURE;
    }

    locus_id = 0;
    snp_buffer = (snp_t *) malloc( pio_row_size( &plink_file ) );
    while( pio_next_row( &plink_file, snp_buffer ) == PIO_OK )
    {
        for( sample_id = 0; sample_id < pio_num_samples( &plink_file ); sample_id++)
        {
            struct pio_sample_t *sample = pio_get_sample( &plink_file, sample_id );
            struct pio_locus_t *locus = pio_get_locus( &plink_file, locus_id );
            printf( "Individual %s has genotype %d for snp %s.\n", sample->iid, snp_buffer[ sample_id ], locus->name );
        }

        locus_id++;
    }

    free( snp_buffer );
    pio_close( &plink_file );
    
    return EXIT_SUCCESS;
}


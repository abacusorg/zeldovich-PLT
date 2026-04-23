#include "fft_wisdom.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
Each rank calls FFTW_IMPORT_WISDOM_FROM_FILENAME(FFTW_WISDOM_FILENAME) with the same string, from same file on shared disk. There is no MPI broadcast; the disk file is the shared copy!
*/

#ifndef FFTW_WISDOM_PREVIEW_MAX
#define FFTW_WISDOM_PREVIEW_MAX 512u // cap on printed wisdom string (during debug)
#endif

void fft_wisdom_import_from_file(int rank)
{
    const int imported = FFTW_IMPORT_WISDOM_FROM_FILENAME(FFTW_WISDOM_FILENAME);

    if (!imported) {
        FFTW_FORGET_WISDOM();
        printf("[FFTW-WISDOM] Rank %d: no usable wisdom at '%s' (%s)\n",
               rank, FFTW_WISDOM_FILENAME, PRECISION_NAME);
        fflush(stdout);
    } else {
        printf("[FFTW-WISDOM] Rank %d: imported wisdom from '%s' (%s precision).\n", rank,
               FFTW_WISDOM_FILENAME, PRECISION_NAME);
        fflush(stdout);
    }
}

void fft_wisdom_export_rank0(int rank)
{
    if (rank != 0) {
        return;
    }

    char *pre = FFTW_EXPORT_WISDOM_TO_STRING();
    if (pre != NULL) {
        const size_t n = strlen(pre);
        printf("[FFTW-WISDOM] Rank 0: pre-export wisdom length %zu bytes (%s precision)\n", n,
               PRECISION_NAME);
        if (n <= FFTW_WISDOM_PREVIEW_MAX) {
            printf("[FFTW-WISDOM] Rank 0: wisdom (full):\n%.*s\n", (int)n, pre);
        } else {
            printf("[FFTW-WISDOM] Rank 0: wisdom (first %u bytes):\n%.*s\n... (%zu more bytes)\n",
                   (unsigned)FFTW_WISDOM_PREVIEW_MAX, (int)FFTW_WISDOM_PREVIEW_MAX, pre,
                   n - (size_t)FFTW_WISDOM_PREVIEW_MAX);
        }
        fflush(stdout);
        FFTW_FREE(pre);
    } else {
        fprintf(stderr,
                "[FFTW-WISDOM] Rank 0: FFTW_EXPORT_WISDOM_TO_STRING returned NULL (%s precision)\n",
                PRECISION_NAME);
        fflush(stderr);
    }

    const int ok = FFTW_EXPORT_WISDOM_TO_FILENAME(FFTW_WISDOM_FILENAME);
    if (!ok) {
        fprintf(stderr,
                "[FFTW-WISDOM] Rank 0: failed to export wisdom to '%s' (%s precision)\n",
                FFTW_WISDOM_FILENAME, PRECISION_NAME);
        fflush(stderr);
    } else {
        printf("[FFTW-WISDOM] Rank 0: exported wisdom to '%s' (%s precision)\n",
               FFTW_WISDOM_FILENAME, PRECISION_NAME);
        fflush(stdout);
    }
}

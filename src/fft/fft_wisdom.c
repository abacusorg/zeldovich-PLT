#include "fft_wisdom.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void fft_wisdom_import_broadcast(int rank, MPI_Comm comm)
{
    int imported = 0;
    char *wisdom_str = NULL;
    size_t wisdom_len = 0;

    if (rank == 0) {
        imported = FFTW_IMPORT_WISDOM_FROM_FILENAME(FFTW_WISDOM_FILENAME);

        if (!imported) {
            // Ensure we start from a clean state on all ranks.
            FFTW_FORGET_WISDOM();
            printf("[FFTW-WISDOM] No existing '%s' for %s precision; "
                   "planning from scratch, will create it later.\n",
                   FFTW_WISDOM_FILENAME, PRECISION_NAME);
            fflush(stdout);
        } else {
            printf("[FFTW-WISDOM] Imported wisdom from '%s' (%s precision).\n",
                   FFTW_WISDOM_FILENAME, PRECISION_NAME);
            fflush(stdout);
        }

        // Export current (possibly empty) wisdom state to a string to share
        // exactly the same wisdom with all ranks.
        wisdom_str = FFTW_EXPORT_WISDOM_TO_STRING();
        if (wisdom_str != NULL) {
            wisdom_len = strlen(wisdom_str);
        } else {
            wisdom_len = 0;
        }
    }

    // Broadcast the length first.
    MPI_Bcast(&wisdom_len, 1, MPI_UNSIGNED_LONG_LONG, 0, comm);

    if (wisdom_len > 0) {
        if (rank != 0) {
            wisdom_str = (char *)malloc(wisdom_len + 1);
            if (!wisdom_str) {
                fprintf(stderr,
                        "[FFTW-WISDOM] Rank %d: failed to allocate buffer for wisdom broadcast\n",
                        rank);
                MPI_Abort(comm, 1);
            }
        }

        MPI_Bcast(wisdom_str, (int)wisdom_len, MPI_CHAR, 0, comm);
        wisdom_str[wisdom_len] = '\0';

        // All ranks (including 0) import the shared wisdom string.
        FFTW_FORGET_WISDOM();
        FFTW_IMPORT_WISDOM_FROM_STRING(wisdom_str);

        if (rank == 0) {
            printf("[FFTW-WISDOM] Broadcasted wisdom to all ranks (len=%zu).\n",
                   wisdom_len);
            fflush(stdout);
        }
    } else {
        // No wisdom available; all ranks already have a clean state.
        if (rank == 0) {
            printf("[FFTW-WISDOM] No wisdom to broadcast; all ranks start fresh.\n");
            fflush(stdout);
        }
    }

    if (wisdom_str != NULL) {
        // Per FFTW docs, the string from FFTW_EXPORT_WISDOM_TO_STRING is
        // allocated with fftw_malloc/fftwf_malloc. Use FFTW_FREE for it.
        FFTW_FREE(wisdom_str);
    }
}

void fft_wisdom_export_rank0(int rank)
{
    if (rank != 0) {
        return;
    }

    int ok = FFTW_EXPORT_WISDOM_TO_FILENAME(FFTW_WISDOM_FILENAME);
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


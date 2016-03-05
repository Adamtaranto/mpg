/*
 * ============================================================================
 *
 *       Filename:  libmpg.cc
 *    Description:  Utilities for mpg
 *        License:  GPLv3+
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */


#include <libmpg.hh>

char NT[] = "ACGT";

std::string
random_dna(size_t size, uint64_t seed)
{
    std::mt19937_64 rand(seed);
    std::string seq;

    for (size_t i = 0; i < (size / 32) + 1; i++) {
        uint64_t r = rand();
        for (size_t j = 0; j < 32; j++) {
            // We use the j'th 2 LSBs to determine which of the 4 nucelotides
            // to chose. The following just extracts this as an index.
            size_t nt = (r >> (j * 2)) & 0x03;
            seq += NT[nt];
        }
    }
    // We normally generate more than we need, so nuke the rest.
    seq.erase(size);
    return seq;
}


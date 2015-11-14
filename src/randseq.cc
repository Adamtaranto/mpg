/*
 * ============================================================================
 *
 *       Filename:  randseq.cc
 *    Description:  mersenne twister DNA seqs
 *        License:  GPLv3+
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */


#include <libmpg.hh>

int
main (int argc, char *argv[])
{
    std::string s = random_dna(1<<20, 1);
    std::cout << ">1\n"<< s << std::endl;
    return EXIT_SUCCESS;
}

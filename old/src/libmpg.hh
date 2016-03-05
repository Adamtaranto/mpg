/*
 * ============================================================================
 *
 *       Filename:  libmpg.hh
 *    Description:  Utilities for mpg
 *        License:  GPLv3+
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */

#ifndef LIBMPG_HH
#define LIBMPG_HH

#include <random>
#include <string>
#include <iostream>
#include <map>
#include <cassert>


extern char NT[];

std::string random_dna(size_t size, uint64_t seed);


#endif /* LIBMPG_HH */

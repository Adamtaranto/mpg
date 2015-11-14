/*
 * ============================================================================
 *
 *       Filename:  mpg-burnin.cc
 *    Description:  markov process with burnin
 *        License:  GPLv3+
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */


#include <libmpg.hh>

#define runif(x) std::generate_canonical<double, 64>(x)
typedef std::map<char, double> base_prob;
typedef std::map<std::string, base_prob> transition_prob;

char
weighted_rand_base(base_prob &probs, double rand)
{
    double culm = 0.;
    for (const auto &bp: probs) {
        culm += bp.second;
        if (culm > rand) {
            return bp.first;
        }
    }
    return 'N';
}

void
normalise_transition_probs(transition_prob &transprob) {
    for (auto &pair: transprob) {
        double sum = 0;
        for (auto &bases: pair.second) {
            sum += bases.second;
        }
        for (auto &bases: pair.second) {
            bases.second /= sum;
        }
    }
}

std::string
mpg_burnin(size_t size, size_t k, transition_prob transprob, size_t burnin_nts,
           uint64_t seed)
{
    std::mt19937_64 rand(seed);
    std::string seq = random_dna(k, seed);

    normalise_transition_probs(transprob);

    for (auto &pair: transprob) {
        std::cout << "\t" << pair.first;
    }
    std::cout << std::endl;
    for (auto &pair: transprob) {
        for (auto &bases: pair.second) {
            std::cout << "\t" << pair.first << "\t" << bases.second << std::endl;
        }
    }

    for (size_t i = k; i < size + burnin_nts + k; i++) {
        std::string last = seq.substr(i-k);
        seq += weighted_rand_base(transprob[last], runif(rand));
    }

    seq.erase(0, k + burnin_nts);
    return seq;
}

int
main (int argc, char *argv[])
{
    /*
    const transition_prob p {
        {'A', {{'A', 0.1}, {'C', 0.4}, {'G', 0.4}, {'T', 0.1}}},
        {'C', {{'A', 0.4}, {'C', 0.1}, {'G', 0.1}, {'T', 0.4}}},
        {'G', {{'A', 0.4}, {'C', 0.1}, {'G', 0.1}, {'T', 0.4}}},
        {'T', {{'A', 0.1}, {'C', 0.4}, {'G', 0.4}, {'T', 0.1}}},
    };
    */

    std::mt19937_64 rand(123);
    transition_prob p;
    for (size_t i = 0; i< 4; i++) {
        std::string s(1, NT[i]);
        base_prob bp;
        for (size_t j = 0; j < 4; j++) {
            bp[NT[j]] = runif(rand);
        }
        p[s] = bp;
    }

    std::string s = mpg_burnin(128, 1, p, 1000, 1);
    std::cout << ">1\n"<< s << std::endl;
    return EXIT_SUCCESS;
}

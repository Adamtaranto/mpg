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

#include <iterator>

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

void
print_tp(transition_prob &tp)
{
    for (auto &pair: tp) {
        for (auto &bases: pair.second) {
            std::cerr << "\t" << bases.first;
        }
        break;
    }
    std::cerr << std::endl;
    for (auto &pair: tp) {
        std::cerr << pair.first;
        for (auto &bases: pair.second) {
            std::cerr << "\t" << bases.second;
        }
        std::cerr << std::endl;
    }
}

std::string
mpg_burnin(size_t size, size_t k, transition_prob transprob, size_t burnin_nts,
           uint64_t seed)
{
    std::mt19937_64 rand(seed);
    std::string seq = random_dna(k, seed);

    normalise_transition_probs(transprob);

    print_tp(transprob);

    for (size_t i = k; i < size + burnin_nts + k; i++) {
        std::string last = seq.substr(i-k);
        seq += weighted_rand_base(transprob[last], runif(rand));
    }

    seq.erase(0, k + burnin_nts);
    return seq;
}

class KmerFactory
{
public:
    typedef std::iterator<std::random_access_iterator_tag,
                          std::string> iterator;

    KmerFactory(size_t k, size_t alpha=4)
        : K(k)
        , alphasize(alpha)
        , n_k(pow(alphasize, K))
        , state(0)
    {}


protected:
    std::string
    n2k             (size_t n)
    {
        std::string kmer;
        // Bits are reversed, so go backwards
        for (size_t i = K; i != 0; i--) {
            // for each 2 bits of the kmer's number, i.e. for each base.
            size_t nt = (n >> (i - 1) * 2) & 0x03;
            kmer += NT[nt];
        }
        return kmer;
    }

private:
    size_t K;
    size_t alphasize;
    size_t n_k;
    size_t state;
};

std::vector<std::string>
all_kmers_old(size_t k)
{
        size_t n_k = pow(4, k);
    std::vector<std::string> kmers(n_k);

    for (size_t i = k; i > 0; i--) {
        size_t each_nt = pow(4, i - 1);
        size_t reps = pow(4, k - i);
        for (size_t r = 0; r < reps; r++) {
            for (size_t nt = 0; nt < 4; nt++) {
                for (size_t j = 0; j < each_nt; j++) {
                    // This index is quite complex, let's break it down:
                    //   Our chunk (rep) starts at 4 * number of nucs per rep,
                    //     which is (r * 4 * each_nt).
                    //   Our position in the chunk is (nt * each_nt + j), i.e.
                    //     nucl block `nt` * how many of each nuc, + j for the
                    //     acutal position of this one.
                    //   The max value of the second half (nt * each_nt + j) is
                    //     (4 * each_nt - 1), i.e. there are (4 * each_nt)
                    //     values per rep.
                    size_t idx =  r * 4 * each_nt + nt * each_nt + j;
                    kmers[idx] += NT[nt];
                }
            }
        }
    }
    return kmers;
}

std::vector<std::string>
all_kmers(size_t k)
{
    size_t n_k = pow(4, k);
    std::vector<std::string> kmers(n_k);

    for (size_t i = 0; i < n_k; i++) {  // for each kmer
        for (ssize_t j = k - 1; j >= 0; j--) {
            // for each 2 bits of the kmer's number, i.e. for each base.
            size_t nt = (i >> j * 2) & 0x03;
            kmers[i] += NT[nt];
        }
    }
    return kmers;
}

int
main (int argc, char *argv[])
{
    std::mt19937_64 rand(123);
    transition_prob p;
    size_t K=4;
    auto kmers = all_kmers(K);
    for (auto &s: kmers) {
        std::cout << s << std::endl;
        base_prob bp;
        for (size_t j = 0; j < 4; j++) {
            bp[NT[j]] = runif(rand);
        }
        p[s] = bp;
    }

    //std::string s = mpg_burnin(2500000, K, p, 100000, (size_t)(void *)(argv));
    //std::cout << ">1\n"<< s << std::endl;
    return EXIT_SUCCESS;
}

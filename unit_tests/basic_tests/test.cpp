//
// Created by Wulfix on 16/11/2017.
//

#include "HubbardClass/include/bitset.h"

#include <boost/math/special_functions.hpp>


/** Give the next bitset permutation in a lexicographical sense.
 *
 *      examples:
 *          011 -> 101
 *          101 -> 110
 *
 * Taken from (http://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation)
 * Credits: Dario Sneidermanis of Argentina, who provided this on November 28, 2009.
 */
boost::dynamic_bitset<> next_bitset_permutation(const boost::dynamic_bitset<>& v) {
    // The given algorithm was constructed for unsigned ints. Since we're accepting boost::dynamic_bitsets, we will first convert to an ulong.
    unsigned long v_ul = v.to_ulong();

    // t gets v's least significant 0 bits set to 1
    unsigned long t = v_ul | (v_ul - 1);
    // Next set to 1 the most significant bit to change and set to 0 the least significant ones, and add the necessary 1 bits.
    unsigned long w = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(v_ul) + 1));

    // Re-convert from ulong to boost::dynamic_bitset
    return boost::dynamic_bitset<> (v.size(), w);
}


/** Generate all bitset permutations in a lexicographical sense, with the smallest being given as an argument
 *
 *      example:
 *          011 -> std::vector {011, 101, 110}
 */
std::vector<boost::dynamic_bitset<>> all_bitset_permutations(const boost::dynamic_bitset<>& v) {
    // First, create an empty set of all possible bitset permutations
    //      The number of different permutations of the 1's in the bitset is equal to
    //          L choose N
    //      with L being the size of the bitset, and N being the number of ones in the bitset
    auto nperm = static_cast<size_t>(boost::math::binomial_coefficient<double>(static_cast<unsigned>(v.size()), static_cast<unsigned>(v.count()))); // total number of permutations
    std::vector<boost::dynamic_bitset<>> permutations (nperm);  // std::vector only accepts unsigned long

    // The first permutation is the smallest lexicographical one
    assert(is_smallest_bitset(v));
    permutations[0] = v;

    for (size_t i = 1; i < nperm; i++) {
        permutations[i] = next_bitset_permutation(permutations[i-1]);
    }

    return permutations;
}


/** Generate the smallest (in lexicographical sense) bitset, given a total length and number of ones (nones)
 *
 *      example:
 *          length=3, nones=2 -> 011
 */
boost::dynamic_bitset<> smallest_bitset(size_t length, size_t nones) {
    boost::dynamic_bitset<> u (length, 0); // Create a zero bitset with length length

    for (size_t i = 0; i < nones; i++) {
        u[i] = true;    // flip the nones least significant bits 0->1
        // operator[0] returns the least significant bit
    }

    return u;
}


/** Checks if the given bitset is the lexicographical smallest one
 */
bool is_smallest_bitset(const boost::dynamic_bitset<>& v) {
    return v == smallest_bitset(v.size(), v.count());
}


/** Concatenate two sets of bitsets, producing every possible combination between both. set1 is major.
 *
 *      example:
 *          concatenate_bitset_sets {000}, {011, 101, 110} -> {000011, 000101, 000110}
 *
 *
 * @param set1
 * @param set2
 * @return
 */
std::vector<boost::dynamic_bitset<>> concatenate_bitset_sets(const std::vector<boost::dynamic_bitset<>>& set1, const std::vector<boost::dynamic_bitset<>>& set2) {
    // Allocate memory for the concatenated sets
    std::vector<boost::dynamic_bitset<>> bitsets (set1.size() * set2.size());

    // i will be the index in the concatenated_bitset_sets
    unsigned long i = 0;
    for (const auto& bitset1 : set1) {
        for (const auto& bitset2: set2) {
            // The new bitset must have the combined length of both bitsets
            boost::dynamic_bitset<> concatenated_bitset (bitset2); // when resizing, the least significant bits stay the same
            concatenated_bitset.resize(bitset1.size() + bitset2.size());

            // Insert bitset1 in big bitset
            for (size_t j = 0; j < bitset1.size(); j++) {
                concatenated_bitset[j + bitset2.size()] = bitset1[j];
            }

            // Add the concatenated bitsets and update the index
            bitsets[i] = concatenated_bitset;
            i++;
        }
    }

    return bitsets;
}


/** Count the number of differences between two bitsets.
 *
 *      example:
 *          count_differences(001, 100) -> 2
 *          count_differences(001, 011) -> 1
 */
unsigned long count_differences(const boost::dynamic_bitset<>& u, const boost::dynamic_bitset<>& v) {
    return (u ^ v).count(); // bitwise XOR
}


/** Slice a given bitset from least significant to most significant bit (i.e. from right to left) This is NOT as they are printed!
 *
 *     example:
 *          u: 0100011101
 *          bitslice(u, 2, 5) is NOT equal to 0001
 *          bitslice(u, 2, 5) => 0100[0111]01 -> 0111
 */
boost::dynamic_bitset<> bitslice(const boost::dynamic_bitset<>& u, unsigned long i, unsigned long j) {
    // The implementation is based on https://stackoverflow.com/a/18316482/7930415.

    if (j < i) {
        throw std::invalid_argument("Check the boundaries of the slicing: given j < i, but it should be j > i.");
    }

    auto slice_length = j - i + 1;  // number of bits the slice will have
    auto shift = i;                 // the number of places the mask will have to shift
    // since we're working with least significant bits, this is just the number of bits to the right of the slice, i.e. the index i

    // Create a mask that will pick out the relevant bits
    boost::dynamic_bitset<> mask (u.size(), (static_cast<unsigned long>(std::pow(2, slice_length)) - 1) << shift);

    // Use the mask on the input bitset, and shift the result to the right i places.
    boost::dynamic_bitset<> v = ((u & mask) >> shift);

    // We're only interested in the slice_length least significant bits
    v.resize(slice_length);
    return v;
}


/** Count the number of double occupancies between two bitsets
 *
 *      example:
 *          count_double_occupancies(00110, 10101) -> 1
 *          count_double_occupancies(1110, 1100) -> 2
 */
unsigned long count_double_occupancies(const boost::dynamic_bitset<>& u, const boost::dynamic_bitset<>& v) {
    return (u & v).count(); // bitwise AND
}

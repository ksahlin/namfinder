//
//  index.hpp
//  cpp_course
//
//  Created by Kristoffer Sahlin on 4/21/21.
//
#ifndef index_hpp
#define index_hpp

#include <chrono>
#include <stdio.h>
#include <string>
#include <vector>
#include <deque>
#include <tuple>
#include <cmath>
#include <iostream>
#include <cassert>
#include "robin_hood.h"
#include "exceptions.hpp"
#include "refs.hpp"
#include "randstrobes.hpp"
#include "indexparameters.hpp"

/*
 * This describes where a randstrobe occurs. Info stored:
 * - reference index
 * - position of the first strobe
 * - offset of the second strobe
*/

class RefRandstrobe {
public:
    RefRandstrobe() { }  // TODO should not be needed
    RefRandstrobe(uint32_t position, uint32_t packed) : position(position), m_packed(packed) {
    }
    uint32_t position;

    int reference_index() const {
        return m_packed >> bit_alloc;
    }

    int strobe2_offset() const {
        return m_packed & mask;
    }

    RefRandstrobeWithHash::packed_t packed() const {
        return m_packed;
    }

private:
    static const int bit_alloc = 8;
    static const int mask = (1 << bit_alloc) - 1;
    RefRandstrobeWithHash::packed_t m_packed;
};

using RefRandstrobeVector = std::vector<RefRandstrobe>;

/*
 * An entry in the randstrobe map that allows retrieval of randstrobe
 * occurrences. To save memory, the entry is either a "direct" or an
 * "indirect" one.
 *
 * - A direct entry is used if the randstrobe occurs only once in the reference.
 *   Then that single occurrence itself is stored and can be retrieved by the
 *   as_ref_randstrobe() method.
 * - An indirect entry is used if the randstrobe has multiple occurrences.
 *   In that case, offset() and count() point to an interval within a second
 *   table (RandstrobeVector).
 */
class RandstrobeMapEntry {
public:
    RandstrobeMapEntry() { }
    RandstrobeMapEntry(unsigned int offset, unsigned int count) : m_offset(offset), m_count(count) { }

    unsigned int count() const {
        if (is_direct()) {
            return 1;
        } else {
            return m_count;
        }
    }

    unsigned int offset() const{
        assert(!is_direct());
        return m_offset;
    }

    bool is_direct() const {
        return m_count & 0x8000'0000;
    }

    RefRandstrobe as_ref_randstrobe() const {  // 将m_count 首位置为0
        assert(is_direct());
        return RefRandstrobe{m_offset, m_count & 0x7fff'ffff};
    }

    void set_count(unsigned int count) {
        m_count = count;
    }

    void set_offset(unsigned int offset) {
        assert(!is_direct());
        m_offset = offset;
    }

private:
    unsigned int m_offset;
    unsigned int m_count;
};

using RandstrobeMap = robin_hood::unordered_map<randstrobe_hash_t, RandstrobeMapEntry>;

typedef std::vector<uint64_t> hash_vector; //only used during index generation

struct IndexCreationStatistics {
    unsigned int flat_vector_size = 0;
    unsigned int tot_strobemer_count = 0;
    unsigned int tot_occur_once = 0;
    float frac_unique = 0;
    unsigned int tot_high_ab = 0;
    unsigned int tot_mid_ab = 0;
    unsigned int tot_distinct_strobemer_count = 0;
    unsigned int index_cutoff = 0;
    unsigned int filter_cutoff = 0;
    uint64_t unique_mers = 0;

    std::chrono::duration<double> elapsed_hash_index;
    std::chrono::duration<double> elapsed_generating_seeds;
    std::chrono::duration<double> elapsed_unique_hashes;
    std::chrono::duration<double> elapsed_sorting_seeds;
};

struct StrobemerIndex {
    StrobemerIndex(const References& references, const IndexParameters& parameters)
        : filter_cutoff(parameters.filter_cutoff)
        , parameters(parameters)
        , references(references) {}
    unsigned int filter_cutoff = parameters.filter_cutoff; //This also exists in mapping_params
    RefRandstrobeVector flat_vector;
    mutable IndexCreationStatistics stats;

    void write(const std::string& filename) const;
    void read(const std::string& filename);
    void populate(int filter_cutoff, size_t n_threads);
    void print_diagnostics(const std::string& logfile_name, int k) const;
    unsigned int find(uint64_t key) const;
    static const unsigned int N = 28;  // store N bits in the 
    static const uint64_t hash_mask = (((uint64_t)1 << (64 - N)) - 1);
    // RandstrobeMap::const_iterator find(uint64_t key) const {
    //     return randstrobe_map.find(key);
    // }

    uint64_t get_hash(unsigned int position) const {
        if (position < randstrobes_vector.size()){
            return randstrobes_vector[position].hash;
        }else{
            return -1;
        }
    } 

    unsigned int get_strob1_position(unsigned int position) const {
        return randstrobes_vector[position].position;
    }

    int strobe2_offset(unsigned int position) const {
        return randstrobes_vector[position].packed & mask;
    }

    int reference_index(unsigned int position) const {
        return randstrobes_vector[position].packed >> bit_alloc;
    }

    unsigned int get_count(unsigned int position) const {
        unsigned int count = randstrobes_vector[position].hash >> (64 - N);
        return count;
    }

    RandstrobeMap::const_iterator end() const {
        return randstrobe_map.cend();
    }

    void add_entry(uint64_t key, unsigned int offset, unsigned int count) {
        randstrobe_map[key] = RandstrobeMapEntry{offset, count};
    }

    int k() const {
        return parameters.k;
    }

private:
    // std::vector<RefRandstrobeWithHash> add_randstrobes_to_hash_table();
    void add_randstrobes_to_vector(int randstrobe_hashes);
    const IndexParameters& parameters;
    const References& references;
    std::vector<RefRandstrobeWithHash> randstrobes_vector;
    unsigned int* hash_positions = new unsigned int[1 << N](); // the position array used to store the position of the hash in the hash vector;
    RandstrobeMap randstrobe_map; // k-mer -> (offset in flat_vector, occurence count )
    static const int bit_alloc = 8;
    static const int mask = (1 << bit_alloc) - 1;
};

#endif

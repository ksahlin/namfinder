//
// Created by Kristoffer Sahlin on 3/22/22.
//

// Using initial base format of Buffer classed from: https://andrew128.github.io/ProducerConsumer/

#include "pc.hpp"
#include <mutex>
#include <iostream>
#include <chrono>
#include <queue>

#include "timer.hpp"
#include "robin_hood.h"
#include "index.hpp"
#include "kseq++.hpp"

// distribute_interleaved implements the 'interleaved' format:
// If two consequent reads have the same name, they are considered to be a pair.
// Otherwise, they are considered to be single-end reads.

void distribute_interleaved(
    std::vector<klibpp::KSeq>& records,
    std::vector<klibpp::KSeq>& records1,
    std::vector<klibpp::KSeq>& records2,
    std::vector<klibpp::KSeq>& records3,
    std::optional<klibpp::KSeq>& lookahead1
) {
    auto it = records.begin();
    if (lookahead1) {
        if (it != records.end() && lookahead1->name == it->name) {
            records1.push_back(*lookahead1);
            records2.push_back(*it);
            ++it;
        } else {
            records3.push_back(*lookahead1);
        }
        lookahead1 = std::nullopt;
    }
    for (; it != records.end(); ++it) {
        if (it + 1 != records.end() && it->name == (it + 1)->name) {
            records1.push_back(*it);
            records2.push_back(*(it + 1));
            ++it;
        } else {
            records3.push_back(*it);
        }
    }
    if (it != records.end()) {
        lookahead1 = *it;
    }
}


size_t InputBuffer::read_records(std::vector<klibpp::KSeq> &records3,
        int to_read) {
    Timer timer;
//    records1.clear();
//    records2.clear();
    records3.clear();
    // Acquire a unique lock on the mutex
    std::unique_lock<std::mutex> unique_lock(mtx);
    if (to_read == -1) {
        to_read = chunk_size;
    }

    records3 = ks1->stream().read(to_read);

    size_t current_chunk_index = chunk_index;
    chunk_index++;

    if (records3.empty()) {
        finished_reading = true;
    }

    unique_lock.unlock();
    return current_chunk_index;
}

void InputBuffer::rewind_reset() {
    std::unique_lock<std::mutex> unique_lock(mtx);
    ks1->rewind();
    if (ks2) {
        ks2->rewind();
    }
    finished_reading = false;
    chunk_index = 0;
}

void OutputBuffer::output_records(std::string chunk, size_t chunk_index) {
    std::unique_lock<std::mutex> unique_lock(mtx);

    // Ensure we print the chunks in the order in which they were read
    assert(chunks.count(chunk_index) == 0);
    chunks.emplace(std::make_pair(chunk_index, chunk));
    while (true) {
        const auto& item = chunks.find(next_chunk_index);
        if (item == chunks.end()) {
            break;
        }
        out << item->second;
        chunks.erase(item);
        next_chunk_index++;
    }
    unique_lock.unlock();
}


void perform_task(
    InputBuffer &input_buffer,
    OutputBuffer &output_buffer,
    AlignmentStatistics& statistics,
    int& done,
//    const alignment_params &aln_params,
    const mapping_params &map_param,
    const IndexParameters& index_parameters,
    const References& references,
    const StrobemerIndex& index
) {
    bool eof = false;
//    Aligner aligner{aln_params};
    int temp_index = 0;
    while (!eof) {
        temp_index += 1;
//        std::vector<klibpp::KSeq> records1;
//        std::vector<klibpp::KSeq> records2;
        std::vector<klibpp::KSeq> records3;
        auto chunk_index = input_buffer.read_records(records3);
//        assert(records1.size() == records2.size());
//        i_dist_est isize_est;
        if ( records3.empty()
                && input_buffer.finished_reading){
            break;
        }

        std::string nam_out;
        nam_out.reserve(100 * (2* records3.size()));
        for (size_t i = 0; i < records3.size(); ++i) {
            auto record = records3[i];
            align_SE_read(record, nam_out, statistics, map_param, index_parameters, references, index);
        }
        output_buffer.output_records(std::move(nam_out), chunk_index);
        assert(nam_out == "");
    }
    done = true;
}

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <thread>
#include <cassert>
#include <iomanip>
#include <chrono>
#include <unistd.h>

#include "refs.hpp"
#include "exceptions.hpp"
#include "cmdline.hpp"
#include "index.hpp"
#include "pc.hpp"
// #include "aln.hpp"
#include "logger.hpp"
#include "timer.hpp"
#include "version.hpp"

static Logger& logger = Logger::get();


/*
 * Return formatted SAM header as a string
 */
std::string sam_header(const References& references, const std::string& read_group_id, const std::vector<std::string>& read_group_fields) {
    std::stringstream out;
    out << "@HD\tVN:1.6\tSO:unsorted\n";
    for (size_t i = 0; i < references.size(); ++i) {
        out << "@SQ\tSN:" << references.names[i] << "\tLN:" << references.lengths[i] << "\n";
    }
    if (!read_group_id.empty()) {
        out << "@RG\tID:" << read_group_id;
        for (const auto& field : read_group_fields) {
           out << '\t' << field;
        }
        out << '\n';
    }
    out << "@PG\tID:strobealign\tPN:strobealign\tVN:" << version_string() << "\tCL:strobealign\n";
    return out.str();
}

void warn_if_no_optimizations() {
    if (std::string(CMAKE_BUILD_TYPE) == "Debug") {
        logger.info() << "\n    ***** Binary was compiled without optimizations - this will be very slow *****\n\n";
    }
}

void log_parameters(const IndexParameters& index_parameters, const mapping_params& map_param) {
    logger.debug() << "Using" << std::endl
        << "k: " << index_parameters.k << std::endl
        << "s: " << index_parameters.s << std::endl
        << "w_min: " << index_parameters.w_min << std::endl
        << "w_max: " << index_parameters.w_max << std::endl
        << "Read length (r): " << map_param.r << std::endl
        << "Maximum seed length: " << index_parameters.max_dist + index_parameters.k << std::endl
        << "C: " << map_param.C << std::endl
        << "L: " << map_param.L << std::endl
        << "Expected [w_min, w_max] in #syncmers: [" << index_parameters.w_min << ", " << index_parameters.w_max << "]" << std::endl
        << "Expected [w_min, w_max] in #nucleotides: [" << (index_parameters.k - index_parameters.s + 1) * index_parameters.w_min << ", " << (index_parameters.k - index_parameters.s + 1) * index_parameters.w_max << "]" << std::endl;
}

bool avx2_enabled() {
#ifdef __AVX2__
    return true;
#else
    return false;
#endif
}

InputBuffer get_input_buffer(const CommandLineOptions& opt) {
        return InputBuffer(opt.reads_filename1, opt.chunk_size);
}


int run_strobealign(int argc, char **argv) {
    auto opt = parse_command_line_arguments(argc, argv);

//    logger.set_level(opt.verbose ? LOG_DEBUG : LOG_INFO);
    logger.info() << std::setprecision(2) << std::fixed;
    logger.info() << "This is namfinder " << version_string() << '\n';
    logger.debug() << "Build type: " << CMAKE_BUILD_TYPE << '\n';
    warn_if_no_optimizations();
    logger.debug() << "AVX2 enabled: " << (avx2_enabled() ? "yes" : "no") << '\n';

//    if (opt.c >= 64 || opt.c <= 0) {
//        throw BadParameter("c must be greater than 0 and less than 64");
//    }

    InputBuffer input_buffer = get_input_buffer(opt);
    input_buffer.rewind_reset();
    IndexParameters index_parameters = IndexParameters(
            opt.k, opt.s,opt.l,opt.u, opt.max_seed_len, opt.filter_cutoff  );
    logger.debug() << index_parameters << '\n';
//    alignment_params aln_params;
//    aln_params.match = opt.A;
//    aln_params.mismatch = opt.B;
//    aln_params.gap_open = opt.O;
//    aln_params.gap_extend = opt.E;
//    aln_params.end_bonus = opt.end_bonus;

    mapping_params map_param;
    map_param.r = opt.r;
//    map_param.max_secondary = opt.max_secondary;
//    map_param.dropoff_threshold = opt.dropoff_threshold;
//    map_param.R = opt.R;
//    map_param.maxTries = opt.maxTries;
//    map_param.is_sam_out = opt.is_sam_out;
//    map_param.output_unmapped = opt.output_unmapped;
    map_param.filter_cutoff = opt.filter_cutoff;

    log_parameters(index_parameters, map_param);
    logger.debug() << "Threads: " << opt.n_threads << std::endl;

//    assert(k <= (w/2)*w_min && "k should be smaller than (w/2)*w_min to avoid creating short strobemers");

    // Create index
    References references;
    Timer read_refs_timer;
    references = References::from_fasta(opt.ref_filename);
    logger.info() << "Time reading reference: " << read_refs_timer.elapsed() << " s\n";

    logger.info() << "Reference size: " << references.total_length() / 1E6 << " Mbp ("
        << references.size() << " contig" << (references.size() == 1 ? "" : "s")
        << "; largest: "
        << (*std::max_element(references.lengths.begin(), references.lengths.end()) / 1E6) << " Mbp)\n";
    if (references.total_length() == 0) {
        throw InvalidFasta("No reference sequences found");
    }

    StrobemerIndex index(references, index_parameters);
        logger.info() << "Indexing ...\n";
        Timer index_timer;
        logger.debug() << "FILTER CUTOFF: " << std::to_string(opt.filter_cutoff) << std::endl;

        index.populate(opt.filter_cutoff, opt.n_threads);
        
        logger.info() << "  Time generating seeds: " << index.stats.elapsed_generating_seeds.count() << " s" <<  std::endl;
        logger.info() << "  Time estimating number of unique hashes: " << index.stats.elapsed_unique_hashes.count() << " s" <<  std::endl;
        logger.info() << "  Time sorting non-unique seeds: " << index.stats.elapsed_sorting_seeds.count() << " s" <<  std::endl;
        logger.info() << "  Time generating hash table index: " << index.stats.elapsed_hash_index.count() << " s" <<  std::endl;
        logger.info() << "Total time indexing: " << index_timer.elapsed() << " s\n";

        logger.debug()
        << "Unique strobemers: " << index.stats.unique_mers << std::endl
        << "Total strobemers count: " << index.stats.tot_strobemer_count << std::endl
        << "Total strobemers occur once: " << index.stats.tot_occur_once << std::endl
        << "Fraction Unique: " << index.stats.frac_unique << std::endl
        << "Total strobemers highly abundant > 100: " << index.stats.tot_high_ab << std::endl
        << "Total strobemers mid abundance (between 2-100): " << index.stats.tot_mid_ab << std::endl
        << "Total distinct strobemers stored: " << index.stats.tot_distinct_strobemer_count << std::endl;
        if (index.stats.tot_high_ab >= 1) {
            logger.debug() << "Ratio distinct to highly abundant: " << index.stats.tot_distinct_strobemer_count / index.stats.tot_high_ab << std::endl;
        }
        if (index.stats.tot_mid_ab >= 1) {
            logger.debug() << "Ratio distinct to non distinct: " << index.stats.tot_distinct_strobemer_count / (index.stats.tot_high_ab + index.stats.tot_mid_ab) << std::endl;
        }
        logger.debug() << "Filtered cutoff index: " << index.stats.index_cutoff << std::endl;
        logger.debug() << "Filtered cutoff count: " << index.stats.filter_cutoff << std::endl;
        
        if (!opt.logfile_name.empty()) {
            index.print_diagnostics(opt.logfile_name, index_parameters.k);
            logger.debug() << "Finished printing log stats" << std::endl;
        }


    // Map/align reads
        
    Timer map_align_timer;
//    map_param.rescue_cutoff = map_param.R < 100 ? map_param.R * index.filter_cutoff : 1000;
//    logger.debug() << "Using rescue cutoff: " << map_param.rescue_cutoff << std::endl;

    std::streambuf* buf;
    std::ofstream of;

    if (!opt.write_to_stdout) {
        of.open(opt.output_file_name);
        buf = of.rdbuf();
    }
    else {
        buf = std::cout.rdbuf();
    }

    std::ostream out(buf);

    std::vector<AlignmentStatistics> log_stats_vec(opt.n_threads);

    logger.info() << "Running in " << (opt.is_SE ? "single-end" : "paired-end") << " mode" << std::endl;

    OutputBuffer output_buffer(out);

    std::vector<std::thread> workers;
    std::vector<int> worker_done(opt.n_threads);  // each thread sets its entry to 1 when it’s done
    for (int i = 0; i < opt.n_threads; ++i) {
        std::thread consumer(perform_task, std::ref(input_buffer), std::ref(output_buffer),
            std::ref(log_stats_vec[i]), std::ref(worker_done[i]),
            std::ref(map_param), std::ref(index_parameters), std::ref(references),
            std::ref(index), std::ref(opt.read_group_id));
        workers.push_back(std::move(consumer));
    }

    for (auto& worker : workers) {
        worker.join();
    }
    logger.info() << "Done!\n";

    AlignmentStatistics tot_statistics;
    for (auto& it : log_stats_vec) {
        tot_statistics += it;
    }

    logger.info() << "Total mapping sites tried: " << tot_statistics.tot_all_tried << std::endl
        << "Total calls to ssw: " << tot_statistics.tot_aligner_calls << std::endl
        << "Calls to ksw (rescue mode): " << tot_statistics.tot_rescued << std::endl
        << "Did not fit strobe start site: " << tot_statistics.did_not_fit << std::endl
        << "Tried rescue: " << tot_statistics.tried_rescue << std::endl
        << "Total time mapping: " << map_align_timer.elapsed() << " s." << std::endl
        << "Total time reading read-file(s): " << tot_statistics.tot_read_file.count() / opt.n_threads << " s." << std::endl
        << "Total time creating strobemers: " << tot_statistics.tot_construct_strobemers.count() / opt.n_threads << " s." << std::endl
        << "Total time finding NAMs (non-rescue mode): " << tot_statistics.tot_find_nams.count() / opt.n_threads << " s." << std::endl
        << "Total time finding NAMs (rescue mode): " << tot_statistics.tot_time_rescue.count() / opt.n_threads << " s." << std::endl;
    //<< "Total time finding NAMs ALTERNATIVE (candidate sites): " << tot_find_nams_alt.count()/opt.n_threads  << " s." <<  std::endl;
    logger.info() << "Total time sorting NAMs (candidate sites): " << tot_statistics.tot_sort_nams.count() / opt.n_threads << " s." << std::endl
        << "Total time base level alignment (ssw): " << tot_statistics.tot_extend.count() / opt.n_threads << " s." << std::endl
        << "Total time writing alignment to files: " << tot_statistics.tot_write_file.count() << " s." << std::endl;
    return EXIT_SUCCESS;
}

int main(int argc, char **argv) {
    try {
        return run_strobealign(argc, argv);
    } catch (BadParameter& e) {
        logger.error() << "A mapping or seeding parameter is invalid: " << e.what() << std::endl;
    } catch (const std::runtime_error& e) {
        logger.error() << "strobealign: " << e.what() << std::endl;
    }
    return EXIT_FAILURE;
}

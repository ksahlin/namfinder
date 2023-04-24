#include "nam.hpp"
#include "logger.hpp"

static Logger& logger = Logger::get();
namespace {

struct Hit {
    int query_s;
    int query_e;
    int ref_s;
    int ref_e;
    bool is_rc = false;
};

void add_to_hits_per_ref(
    robin_hood::unordered_map<unsigned int, std::vector<Hit>>& hits_per_ref,
    int query_s,
    int query_e,
    bool is_rc,
    const StrobemerIndex& index,
    // RandstrobeMapEntry randstrobe_map_entry,
    unsigned int position,
    int min_diff,
    int& tot_hits
) {
    // Determine whether the hash table’s value directly represents a
    // ReferenceMer (this is the case if count==1) or an offset/count
    // pair that refers to entries in the flat_vector.
    unsigned int count = index.get_count(position);
    if (count == 1) {
        // auto r = randstrobe_map_entry.as_ref_randstrobe();
        int ref_s = index.get_strob1_position(position);
        int ref_e = ref_s + index.strobe2_offset(position) + index.k();
        int diff = std::abs((query_e - query_s) - (ref_e - ref_s));
        if (diff <= min_diff) {
            hits_per_ref[index.reference_index(position)].push_back(Hit{query_s, query_e, ref_s, ref_e, is_rc});
            min_diff = diff;
            tot_hits ++;
        }
    } else {
//        logger.debug() << "add_to_hits_per_ref count: " << std::to_string(count) << std::endl;
        for (unsigned int j = position; j < position + count; ++j) {
            int ref_s = index.get_strob1_position(j);
            int ref_e = ref_s + index.strobe2_offset(j) + index.k();
            int diff = std::abs((query_e - query_s) - (ref_e - ref_s));
            if (diff <= min_diff) {
                hits_per_ref[index.reference_index(j)].push_back(Hit{query_s, query_e, ref_s, ref_e, is_rc});
                min_diff = diff;
                tot_hits ++;
            }
        }
    }
}

std::vector<Nam> merge_hits_into_nams(
    robin_hood::unordered_map<unsigned int, std::vector<Hit>>& hits_per_ref,
    int k,
    bool sort
) {
    std::vector<Nam> nams;
    int nam_id_cnt = 0;
    for (auto &[ref_id, hits] : hits_per_ref) {
        if (sort) {
            std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) -> bool {
                    // first sort on query starts, then on reference starts
                    return (a.query_s < b.query_s) || ( (a.query_s == b.query_s) && (a.ref_s < b.ref_s) );
                }
            );
        }

        std::vector<Nam> open_nams;
        unsigned int prev_q_start = 0;
        for (auto &h : hits) {
            bool is_added = false;
            for (auto & o : open_nams) {

//                // Extend NAM
//                if (( o.is_rc == h.is_rc) && (o.query_prev_hit_startpos < h.query_s) && (h.query_s <= o.query_e ) && (o.ref_prev_hit_startpos < h.ref_s) && (h.ref_s <= o.ref_e) ){
//                    if ( (h.query_e > o.query_e) && (h.ref_e > o.ref_e) ) {
//                        o.query_e = h.query_e;
//                        o.ref_e = h.ref_e;
//                        o.query_prev_hit_startpos = h.query_s; // ` the last strobemer hit in case of outputting paf
//                        o.ref_prev_hit_startpos = h.ref_s; // log the last strobemer hit in case of outputting paf
//                        o.n_hits ++;
//                        is_added = true;
//                        break;
//                    }
//                    else if ((h.query_e <= o.query_e) && (h.ref_e <= o.ref_e)) {
//                        o.query_prev_hit_startpos = h.query_s; // log the last strobemer hit in case of outputting paf
//                        o.ref_prev_hit_startpos = h.ref_s; // log the last strobemer hit in case of outputting paf
//                        o.n_hits ++;
//                        is_added = true;
//                        break;
//                    }
//                }

                // Extend NAM (version from StrobeMap)
                if ( ( o.query_prev_hit_startpos < h.query_s) && (h.query_s <= o.query_e ) && ( o.ref_prev_hit_startpos <= h.ref_s) && (h.ref_s <= o.ref_e) ){

                    if (h.query_e > o.query_e) {
                        o.query_e = h.query_e;
                    }
                    if (h.ref_e > o.ref_e) {
                        o.ref_e = h.ref_e;
                    }
                    o.query_prev_hit_startpos = h.query_s;
                    o.ref_prev_hit_startpos = h.ref_s; // keeping track so that we don't . Can be caused by interleaved repeats.
                    o.n_hits ++;
                    is_added = true;
                    break;
                }

            }
            // Add the hit to open matches
            if (!is_added){
                Nam n;
                n.nam_id = nam_id_cnt;
                nam_id_cnt ++;
                n.query_s = h.query_s;
                n.query_e = h.query_e;
                n.ref_s = h.ref_s;
                n.ref_e = h.ref_e;
                n.ref_id = ref_id;
//                n.previous_query_start = h.query_s;
//                n.previous_ref_start = h.ref_s;
                n.query_prev_hit_startpos = h.query_s;
                n.ref_prev_hit_startpos = h.ref_s;
                n.n_hits = 1;
                n.is_rc = h.is_rc;
//                n.score += (float)1 / (float)h.count;
                open_nams.push_back(n);
            }

            // Only filter if we have advanced at least k nucleotides
            if (h.query_s > prev_q_start + k) {

                // Output all NAMs from open_matches to final_nams that the current hit have passed
                for (auto &n : open_nams) {
                    if (n.query_e < h.query_s) {
                        int n_max_span = std::max(n.query_span(), n.ref_span());
                        int n_min_span = std::min(n.query_span(), n.ref_span());
                        float n_score;
                        n_score = ( 2*n_min_span -  n_max_span) > 0 ? (float) (n.n_hits * ( 2*n_min_span -  n_max_span) ) : 1;   // this is really just n_hits * ( min_span - (offset_in_span) ) );
//                        n_score = n.n_hits * n.query_span();
                        n.score = n_score;
                        nams.push_back(n);
                    }
                }

                // Remove all NAMs from open_matches that the current hit have passed
                auto c = h.query_s;
                auto predicate = [c](decltype(open_nams)::value_type const &nam) { return nam.query_e < c; };
                open_nams.erase(std::remove_if(open_nams.begin(), open_nams.end(), predicate), open_nams.end());
                prev_q_start = h.query_s;
            }
        }

        // Add all current open_matches to final NAMs
        for (auto &n : open_nams) {
            int n_max_span = std::max(n.query_span(), n.ref_span());
            int n_min_span = std::min(n.query_span(), n.ref_span());
            float n_score;
            n_score = ( 2*n_min_span -  n_max_span) > 0 ? (float) (n.n_hits * ( 2*n_min_span -  n_max_span) ) : 1;   // this is really just n_hits * ( min_span - (offset_in_span) ) );
//            n_score = n.n_hits * n.query_span();
            n.score = n_score;
            nams.push_back(n);
        }
    }
    logger.debug() << "NAMS: " << std::to_string(std::size(nams)) << std::endl;
    return nams;
}

} // namespace

/*
 * Find a query’s NAMs, ignoring randstrobes that occur too often in the
 * reference (have a count above filter_cutoff).
 *
 * Return the fraction of nonrepetitive hits (those not above the filter_cutoff threshold)
 */
std::pair<float, std::vector<Nam>> find_nams(
    const QueryRandstrobeVector &query_randstrobes,
    const StrobemerIndex& index
) {
    robin_hood::unordered_map<unsigned int, std::vector<Hit>> hits_per_ref;
    hits_per_ref.reserve(100);

    /*
    1. Find the hash in the vector
    2. the occur times of the hash value, use a flag 
    3. need to know reference index, strobe1 position, storbe2 - strobe1
    */
    int nr_good_hits = 0, total_hits = 0, tot_hits = 0;
    for (const auto &q : query_randstrobes) {
        unsigned int position = index.find(q.hash);
        if (position != -1){
            total_hits++;
            unsigned int count = index.get_count(position);
//            logger.debug() << "COUNT: " << std::to_string(count)  << "FILTER CUTOFF: " << std::to_string(index.filter_cutoff) << std::endl;
            if (count > index.filter_cutoff){
                continue;
            } 
            nr_good_hits++;
            add_to_hits_per_ref(hits_per_ref, q.start, q.end, q.is_reverse, index, position, 100'000, tot_hits);
        }
    }
    logger.debug() << "add_to_hits_per_ref DONE: " << std::to_string(hits_per_ref.size()) << std::endl;
    logger.debug() << "add_to_hits_per_ref TOT count: " << std::to_string(tot_hits) << std::endl;

    float nonrepetitive_fraction = total_hits > 0 ? ((float) nr_good_hits) / ((float) total_hits) : 1.0;
    auto nams = merge_hits_into_nams(hits_per_ref, index.k(), false);
    logger.debug() << "merge_hits_into_nams DONE: " << std::to_string(nams.size()) << std::endl;

    return make_pair(nonrepetitive_fraction, nams);
}


std::ostream& operator<<(std::ostream& os, const Nam& n) {
    os << "Nam(query: " << n.query_s << ".." << n.query_e << ", ref: " << n.ref_s << ".." << n.ref_e << ", score=" << n.score << ")";
    return os;
}

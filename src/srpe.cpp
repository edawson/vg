#include "srpe.hpp"

using namespace std;
namespace vg{

    bool CompareBreakpoint(const BREAKPOINT& a, const BREAKPOINT& b){
        id_t a_node = std::get<0>(a.position);
        off_t a_off = std::get<2>(a.position);

        id_t b_node = std::get<0>(b.position);
        off_t b_off = std::get<2>(b.position);

        if (a_node == b_node){
            return a_off < b_off;
        }
        return a_node < b_node;
    };

    void SRPE::merge_breakpoints(vector<BREAKPOINT>& bps, int dist){
        std::sort(bps.begin(), bps.end(), CompareBreakpoint);

        vector<BREAKPOINT> merged;

        if (bps.size() > 2){
            merged.push_back(bps[0]);
            for (int i = 1; i < bps.size(); ++i){
                if (overlap(merged.back(), bps[i], dist)){
                    merged.back().merge(bps[i]);
                }
                else{
                    merged.push_back(bps[i]);
                }
            }

        }

        
        
    }

    bool SRPE::overlap(BREAKPOINT a, BREAKPOINT p, int dist){

            if (a.start > -1 && p.start > -1){
                if ( abs(a.start - p.start) < dist){
                    return true;
                }
            }
            else{
                // if (std::get<0>(a.position) == std::get<0>(p.position) && abs(std::get<2>(a.position) - std::get<2>(p.position)) < dist){
                //     return true;
                // }
                // else{
                    LRUCache<id_t,Node> ncache (1000);
                    LRUCache<id_t,vector<Edge>> ecache(1000);
                    int xgdist = xg_cached_distance(a.position, p.position, 100000, ff.my_xg_index, ncache, ecache);
                    return xgdist < dist;
                // }
            }
            return false;
        }

    double SRPE::discordance_score(vector<pair<Alignment&, Alignment&>> alns, VG* subgraph){
    // Sum up the mapping scores
    // subtract the soft clips
    // the discordant insert size reads
    // and a flat penalty for discordant orientation

    double ret = 0.0;
    for (int i = 0; i < alns.size(); ++i){

    }

    return ret;
    }

    void SRPE::call_svs_paired_end(vg::VG* graph, istream& gamstream, vector<BREAKPOINT>& bps, string refpath){

        vector<pair<Alignment, Alignment>> oeas_to_assemble;
        vector<pair<Alignment, Alignment>> unmapped_to_assemble;

        std::function<void(Alignment&, Alignment&)> pe_func = [&](Alignment& a, Alignment& b){
            bool supports_sv = ff.mark_sv_alignments(a, b);
            if (supports_sv){
                BREAKPOINT bpoint;
                BREAKPOINT secondary_break;
                bool break_is_good = true;
                if (!a.read_mapped() && a.mate_unmapped()){
                    // Unmapped
                    break_is_good = false;
                }
                else if (a.read_mapped() == a.mate_unmapped()){
                    // One end anchored
                    // possible insertion or deletion anchor
                }
                if (ff.pair_orientation_filter(a, b)){
                    // Inversion hint
                    if (a.read_on_reverse_strand()){
                        Path p = trim_hanging_ends(a.path());
                        Mapping rmap = p.mapping(p.mapping_size() - 1);
                        Position pos = rmap.position();
                        bpoint.position = make_pos_t(pos);
                        bpoint.add_pe_support(a, b);

                        Path p2 = trim_hanging_ends(b.path());
                        Position pos2 = p2.mapping(0).position();
                        secondary_break.position = make_pos_t(pos2);
                        bpoint.mates.push_back(secondary_break);



                    }
                    else{
                        Path p = trim_hanging_ends(a.path());
                        Mapping rmap = p.mapping(0);
                        Position pos = rmap.position();
                        bpoint.position = make_pos_t(pos);
                        bpoint.add_pe_support(a, b);
                        bpoint.sv_type_counts[3]++;

                        p = trim_hanging_ends(b.path());
                        pos = p.mapping(p.mapping_size() - 1).position();
                        secondary_break.position = make_pos_t(pos);
                        bpoint.mates.push_back(secondary_break);
                    }
                }

                if (a.discordant_insert_size()){
                    // Insertion or deletion
                    Path p = trim_hanging_ends(a.path());
                    Path p2 = trim_hanging_ends(b.path());

                    Position pos = p.mapping(p.mapping_size() - 1).position();
                    bpoint.position = make_pos_t(pos);

                    Position pos2 = p2.mapping(0).position();
                    secondary_break.position = make_pos_t(pos2);
                    bpoint.mates.push_back(secondary_break);
                }

                if (a.soft_clipped() | b.soft_clipped()){

                }

                if (break_is_good){
                    #pragma omp critical
                    bps.push_back(bpoint);
                }
                

            }
        };

        stream::for_each_interleaved_pair_parallel(gamstream, pe_func);

        

        merge_breakpoints(bps, 100);

        cerr << bps.size() << " breakpoints found and merged." << endl;



    }

    void SRPE::call_svs_split_read(vg::VG* graph, istream& gamstream, vector<BREAKPOINT>& bps, string refpath){
        
        std::function<void(Alignment&)> sr_func = [&](Alignment& a){

        };

        stream::for_each_parallel(gamstream, sr_func);

    }


    void SRPE::call_svs(string graphfile, string gamfile, string refpath){
        vg::VG* graph;
        if (!graphfile.empty()){
            ifstream in(graphfile);
            graph = new VG(in, false);
        }
        ifstream gamstream;
        gamstream.open(gamfile);
        // Set up path index
        ff.set_my_vg(graph);
        ff.soft_clip_limit = 20;
        ff.fill_node_to_position(refpath);
    
    //     std::function<vector<BREAKPOINT> (vector<BREAKPOINT>)> merge_breakpoints = [](vector<BREAKPOINT> bps){
    //     vector<BREAKPOINT> ret;
    //     BREAKPOINT sent;
    //     sent.start = -100;
    //     ret.push_back(sent);
    //     for (int i = 0; i < bps.size(); i++){
    //         BREAKPOINT a = bps[i];
    //         bool merged = false;
    //         for (int j = 0; j < ret.size(); j++){
    //             if (ret[j].overlap(a, 20)){
    //                 ret[j].other_supports += 1;
    //                 merged = true;
    //             }
    //         }
    //         if (!merged){
    //             ret.push_back(a);
    //         }
    //     }
    //     return ret;
    // };

    vector<BREAKPOINT> pe_bps;
    vector<BREAKPOINT> sr_bps;

    

    }

    void SRPE::aln_to_bseq(Alignment& a, bseq1_t* read){
        read->seq = (char*) a.sequence().c_str();
        read->qual = (char*) a.quality().c_str();
        read->l_seq = a.sequence().length();
    }

    /*
        typedef struct {
	    int32_t len;      // length of sequence
	    int32_t nsr;      // number of supporting reads
	    char *seq;        // unitig sequence
	    char *cov;        // cov[i]-33 gives per-base coverage at i
	    int n_ovlp[2];    // number of 5'-end [0] and 3'-end [1] overlaps
	    fml_ovlp_t *ovlp; // overlaps, of size n_ovlp[0]+n_ovlp[1]
        } fml_utg_t;
    */

    /**
    * function assemble
    * inputs: a vector of Alignments to be assembled (based on their sequences)
    * outputs: 
    */
    void SRPE::assemble(vector<Alignment> alns, vector<fml_utg_t>& unitigs){
        int n_seqs, n_utgs;
        n_seqs = alns.size();
        bseq1_t* mr_bseqs = new bseq1_t [alns.size()];
        for (int i = 0; i < n_seqs; ++i){
            aln_to_bseq( alns[i], mr_bseqs + i );
        }
        fml_utg_t *utgs;
        fml_opt_t opt;
        fml_opt_init(&opt);
        utgs = fml_assemble(&opt, n_seqs, mr_bseqs, &n_utgs);
        for (int i = 0; i < n_utgs; ++i){
            unitigs.push_back( *(utgs + i) );
        }

        fml_utg_destroy(n_utgs, utgs);
    }

    void SRPE::assemble(string refpath, int64_t start_pos, int64_t end_pos, vector<fml_utg_t>& unitigs){
        // Get all alignments on <refpath> from <startpos> to <endpos>
    }
    void SRPE::assemble(int64_t node_id, int64_t pos, int window_size){

    }



}


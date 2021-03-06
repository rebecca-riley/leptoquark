#include <fstream>
#include <iterator>
#include <chrono>
#include "fastjet/ClusterSequence.hh"
using namespace std;
using namespace fastjet;

// --------- GLOBAL FLAGS/CONSTANTS --------- //
// true = display event processing messages on terminal
const bool WRITE_TO_TERM = true;
// true = write runlog to file
const bool WRITE_TO_FILE = false;
ofstream runlog;
// true = don't print failed events to terminal or output file
const bool SUPPRESS_FAILURE_OUTPUT = true;
// true = keep count of how many final state particles/neutrinos, check for all
//        particles clustered, prints all failure messages regardless of suppression
const bool DEBUG_MODE = false;
// indices
const int PX = 3, PY = 4, PZ = 5, E = 6;
// colors
const int RED = 31, GREEN = 32, YELLOW = 33, BLUE = 34, PINK = 35, CYAN = 36, GREY = 37;

// --------- HELPER FUNCTION DECLARATIONS --------- //
PseudoJet get_jet(vector<string> delimited);
void print_jet(PseudoJet jet, string identifier = "");
void _print_event(int event_number, string message, int color_message = GREY,
                  string other_info = "", int color_other_info = GREY);
void print_success(int event_number, string message, string other_info = "",
                   int color_message = GREEN, int color_other_info = GREY);
void print_error(int event_number, string message, int color = RED);
void print_warning(int event_number, string message, int color = PINK);
vector<string> split_line(string line);
double get_spatial_separation(PseudoJet jet1, PseudoJet jet2);

// --------- CUT LIST --------- //
// (1) select final state with two taus and two b-jets
// (2) T(l): pt > 50 GeV, |eta| < 2.1
// (3) T(h): pt > 50 GeV, |eta| < 2.3
// (4) T(l), T(h) must originate from the same vertex, have opposite charge, and have
//     spatial separation sqrt((del_phi)^2 + (del_eta)^2) > 0.5
// (5) b-jet reconstruction: anti-kt algorithm, R = 0.4
//     jet requirements: pt > 50 GeV, |eta| < 2.4, spatial sep. from T(l),T(h) > 0.5
// (6) at least 1 of 2 leading (highest pt?) jets is "b-tagged"
// (7) select T(l)/T(h) + b_grtr/b_less combination that minimizes mass difference
//     between pair
// (8) check that T(h) + jet > 250 GeV


int main() {
    // --------- RUNTIME --------- //
    auto start = std::chrono::high_resolution_clock::now();

    // --------- CONSTANTS --------- //
    // input, output file names
    const string input_filename = "ditop_experiment4.hepmc",
                 runlog_filename = "runlog.txt";
    // indices
    const int barcode = 1, pdg_code = 2, px = 3, py = 4, pz = 5, E = 6, gen_mass = 7,
              status = 8;
    // particle codes
    const string nu_e = "12", nu_mu = "14", nu_tau = "16", nu_tau_pr = "18";
    const string tau_p = "15", tau_m = "-15", b_p = "5", b_m = "-5";
    // status flags
    const string final = "1";
    // jet clustering parameters
    const double R = 0.4;
    const int jet_pt_cutoff = 50;
    // null jet
    const PseudoJet null_jet = PseudoJet(0,0,0,0);


    /// --------- VARIABLES --------- //
    // global counters
    int events = 0, num_fail = 0;
    int total_final_state = 0, total_neutrinos = 0;
    // event counters
    int final_state = 0, neutrinos = 0, taus = 0, num_particles_clustered = 0;
    // event vertex
    PseudoJet current_vertex;     // bastardization of PseudoJet to use its functions
    // vectors
    vector<PseudoJet> particles;
    PseudoJet tau_candidates[2];
    vector<vector<PseudoJet>> passing_events;
    // helper struct
    struct Tau {
        PseudoJet tau;
        bool is_tau_plus;
        PseudoJet vertex;
    };
    vector<Tau> vec_taus;
    // FastJet setup
    JetDefinition jet_def(antikt_algorithm, R);

    // --------- FILE IO SETUP --------- //
    ifstream hepmc_file;
    hepmc_file.open(input_filename);
    if (WRITE_TO_FILE) runlog.open(runlog_filename);


    // --------- HEPMC PARSING --------- //
    string line;
    getline(hepmc_file,line);
    while(line[0]!='E') {                       // get rid of lines until first event
        getline(hepmc_file,line);
    } // found first event

    while(!hepmc_file.eof()) {                  // process full event
        events += 1;
        // reset event counters
        final_state = 0;
        neutrinos = 0;
        taus = 0;
        num_particles_clustered = 0;
        // reset vectors, arrays
        particles.clear();
        vec_taus.clear();
        tau_candidates[0] = null_jet; tau_candidates[1] = null_jet;

        getline(hepmc_file,line);               // get rid of event header line

        while(!hepmc_file.eof() && line[0]!='E') {  // process lines within event
            if(line[0]=='P'){                       // only process particles
                vector<string> delimited = split_line(line);

                // only want final state particles
                if(delimited[status]==final) {
                    // counters //
                    final_state += 1;
                    if (DEBUG_MODE) total_final_state += 1;

                    // skip storing neutrinos
                    if(delimited[pdg_code] == nu_e || delimited[pdg_code] == nu_mu ||
                       delimited[pdg_code] == nu_tau || delimited[pdg_code] == nu_tau_pr)
                    {
                        // counters //
                        neutrinos += 1;
                        if (DEBUG_MODE) total_neutrinos += 1;
                        goto NextItem;
                    }

                    // skip storing taus
                    // store instead tau,charge,vertex info in Tau vector 'vec_taus'
                    if(delimited[pdg_code] == tau_p || delimited[pdg_code] == tau_m) {
                        vec_taus.push_back( Tau{get_jet(delimited),
                                            (delimited[pdg_code]==tau_p ? true:false),
                                            current_vertex} );
                        taus += 1;
                        goto NextItem;
                    }

                    // store all other final state particles in vector 'particles'
                    // if particle is b, tag before storing
                    if(delimited[pdg_code] == b_p || delimited[pdg_code] == b_m) {
                        PseudoJet b = get_jet(delimited);
                        b.set_user_index(1);
                        particles.push_back(b);
                    }
                    else particles.push_back(get_jet(delimited));
                }
            }
            else if(line[0]=='V') current_vertex.reset(get_jet(split_line(line)));

            NextItem:
            getline(hepmc_file,line);
        } // next event reached or eof


        // --------- CUTS ON TAUS --------- //
        // cut 1 -- at least two taus
        if (vec_taus.size() < 2) {
            num_fail++;
            print_error(events,"failed cut 1: lacks tau+ or tau-");
            continue;                           // continue => process next event
        }

        // cuts 2-4
        bool vertex_match = false, opposite_charge = false, pt_pass = false,
             eta_pass = false, spatially_separated = false;
        for (int i = 0; i < vec_taus.size() - 1; i++) {     // compare all combinations of taus
            for (int j = i + 1; j < vec_taus.size(); j++) {
                // cut 4 -- same vertex
                // check vertex first, other criteria only checked if vertex matches
                if (vec_taus[i].vertex == vec_taus[j].vertex) {
                    // if more than one vertex match in event (highly unusual), print warning
                    if (vertex_match == true)
                        print_warning(events,"ANALYSIS NOT VALID FOR THIS EVENT -- "
                              "more than one tau production vertex, check manually");

                    vertex_match = true;
                    // cut 4 -- opposite charge
                    if (vec_taus[i].is_tau_plus == !vec_taus[j].is_tau_plus)
                        opposite_charge = true;
                    // cut 2 -- pt > 50
                    if (vec_taus[i].tau.pt() > 50 && vec_taus[j].tau.pt() > 50)
                        pt_pass = true;
                    // cut 3 -- |eta| < 2.3
                    if (abs(vec_taus[i].tau.eta()) < 2.3 && abs(vec_taus[j].tau.eta()) < 2.3)
                        eta_pass = true;
                    // cut 4 -- spatial separation
                    if (get_spatial_separation(vec_taus[i].tau,vec_taus[j].tau) > 0.5)
                        spatially_separated = true;

                    // if taus are from same vertex, they go on array as candidates
                    // if there is more than one production vertex per event, this
                    //    will cause problems since only taus from the last vertex
                    //    will go on the array (very unlikely to happen)
                    tau_candidates[0] = vec_taus[i].tau;
                    tau_candidates[1] = vec_taus[j].tau;
                }
            }
        }
        if (!vertex_match) {
            num_fail++;
            print_error(events,"failed cut 4: taus originate from different vertices");
            continue;
        }
        if (!opposite_charge) {
            num_fail++;
            print_error(events,"failed cut 4: taus have same opposite charges");
            continue;
        }
        if (!pt_pass) {
            num_fail++;
            print_error(events,"failed cut 2: taus have pt <= 50");
            continue;
        }
        if (!eta_pass) {
            num_fail++;
            print_error(events,"failed cut 3: taus have |eta| >= 2.3");
            continue;
        }
        if (!spatially_separated) {
            num_fail++;
            print_error(events,"failed cut 4: taus have separation <= 0.5");
            continue;
        }


        // --------- JET CLUSTERING, CUTS ON JETS --------- //
        ClusterSequence cs(particles, jet_def);
        // cut 5 -- pt_jet > 50
        vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(jet_pt_cutoff));

        // cut 5 -- |eta_jet| < 2.4, spatial separation from taus > 0.5
        for (unsigned i = 0; i < jets.size(); i++) {
            if (    abs(jets[i].eta()) >= 2.4 ||
                    get_spatial_separation(jets[i],tau_candidates[0]) <= 0.5 ||
                    get_spatial_separation(jets[i],tau_candidates[1]) <= 0.5    )
                jets.erase(jets.begin()+i);
        }

        // cut 5 -- at least two jets must remain for possible b-tagging
        if (jets.size() < 2) {
            num_fail++;
            print_error(events,"failed cut 5: less than two jets remaining after"
                                "eta, separation cuts");
            continue;
        }

        // cut 6 -- at least two jets contain b's
        int num_jets_with_b = 0;
        for (int i = 0; i < jets.size(); i++) {
            vector<PseudoJet> constits = jets[i].constituents();
            for (int j = 0; j < constits.size(); j++) {
                if (constits[j].user_index() == 1) {
                    num_jets_with_b++;
                    break;
                }
            }
        }

        if (num_jets_with_b < 1) {
            num_fail++;
            print_error(events,"failed cut 6: no jets containing b's");
            continue;
        }


        // --------- JET PAIRING, CUTS ON PAIRS --------- //
        // cut 7 -- select tau, b pairs that minimizes mass difference b/w pairs
        double min_diff = 10000000000000;       // very large number -- will be reset
        struct {                                // immediately by first mass_diff
            PseudoJet jet[2];
            Tau tau[2];
        } min_pairs;

        for (int i = 0; i < jets.size(); i++) {
            for (int j = 0; j < vec_taus.size(); j++) {
                int i_next = (i+1)%jets.size(), j_next = (j+1)%vec_taus.size();
                PseudoJet combo1 = jets[i] + vec_taus[j].tau;
                PseudoJet combo2 = jets[i_next] + vec_taus[j_next].tau;
                double mass_diff = abs(combo1.m() - combo2.m());

                if (mass_diff < min_diff){
                    min_diff = mass_diff;
                    min_pairs.jet[0] = jets[i];
                    min_pairs.jet[1] = jets[i_next];
                    min_pairs.tau[0] = vec_taus[j];
                    min_pairs.tau[1] = vec_taus[j_next];
                }
            }
        }

        // cut 8 -- at least one tau, b pair has inv. mass > 250 GeV
        if ( (min_pairs.jet[0] + min_pairs.tau[0].tau).m() <= 250 &&
             (min_pairs.jet[1] + min_pairs.tau[1].tau).m() <= 250 ) {
            num_fail++;
            print_error(events,"failed cut 8: both tau + b jets have inv. mass < 250 GeV");
            continue;
        }


        // --------- EVENT PROCESSING OUTPUT --------- //
        if (DEBUG_MODE) {
            // find total number of particles clustered
            vector<PseudoJet> all_jets = cs.inclusive_jets();
            for (unsigned i = 0; i < all_jets.size(); i++)
                num_particles_clustered += all_jets[i].constituents().size();

            // verify that all particles were clustered; track progress in terminal
            if (num_particles_clustered == final_state - neutrinos - taus)
                print_success(events,"ALL PARTICLES CLUSTERED",
                              " ("+ to_string(num_particles_clustered) + ")");
            else print_warning(events,"NOT ALL PARTICLES CLUSTERED");
        }
        else print_success(events, "PASS");


        // --------- RECOMBINATORICS COMMENTED OUT FOR TIME BEING --------- //
        // // --------- PARTICLE RECONSTRUCTION --------- //
        // // note, for array:
        // // array[row][column] = array[j][i] ; j identifies the row, i the column
        // // array[row][column][layer] = array[k][j][i]
        // int n = jets.size()-1;
        // vector<double> inv_mass_w_vec;
        // // total two-jet combinatorics (for finding w)
        // double inv_mass_w_arr [n][n] = { };
        // for (int i = 0; i <= n-1; i++) {
        //     for (int j = i+1; j <= n; j++) {
        //         inv_mass_w_arr[j][i] = (jets[i] + jets[j]).m();
        //         inv_mass_w_vec.push_back(inv_mass_w_arr[j][i]);
        //     }
        // }
        // // unique two-jet combinatorics
        // double within_percent_of = 0.2;
        // // determine base 2-jet combinations
        // for (int i = 0; i<= n-1; i++) {                     // i = columns
        //     for (int j = i+1; j <= n; j++) {                // j = rows
        //         // determine comparison 2-jet combinations
        //         for (int x = i+1; x <= n-1; x++) {          // x = columns > i
        //             for (int y = 1; y <= n; y++) {          // y = all rows
        //                 // skip disallowed combinations
        //                 if (!(x<y) || i==x || i==y || j==x || j==y) continue;

        //                 // store only desired combinations
        //                 double m1 = (jets[i] + jets[j]).m();
        //                 double m2 = (jets[x] + jets[y]).m();

        //                 // -- info -- //
        //                 cout << "[" << i << "][" << j << "], ";
        //                 cout << "[" << x << "][" << y << "]: ";
        //                 // -- info -- //
        //                 if (abs(m1-m2) < (within_percent_of*max(m1,m2))) {
        //                     inv_mass_w_vec.push_back(m1);
        //                     inv_mass_w_vec.push_back(m2);
        //                 // -- info -- //
        //                     cout << "\033[32m" << m1 << "\t" << m2 << "\033[0m" << endl;
        //                 }
        //                 else {
        //                     cout << m1 << "\t" << m2 << endl;
        //                 // -- info -- //
        //                 }
        //             }
        //         }
        //     }
        // }

        // // three-jet combinatorics (for finding t)
        // double inv_mass_t_arr [n][n][n] = { };
        // vector<double> inv_mass_t_vec;
        // for (int i = 0; i <= n-2; i++) {
        //     for (int j = i+1; j <= n-1; j++) {
        //         for (int k = j+1; k <= n; k++) {
        //             inv_mass_t_arr[k][j][i] = (jets[i] + jets[j] + jets[k]).m();
        //             inv_mass_t_vec.push_back(inv_mass_t_arr[k][j][i]);
        //         }
        //     }
        // }

        // break;                                   // DEBUG -- run single event
    }

    // --------- FINAL STATUS INFO --------- //
    // write to terminal even if terminal output is suppressed
    cout << "Number of events processed: " << events << endl;
    cout << "Number of events passing cuts: " << (events - num_fail) << endl;
    cout << "Percent pass: " << (((events - num_fail))/double(events)*100) << endl;

    if (WRITE_TO_FILE) {
        runlog << "Number of events processed: " << events << endl;
        runlog << "Number of events passing cuts: " << (events - num_fail) << endl;
        runlog << "Percent pass: " << (((events - num_fail))/double(events)*100) << endl;

        cout << "Jet information written to " << runlog_filename << endl;
    }

    if (DEBUG_MODE) {
        cout << "Number of final state particles: " << total_final_state << endl;
        cout << "Number of final state neutrinos: " << total_neutrinos << endl;

        if (WRITE_TO_FILE) {
            runlog << "Number of final state particles: " << total_final_state << endl;
            runlog << "Number of final state neutrinos: " << total_neutrinos << endl;
        }
    }

    // --------- RUNTIME --------- //
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    cout << "Runtime (seconds): " << elapsed.count() << endl;
    if (WRITE_TO_FILE) runlog << "Runtime (seconds): " << elapsed.count() << endl;

    // --------- CLEANUP --------- //
    hepmc_file.close();
    if (WRITE_TO_FILE) runlog.close();

    return 0;
}


// --------- HELPER FUNCTIONS --------- //
PseudoJet get_jet(vector<string> delimited) {
    return PseudoJet(stof(delimited[PX]),stof(delimited[PY]),
                     stof(delimited[PZ]),stof(delimited[E]));
}

// debugging function
void print_jet(PseudoJet jet, string identifier) {
    int buffer = 15;
    string print = identifier;
    print.resize(buffer/2, ' ');
    string temp = "";
    for (int i = 0; i < 4; i++) {
        temp = to_string(jet[i]);
        if (temp[0] != '-') {
            print += ' ';
            temp.resize(buffer-1, ' ');
        }
        else temp.resize(buffer, ' ');
        print += temp;
    }
    cout << print << endl;
}

void _print_event(int event_number, string message, int color_message,
                  string other_info, int color_other_info) {
    if (WRITE_TO_TERM) {
        cout << "EVENT " << event_number << ":\t";
        cout << ("\033[" + to_string(color_message) + "m" + message + "\033[0m")
             << ("\033[" + to_string(color_other_info) + "m" + other_info + "\033[0m")
             << endl;
     }
     if (WRITE_TO_FILE) {
        runlog << "EVENT " << event_number << ":\t";
        runlog << message << other_info << endl;
     }
}

void print_success(int event_number, string message, string other_info,
                   int color_message, int color_other_info) {
    _print_event(event_number,message,color_message,other_info,color_other_info);
}

void print_error(int event_number, string message, int color) {
    if (!DEBUG_MODE && SUPPRESS_FAILURE_OUTPUT) return;
    _print_event(event_number,message,color);
}

void print_warning(int event_number, string message, int color) {
    _print_event(event_number,message,color);
}

vector<string> split_line(string line) {
    istringstream iss(line);
    return vector<string>((istream_iterator<string>(iss)),istream_iterator<string>());
}

double get_spatial_separation(PseudoJet jet1, PseudoJet jet2) {
    return sqrt(pow( (jet1.phi()-jet2.phi()), 2) + pow( (jet1.eta()-jet2.eta()), 2));
}

#include <fstream>
#include <iterator>
#include "fastjet/ClusterSequence.hh"
using namespace std;
using namespace fastjet;

// --------- GLOBAL FLAGS/CONSTANTS --------- //
// supress failure messages
const bool SUPRESS_FAILURE_OUTPUT = true;
// indices
const int PX = 3, PY = 4, PZ = 5, E = 6;

// --------- HELPER FUNCTION DECLARATIONS --------- //
PseudoJet get_jet(vector<string> delimited);
double get_spatial_separation(PseudoJet jet1, PseudoJet jet2);
void print_event(int event_number, string message, int color);
vector<string> split_line(string line);

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
    // --------- CONSTANTS --------- //
    // indices
    const int barcode = 1, pdg_code = 2, px = 3, py = 4, pz = 5, E = 6, gen_mass = 7,
              status = 8;
    // particle codes
    const string nu_e = "12", nu_mu = "14", nu_tau = "16", nu_tau_pr = "18";
    const string tau_p = "15", tau_m = "-15";
    // status flags
    const string final = "1";
    // jet clustering parameters
    const double R = 0.4;
    const int jet_pt_cutoff = 50;
    // null jet
    const PseudoJet null_jet = PseudoJet(0,0,0,0);
    // colors
    const int RED = 31, GREEN = 32, YELLOW = 33, BLUE = 34, PINK = 35, CYAN = 36;
    // input, output file names
    const string input_filename = "ditop_experiment4.hepmc",
                 jet_output_filename = "jet_output.txt",
                 w_output_filename = "w_output.txt",
                 t_output_filename = "t_output.txt";

    /// --------- VARIABLES --------- //
    // global counters
    int events = 0, num_fail = 0;
    int total_final_state = 0, total_neutrinos = 0;
    // event counters
    int final_state = 0, neutrinos = 0, num_particles_clustered = 0;
    // helper struct
    struct Tau {
        PseudoJet tau;
        bool is_tau_plus;
        PseudoJet vertex;
    };
    // event vertex
    PseudoJet current_vertex;     // bastardization of PsuedoJet to use functions
    // vectors
    vector<Tau> taus;
    PseudoJet tau_candidates[2];
    vector<PseudoJet> particles;
    vector<vector<PseudoJet>> passing_events;
    // FastJet setup
    JetDefinition jet_def(antikt_algorithm, R);

    // --------- FILE IO SETUP --------- //
    ifstream hepmc_file;
    ofstream jet_output, w_output, t_output;
    hepmc_file.open(input_filename);
    jet_output.open(jet_output_filename);
    w_output.open(w_output_filename);
    t_output.open(t_output_filename);



    // --------- HEPMC PARSING --------- //
    string line;
    getline(hepmc_file,line);
    while(line[0]!='E') {                       // get rid of lines until first event
        getline(hepmc_file,line);
    }

    while(!hepmc_file.eof()) {
        // -- info -- //
        events += 1;
        final_state = 0;
        neutrinos = 0;
        num_particles_clustered = 0;
        // -- info -- //

        // reset vectors
        particles.clear();
        taus.clear();
        tau_candidates[0] = null_jet; tau_candidates[1] = null_jet;

        getline(hepmc_file,line);               // gets line after event or neutrino
        // cycle through lines until next event
        while(!hepmc_file.eof() && line[0]!='E') {
            if(line[0]=='P'){                       // only want particles
                vector<string> delimited = split_line(line);

                // only want non-neutrino, final state particles
                // push candidates onto PseudoJet vector 'particles'
                if(delimited[status]==final) {
                    // -- info -- //
                    final_state += 1;               // count all final state particles
                    total_final_state += 1;
                    // -- info -- //

                    if(delimited[pdg_code] == nu_e || delimited[pdg_code] == nu_mu ||
                       delimited[pdg_code] == nu_tau || delimited[pdg_code] == nu_tau_pr)
                    {                               // don't store neutrinos
                        // -- info -- //
                        neutrinos += 1;
                        total_neutrinos += 1;
                        // -- info -- //
                        goto NextItem;
                    }

                    // store all taus + vertex, charge info in Tau vector 'taus'
                    if(delimited[pdg_code] == tau_p) {
                        taus.push_back( Tau{get_jet(delimited),true,current_vertex} );
                        goto NextItem;
                    }
                    if(delimited[pdg_code] == tau_m) {
                        taus.push_back( Tau{get_jet(delimited),false,current_vertex});
                        goto NextItem;
                    }

                    particles.push_back(get_jet(delimited));
                }
            }
            if(line[0]=='V') current_vertex.reset(get_jet(split_line(line)));

            NextItem:
            getline(hepmc_file,line);
        } // next event reached or eof


        // --------- CUTS ON TAUS --------- //
        // cut 1 -- at least two taus
        if (taus.size() < 2) {
            num_fail++;
            print_event(events,"failed cut 1 (lacks tau+ or tau-)",RED);
            continue;
        }

        // cuts 2-4
        bool vertex_match = false, opposite_charge = false, pt_pass = false,
             eta_pass = false, spatially_separated = false;
        for (int i = 0; i < taus.size() - 1; i++) {     // compare all combinations of taus
            for (int j = i + 1; j < taus.size(); j++) {
                // cut 4 -- same vertex
                if (taus[i].vertex == taus[j].vertex) {
                    vertex_match = true;
                    // cut 4 -- opposite charge
                    if (taus[i].is_tau_plus == !taus[j].is_tau_plus)
                        opposite_charge = true;
                    // cut 2 -- pt > 50
                    if (taus[i].tau.pt() > 50 && taus[j].tau.pt() > 50)
                        pt_pass = true;
                    // cut 3 -- |eta| < 2.3
                    if (abs(taus[i].tau.eta()) < 2.3 && abs(taus[j].tau.eta()) < 2.3)
                        eta_pass = true;
                    // cut 4 -- spatial separation
                    if (get_spatial_separation(taus[i].tau,taus[j].tau) > 0.5)
                        spatially_separated = true;

                    // usable taus go on array tau_candidates
                    tau_candidates[0] = taus[i].tau;
                    tau_candidates[1] = taus[j].tau;
                }
            }
        }
        if (!vertex_match) {
            num_fail++;
            print_event(events,"failed cut 4 (taus originate from different vertices)",RED);
            continue;
        }
        if (!opposite_charge) {
            num_fail++;
            print_event(events,"failed cut 4 (taus have same opposite charges)",RED);
            continue;
        }
        if (!pt_pass) {
            num_fail++;
            print_event(events,"failed cut 2 (taus have pt <= 50)",RED);
            continue;
        }
        if (!eta_pass) {
            num_fail++;
            print_event(events,"failed cut 3 (taus have |eta| >= 2.3)",RED);
            continue;
        }
        if (!spatially_separated) {
            num_fail++;
            print_event(events,"failed cut 4 (taus have separation <= 0.5)",RED);
            continue;
        }


        // --------- JET CLUSTERING --------- //
        // -- info -- //
        // if (events%50 == 0) cout << events << endl;
        // -- info -- //
        // if (events>5) continue;  // DEBUG -- run jet analysis for specific events

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
            print_event(events,"failed cut 5: less than two jets remaining after"
                                "eta, separation cuts",RED);
            continue;
        }

        // -- info -- //
        // write jets for event to jet_output
        // jet_output << "EVENT " << events
        //            << " (" << jets.size() << " jets w/ pt>" << jet_pt_cutoff << ")"
        //            << endl;
        // for (unsigned i = 0; i < jets.size(); i++) {
        //     jet_output << "jet " << i << " " << jets[i].E() << " " << jets[i].px()
        //                << " " << jets[i].py() << " " << jets[i].pz() << endl;
        // }

        // find total number of particles clustered
        vector<PseudoJet> all_jets = cs.inclusive_jets();
        for (unsigned i = 0; i < all_jets.size(); i++) {
            num_particles_clustered += all_jets[i].constituents().size();
        }

        // verify that all particles were clustered; track progress in terminal
        cout << "EVENT " << events << ": ";
        if (num_particles_clustered == final_state - neutrinos) {
            cout << "\033[32mALL PARTICLES CLUSTERED\033[0m ("
                 << num_particles_clustered << ")" << endl;
            // jet_output << "ALL PARTICLES CLUSTERED" << endl;
        }
        else {
            cout << "\033[31mNOT ALL PARTICLES CLUSTERED\033[0m ("
                 << num_particles_clustered << ")" << endl;
            // jet_output << "NOT ALL PARTICLES CLUSTERED" << endl;
        }
        // -- info -- //


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

        // // -- info -- //
        // // prints invariant masses from w,t searches to w,t output files
        // sort(inv_mass_w_vec.begin(),inv_mass_w_vec.end());
        // for (int i = 0; i < inv_mass_w_vec.size(); i++) {
        //     w_output << inv_mass_w_vec[i] << endl;
        // }

        // sort(inv_mass_t_vec.begin(),inv_mass_t_vec.end());
        // for (int i = 0; i < inv_mass_t_vec.size(); i++) {
        //     t_output << inv_mass_t_vec[i] << endl;
        // }
        // // -- info -- //

        // break;                                   // DEBUG -- run single event
    }

    // --------- FINAL STATUS INFO --------- //
    cout << "Number of events processed: " << events << endl;
    cout << "Number of final state particles: " << total_final_state << endl;
    cout << "Number of final state neutrinos: " << total_neutrinos << endl;
    cout << "Number of events passing cuts: " << (events - num_fail) << endl;
    cout << "Percent pass: " << (((events - num_fail))/double(events)*100) << endl;

    jet_output << "Number of events processed: " << events << endl;
    jet_output << "Number of final state particles: " << total_final_state << endl;
    jet_output << "Number of final state neutrinos: " << total_neutrinos << endl;

    cout << "Jet information written to " << jet_output_filename << endl;


    // --------- CLEANUP --------- //
    hepmc_file.close();
    jet_output.close();
    w_output.close();
    t_output.close();

    return 0;
}


// --------- HELPER FUNCTIONS --------- //
PseudoJet get_jet(vector<string> delimited) {
    return PseudoJet(stof(delimited[PX]),stof(delimited[PY]),
                     stof(delimited[PZ]),stof(delimited[E]));
}

double get_spatial_separation(PseudoJet jet1, PseudoJet jet2) {
    return sqrt(pow( (jet1.phi()-jet2.phi()), 2) + pow( (jet1.eta()-jet2.eta()), 2));
}

void print_event(int event_number, string message, int color = 37) {
    if (SUPRESS_FAILURE_OUTPUT) return;
    cout << "EVENT " << event_number << ": ";
    cout << ("\033[" + to_string(color) + "m" + message + "\033[0m") << endl;
}

vector<string> split_line(string line) {
    istringstream iss(line);
    return vector<string>((istream_iterator<string>(iss)),istream_iterator<string>());
}

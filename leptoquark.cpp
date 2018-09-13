#include <fstream>
#include <iterator>
#include "fastjet/ClusterSequence.hh"
using namespace std;
using namespace fastjet;

// --------- HELPER FUNCTIONS --------- //
void print_event(int event_number, string message, int color);

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
    int pdg_code = 2, px = 3, py = 4, pz = 5, E = 6, gen_mass = 7, status = 8;
    // neutrino codes
    string nu_e = "12", nu_mu = "14", nu_tau = "16", nu_tau_pr = "18";
    string tau_p = "15", tau_m = "-15";
    // status flags
    string final = "1";
    // counters
    int events = 0, final_state = 0, total_final_state = 0, neutrinos = 0,
        total_neutrinos = 0, num_particles_clustered = 0;
    int num_fail = 0;
    // event flags
    bool event_has_tau_p = false, event_has_tau_m = false;
    // jet clustering parameters
    double R = 0.4;
    int jet_pt_cutoff = 50;
    // input, output file names
    string input_filename = "ditop.hepmc",
    jet_output_filename = "jet_output.txt",
    w_output_filename = "w_output.txt",
    t_output_filename = "t_output.txt";
    // colors
    const int RED = 31, GREEN = 32, YELLOW = 33, BLUE = 34, PINK = 35, CYAN = 36;

    /// --------- INITIALIZATION --------- //
    vector<PseudoJet> particles;
    vector<PseudoJet> tau_plus;
    vector<PseudoJet> tau_minus;
    vector<vector<PseudoJet>> passing_events;
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

        // reset flags, vectors
        particles.clear();
        tau_plus.clear();
        tau_minus.clear();
        event_has_tau_p = false;
        event_has_tau_m = false;

        getline(hepmc_file,line);               // gets line after event or neutrino
        // cycle through lines until next event
        while(!hepmc_file.eof() && line[0]!='E') {
            if(line[0]=='P'){                       // only want particles
                // delimit line by space; results go into string vector 'delimited'
                istringstream iss(line);
                vector<string> delimited((istream_iterator<string>(iss)),istream_iterator<string>());

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

                    if(delimited[pdg_code] == tau_p) {
                        event_has_tau_p = true;
                        tau_plus.push_back(PseudoJet(stof(delimited[px]),stof(delimited[py]),
                                                 stof(delimited[pz]),stof(delimited[E])));
                    }
                    if(delimited[pdg_code] == tau_m) {
                        event_has_tau_m = true;
                        tau_minus.push_back(PseudoJet(stof(delimited[px]),stof(delimited[py]),
                                                 stof(delimited[pz]),stof(delimited[E])));
                    }

                    particles.push_back(PseudoJet(stof(delimited[px]),stof(delimited[py]),
                                                  stof(delimited[pz]),stof(delimited[E])));
                }
            }
            NextItem:
            getline(hepmc_file,line);
        } // next event reached or eof

        // CUT 1 -- final state must have two taus //
        if (!event_has_tau_p || !event_has_tau_m) {
            num_fail++;
            print_event(events,"failed cut 1 (lacks tau+ or tau-)",RED);
            continue;
        }

        // CUT 2 -- taus must have pt > 50 GeV, |eta| < 2.3 //
        bool plus_pass = false;
        for(vector<PseudoJet>::iterator it=tau_plus.begin(); it!=tau_plus.end(); ++it) {
            if ((*it).pt() > 50 && abs((*it).eta()) < 2.3) plus_pass = true;

        }
        if (!plus_pass) {
            num_fail++;
            print_event(events,"failed cut 2 (pt of tau+ < 50 GeV or |eta| > 2.3)",RED);
            continue;
        }

        // CUT 3 -- taus must have pt > 50 GeV, |eta| < 2.3 //
        bool minus_pass = false;
        for(vector<PseudoJet>::iterator it=tau_plus.begin(); it!=tau_plus.end(); ++it) {
            if ((*it).pt() > 50 && abs((*it).eta()) < 2.3) minus_pass = true;
        }
        if (!minus_pass) {
            num_fail++;
            print_event(events,"failed cut 3 (pt of tau- < 50 GeV or |eta| > 2.3)",RED);
            continue;
        }


        // --------- JET CLUSTERING --------- //
        // -- info -- //
        if (events%50 == 0) cout << events << endl;
        // -- info -- //
        // if (events>5) continue;  // DEBUG -- run jet analysis for specific events

        ClusterSequence cs(particles, jet_def);
        vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(jet_pt_cutoff));

        // -- info -- //
        // write jets for event to jet_output
        jet_output << "EVENT " << events
                   << " (" << jets.size() << " jets w/ pt>" << jet_pt_cutoff << ")"
                   << endl;
        for (unsigned i = 0; i < jets.size(); i++) {
            jet_output << "jet " << i << " " << jets[i].E() << " " << jets[i].px()
                       << " " << jets[i].py() << " " << jets[i].pz() << endl;
        }

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

void print_event(int event_number, string message, int color = 37) {
    cout << "EVENT " << event_number << ": ";
    cout << ("\033[" + to_string(color) + "m" + message + "\033[0m") << endl;
}

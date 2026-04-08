#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h> // pow, sqrt, lgamma
#include <cmath> //
#include <unistd.h> // getopt
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <algorithm>    // sort
#include "MCMC_v8.h"

using namespace std;

int run_contamination(
    int G, int K, int N_nei, int N, int N_MB, int N_tail, int n_record, int seed,
    const string& output_list,
    const string& data_name,
    const string& nei_name,
    const string& dist_name,
    const string& label_name,
    const string& cell_size_name,
    const string& MB_dir
) {
    cout << "Load all information successfully!" << endl;

    // loading observed data (gene expression)
    ifstream obs_data;
    string data_file;
    obs_data.open(data_name);
    if (!obs_data) {
        cout << "Unable to open the Y!" << endl;
        exit(1); // terminate with error
    }
    double **Y = new double *[G];
    for (int g = 0; g < G; ++g) {
        Y[g] = new double[N]();
        for (int i = 0; i < N; ++i) {
            obs_data >> Y[g][i];
        }
    }

    cout << "Load all expressions data successfully!" << endl;

    // loading data
    ifstream Y_nei;
    Y_nei.open(nei_name);
    if (!Y_nei) {
        cout << "Unable to open the file: nei_list.txt" << endl;
        exit(1); // terminate with error
    }

    // set-up spatial information
    int **neighbor_list = new int *[N];
    for (int i = 0; i < N; ++i) {
        neighbor_list[i] = new int[N_nei]();
        for (int j = 0; j < N_nei; ++j) {
            Y_nei >> neighbor_list[i][j];
        }
    }
    Y_nei.close();

    Y_nei.open(dist_name);
    if (!Y_nei) {
        cout << "Unable to open the file: nei_dist.txt" << endl;
        exit(1); // terminate with error
    }

    auto **neighbor_dist = new double *[N];
    for (int i = 0; i < N; ++i) {
        neighbor_dist[i] = new double[N_nei]();
        for (int j = 0; j < N_nei; ++j) {
            Y_nei >> neighbor_dist[i][j];
        }
    }
    Y_nei.close();

    // loading data
    ifstream size_file;
    size_file.open(cell_size_name);
    if (!size_file) {
        cout << "Unable to open the file: cell_size.txt" << endl;
        exit(1); // terminate with error
    }

    double *cell_size = new double[N]();
    for (int i = 0; i < N; ++i) {
        size_file >> cell_size[i];
    }

    size_file.close();

    int *length_array = new int[N_MB]();
    ifstream length_data;
    length_data.open(MB_dir + "length_vec.txt");
    if (!length_data) {
        cout << "Unable to open the length data!" << endl;
        exit(1); // terminate with error
    }
    for (int i = 0; i < N_MB; ++i) {
        length_data >> length_array[i];
    }

    length_data.close();

    int **MB_list = new int *[N_MB];
    for (int i = 0; i < N_MB; ++i) {
        MB_list[i] = new int[length_array[i]]();
    }
    int *tail_list = new int[N_tail]();
    ifstream MB_data;
    string MB_name;

    for (int i = 0; i < N_MB; ++i) {
        MB_name = MB_dir + "MB_" + to_string(i + 1) + ".txt";

        MB_data.open(MB_name);
        if (!MB_data) {
            cout << "Unable to open the MB data!" << endl;
            exit(1); // terminate with error
        }

        for (int j = 0; j < length_array[i]; ++j) {
            MB_data >> MB_list[i][j];
        }

        MB_data.close();

        cout << "Finish loading " << to_string(i + 1) << endl;
    }

    MB_name = MB_dir + "RMB.txt";
    MB_data.open(MB_name);
    if (!MB_data) {
        cout << "Unable to open the tail data!" << endl;
        exit(1); // terminate with error
    }
    for (int j = 0; j < N_tail; ++j) {
        MB_data >> tail_list[j];
    }

    MB_data.close();
    cout << "Finish loading tails" << endl;

    cout << "Successful loading!" << endl;



    ///////////////////////////////////////////////////////
    // 2. Set hyper-parameters and Initialize parameters //
    ///////////////////////////////////////////////////////
    cout << "Set initial values." << endl;
    int num_threads = omp_get_max_threads();
    cout << "The max number of core is " << num_threads << "." << endl;
    // Use std::vector to avoid manual memory management
    mt19937 *rng_array = new mt19937[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        rng_array[i].seed(seed + i); // Unique seed per thread
    }


    ///////////////////////////////////////////////////////
    // 3. Data Generation //
    ///////////////////////////////////////////////////////


    string out_file;

    //hyper-parameters
    //hyper-parameters
    double sigma_sq_a = 5, tau_mu_sq = 5,
            a_spike = 0.01, b_spike = 10, a_spike_spike = 0.001, a_beta = 0.0001, b_beta = 1,
            a_tau = 1, b_tau = 1,
            a_p = 7, b_p = 3,
            a_s = 1, b_s = 1, alpha = 0.5;



    // initialization: data-dependent initial values

    ifstream Y_label;
    Y_label.open(label_name);
    if (!Y_label) {
        cout << "Unable to open the file: Y_label.txt" << endl;
        exit(1); // terminate with error
    }

    // set-up spatial information
    int *Z_t = new int[N]();
    for (int i = 0; i < N; ++i) {
        Y_label >> Z_t[i];
//        Z_t[i] = rand_cate_uni(K, rng_array[0]);
//        cout << Z_t[i] <<" ";
    }
    Y_label.close();

    cout << "Successfully load labels!" << endl;


    double **gamma_t = new double *[G]();
    for (int g = 0; g < G; ++g) {
        gamma_t[g] = new double[2]();
        gamma_t[g][0] = 0.2;
    }

    double *beta_t = new double[2]();
    beta_t[0] = 0.01;
    beta_t[1] = 0;
    double **nei_exp_dist = new double *[N];
    for (int i = 0; i < N; ++i) {
        nei_exp_dist[i] = new double[N_nei]();
    }

    _update_exp_dist(N, N_nei, neighbor_dist, beta_t[0], nei_exp_dist);

    auto *sigma_sq_t = new double[G]();
    auto *pi_t = new double[K]();
    int *count_Z = new int[K]();

    auto **m = new double *[G];
    double **mu_t = new double *[G];
    for (int g = 0; g < G; ++g) {
        mu_t[g] = new double[K]();
        m[g] = new double[K]();
        for (int i = 0; i < N; ++i) {
            mu_t[g][Z_t[i]] += Y[g][i];
            count_Z[Z_t[i]]++;
        }

        for (int k = 0; k < K; ++k) {
            mu_t[g][k] = mu_t[g][k] / (count_Z[k] + 1);
            m[g][k] = mu_t[g][k];
        }
    }

//    ifstream mu_file;
//    mu_file.open("/lustre/project/Stat/s1155184072/MERFISH_v3/new_normal_dist/simulation/data_irregular/mu_true.txt");
//    if (!mu_file) {
//        cout << "Unable to open the file: mu_true.txt" << endl;
//        exit(1); // terminate with error
//    }
//
//    for (int g = 0; g < G; ++g) {
//        for (int k = 0; k < K; ++k) {
//            mu_file >> mu_t[g][k];
//        }
//    }
//    mu_file.close();
//
//    cout << "Successfully load mu!" << endl;

//    ifstream gamma_file;
//    gamma_file.open("/lustre/project/Stat/s1155184072/MERFISH_v3/new_normal_dist/simulation/data_irregular/gamma_true.txt");
//    if (!gamma_file) {
//        cout << "Unable to open the file: gamma_true.txt" << endl;
//        exit(1); // terminate with error
//    }
//
//    for (int g = 0; g < G; ++g) {
//        gamma_file >> gamma_t[g][0];
//    }
//    gamma_file.close();
//
//    cout << "Successfully load gamma!" << endl;


    cout << "Successful loading!" << endl;

    for (int k = 0; k < K; ++k) {
        pi_t[k] = (count_Z[k] + 0.0) / N;
//        pi_t[k] = 1.0 / K;
    }


#pragma omp parallel for
    for (int g = 0; g < G; ++g) {
        sigma_sq_t[g] = _update_sigma_sq(N, N_nei,
                                         a_s, b_s,
                                         Y[g], neighbor_list, cell_size,
                                         nei_exp_dist,
                                         Z_t,
                                         mu_t[g], gamma_t[g][0],
                                         rng_array[omp_get_thread_num()]);

        sigma_sq_t[g] = 0.01;
    }

    int output_size = 500;
    // allocate memory to store MCMC samplings
    int **Z_record = new int *[output_size / 25];
    double ***mu_record = new double **[output_size];
    int **I_beta_record = new int *[output_size]();
    auto **gamma_record = new double *[output_size];
    auto *beta_record = new double[output_size]();
    int *I_gamma_record = new int[output_size]();
    auto **sigma_sq_record = new double *[output_size];
    auto **pi_record = new double *[output_size];
    cout << "Finish Memories storage!" << endl;


    for (int t = 0; t < output_size; ++t) {

        mu_record[t] = new double *[G];
        for (int g = 0; g < G; ++g) {
            mu_record[t][g] = new double[K]();
        }

        gamma_record[t] = new double[G]();
        I_beta_record[t] = new int[G]();
        sigma_sq_record[t] = new double[G]();
        pi_record[t] = new double[K]();

    }

    for (int t = 0; t < output_size / 25; ++t) {
        Z_record[t] = new int[N]();
    }

    int iter_index = 0;
    int Z_index = 0;
    int *Z_last = new int[N]();
    int *max_index = new int[K]();

    double **con_table = new double *[K];
    for (int k = 0; k < K; ++k) {
        con_table[k] = new double[K]();
    }

    auto **nei_exp_dist_r = new double *[N];
    auto **nei_exp_dist_t = new double *[N];
    for (int i = 0; i < N; ++i) {
        nei_exp_dist_r[i] = new double[N_nei]();
        nei_exp_dist_t[i] = new double[N_nei]();
    }

    cout << "Finish initial values!" << endl;

    //////////////////////
    // 4. MCMC sampling //
    //////////////////////
    double start = omp_get_wtime();
    cout << "Start MCMC sampling." << endl;

    for (int iter = 0; iter < n_record; ++iter) {
#pragma omp parallel for
        for (int g = 0; g < G; ++g) {
            //update mu
            _update_mu(N, K, N_nei,
                       m[g], tau_mu_sq,
                       Y[g], neighbor_list, cell_size,
                       nei_exp_dist,
                       Z_t,
                       gamma_t[g][0], sigma_sq_t[g],
                       rng_array[omp_get_thread_num()],
                       mu_t[g]);
        }


#pragma omp parallel for
        for (int g = 0; g < G; ++g) {
            //update sigma_sq
            sigma_sq_t[g] = _update_sigma_sq(N, N_nei,
                                             a_s, b_s,
                                             Y[g], neighbor_list, cell_size,
                                             nei_exp_dist,
                                             Z_t, mu_t[g], gamma_t[g][0],
                                             rng_array[omp_get_thread_num()]);
        }


#pragma omp parallel for
        for (int g = 0; g < G; ++g) {
            _update_gamma(N, N_nei,
                          1, 1,
                          Y[g], neighbor_list, cell_size,
                          nei_exp_dist,
                          Z_t,
                          mu_t[g], sigma_sq_t[g],
                          rng_array[omp_get_thread_num()],
                          gamma_t[g]);
        }

        _update_beta(N, G, N_nei,
                     a_beta, b_beta, Y, neighbor_list, neighbor_dist, cell_size,
                     nei_exp_dist_t, nei_exp_dist_r,
                     Z_t,
                     mu_t, gamma_t, sigma_sq_t,
                     beta_t,
                     rng_array[0]);

        _update_exp_dist(N, N_nei,
                         neighbor_dist, beta_t[0],
                         nei_exp_dist);

        // update Z
        if (((iter + 1) % 25 == 0)) {
#pragma omp parallel for
            for (int i = 0; i < N; ++i) {
                Z_last[i] = Z_t[i];
            }

            _update_Z_MB(G, K, N_nei,
                         Y, neighbor_list, cell_size,
                         nei_exp_dist,
                         N_MB, length_array, MB_list, N_tail, tail_list,
                         mu_t, gamma_t, sigma_sq_t, pi_t,
                         rng_array,
                         Z_t);

            // clear
            for (int k = 0; k < K; ++k) {
                for (int kk = 0; kk < K; ++kk) {
                    con_table[k][kk] = 0;
                }
            }

            // compute the contingency table
#pragma omp parallel for
            for (int i = 0; i < N; ++i) {
                con_table[Z_last[i]][Z_t[i]]++;
            }

            // find the corresponding index based on last record
            for (int k = 0; k < K; ++k) {
                max_index[k] = my_max(con_table[k], K);
            }

            // reorder and record
#pragma omp parallel for
            for (int i = 0; i < N; ++i) {
                for (int k = 0; k < K; ++k) {
                    if (Z_t[i] == max_index[k]) {
                        Z_t[i] = k;
                    }
                }
 // record Z
                Z_record[Z_index][i] = Z_t[i];
            }

            Z_index++;
        }

        //update pi
        _update_pi(N, K,
                   alpha,
                   Z_t,
                   rng_array[0],
                   pi_t);
        //record samplings

#pragma omp parallel for
        for (int g = 0; g < G; ++g) {
            gamma_record[iter_index][g] = gamma_t[g][0];
            beta_record[iter_index] = beta_t[0];
            I_beta_record[iter_index][g] = int(gamma_t[g][1]);
            I_gamma_record[iter_index] = int(beta_t[1]);

            for (int k = 0; k < K; ++k) {
                mu_record[iter_index][g][k] = mu_t[g][k];
            }
            sigma_sq_record[iter_index][g] = sigma_sq_t[g];
        }


        //record samplings of pi
        for (int k = 0; k < K; ++k) {
            pi_record[iter_index][k] = pi_t[k];
        }

        iter_index++;
//        cout << "Finish all" << endl;

        if ((iter + 1) % output_size == 0) {
            /// Output now
            ofstream output_temp;
            cout << "output samplings now" << endl;


            out_file = output_list + "mu_" + to_string(iter + 1) + ".txt";
            output_temp.open(out_file.c_str(), ios::out | ios::app);
            for (int g = 0; g < G; ++g) {
                for (int k = 0; k < K; ++k) {
                    for (int iter_o = 0; iter_o < output_size; ++iter_o) {
                        output_temp << mu_record[iter_o][g][k];
                        output_temp << " ";
                    }
                    output_temp << endl;
                }
            }
            output_temp.close();
            cout << "mu finished!" << endl;

            out_file = output_list + "gamma_" + to_string(iter + 1) + ".txt";
            output_temp.open(out_file.c_str(), ios::out | ios::app);
            for (int g = 0; g < G; ++g) {
                for (int iter_o = 0; iter_o < output_size; ++iter_o) {
                    output_temp << gamma_record[iter_o][g];
                    output_temp << " ";
                }
                output_temp << endl;

            }
            output_temp.close();
            cout << "gamma finished!" << endl;

            out_file = output_list + "beta_" + to_string(iter + 1) + ".txt";
            output_temp.open(out_file.c_str(), ios::out | ios::app);
            for (int iter_o = 0; iter_o < output_size; ++iter_o) {
                output_temp << beta_record[iter_o];
                output_temp << endl;
            }
            output_temp.close();
            cout << "beta finished!" << endl;

            out_file = output_list + "I_beta_" + to_string(iter + 1) + ".txt";
            output_temp.open(out_file.c_str(), ios::out | ios::app);
            for (int g = 0; g < G; ++g) {
                for (int iter_o = 0; iter_o < output_size; ++iter_o) {
                    output_temp << I_beta_record[iter_o][g];
                    output_temp << " ";
                }
                output_temp << endl;
            }

            output_temp.close();
            cout << "I_beta finished!" << endl;

            out_file = output_list + "I_gamma_" + to_string(iter + 1) + ".txt";
            output_temp.open(out_file.c_str(), ios::out | ios::app);
            for (int iter_o = 0; iter_o < output_size; ++iter_o) {
                output_temp << I_gamma_record[iter_o];
                output_temp << endl;
            }

            output_temp.close();
            cout << "I_gamma finished!" << endl;


            out_file = output_list + "sigma_sq_" + to_string(iter + 1) + ".txt";
            output_temp.open(out_file.c_str(), ios::out | ios::app);
            for (int g = 0; g < G; ++g) {
                for (int iter_o = 0; iter_o < output_size; ++iter_o) {
                    output_temp << sigma_sq_record[iter_o][g];
                    output_temp << " ";
                }
                output_temp << endl;

            }
            output_temp.close();
            cout << "sigma square finished!" << endl;

//            if (iter < n_record / 2) {
//
//            }
            out_file = output_list + "Z_" + to_string((iter + 1) / 25) + ".txt";
            output_temp.open(out_file.c_str(), ios::out | ios::app);
            for (int i = 0; i < N; ++i) {
                for (int iter_o = 0; iter_o < output_size / 25; ++iter_o) {
                    output_temp << Z_record[iter_o][i];
                    output_temp << " ";
                }
                output_temp << endl;
            }

            output_temp.close();
            cout << "Z finished!" << endl;

            out_file = output_list + "pi_" + to_string(iter + 1) + ".txt";
            output_temp.open(out_file.c_str(), ios::out | ios::app);

            for (int k = 0; k < K; ++k) {
                for (int iter_o = 0; iter_o < output_size; ++iter_o) {
                    output_temp << pi_record[iter_o][k];
                    output_temp << " ";
                }
                output_temp << endl;
            }

            output_temp.close();
            cout << "pi finished!" << endl;

            iter_index = 0;
            Z_index = 0;
        }

        if (iter % 100 == 0) {
            double end = omp_get_wtime();
            cout << "Iteration " << iter << endl;
            cout << "Already used " << end - start << endl;
        }

    }

//    end_MCMC = chrono::system_clock::now();
//    elapsed_seconds_MCMC = end_MCMC - start_MCMC;
    double end = omp_get_wtime();
    cout << "Time of " << n_record << " iterations of MCMC sampling is: "
         << (end - start) / 60.0 << "min" << endl;


    cout << "MCMC finished!" << endl;

    for (int g = 0; g < G; ++g) {
        delete[] Y[g];
        delete[] mu_t[g];
    }

    delete[] Y;
    delete[] mu_t;

    delete[] sigma_sq_t;


    delete[] pi_t;
    delete[] count_Z;

    for (int i = 0; i < N; ++i) {
        delete[] neighbor_list[i];

    }

    delete[] neighbor_list;


    delete[] Z_t;


    for (int iter = 0; iter < output_size; ++iter) {
        for (int g = 0; g < G; ++g) {
            delete[] mu_record[iter][g];
        }

        delete[] sigma_sq_record[iter];
        delete[] pi_record[iter];
    }

    for (int iter = 0; iter < output_size; ++iter) {
        delete[] mu_record[iter];
    }

    delete[] mu_record;
    delete[] sigma_sq_record;
    delete[] pi_record;

    return 0;
}


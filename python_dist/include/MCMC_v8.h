//
// Created by s1155184072 on 2025/3/18.
//

#ifndef MERIFISH_MCMC_V6_1_H
#define MERIFISH_MCMC_V6_1_H

#include <cmath> //
#include <random>
#include <omp.h>
#include <iostream>
#include <fstream>

using namespace std;
const double pi = acos(-1.0);

double my_max(const double *_seq, int _length) {
    double temp_max = _seq[0];
    for (int i = 1; i < _length; ++i) {
        if (temp_max < _seq[i]) {
            temp_max = _seq[i];
        }
    }
    return temp_max;
}

int rand_cate(const double *_prop, mt19937 &rng) {
    uniform_real_distribution<double> uniform_dist(0.0, 1.0);

    int res = 0;
    double u = uniform_dist(rng);
    while (u > _prop[res]) {
        u = u - _prop[res];
        res++;
    }
    return res;
}

int rand_cate_uni(int _K, std::mt19937 &rng) {
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    int res = 0;
    auto *prob = new double[_K];
    for (int k = 0; k < _K; ++k) {
        prob[k] = 1.0 / _K;
    }

    double u = uniform_dist(rng);
    while (u > prob[res]) {
        u = u - prob[res];
        res++;
    }
    delete[] prob;
    return res;
}

double rgamma(double alpha, double beta, mt19937 &rng) {
    // generate random number from X~gamma(alpha,beta)
    // f(x) = 1/(gamma(alpha)*beta^alpha) * x^(alpha-1) * exp(-x/beta)
    // x > 0, alpha > 0, beta > 0
    // E(X) = alpha*beta, Var(X) = alpha*beta^2
    uniform_real_distribution<double> uniform_dist(0.0, 1.0);


    double x = 0.0;
    double delta, v0, v1, v2, v3, psi = 1.0, nu = 1.0e10;

    for (int i = 0; i < floor(alpha); i++)
        x = x + -log(uniform_dist(rng));

    if (alpha > floor(alpha)) { // alpha has fractional part
        delta = alpha - floor(alpha);
        while (nu > pow(psi, delta - 1.0) * exp(-psi)) {
            v1 = uniform_dist(rng);
            v2 = uniform_dist(rng);
            v3 = uniform_dist(rng);
            v0 = exp(1.0) / (exp(1.0) + delta);
            if (v1 <= v0) {
                psi = pow(v2, 1.0 / delta);
                nu = v3 * pow(psi, delta - 1.0);
            } else {
                psi = 1.0 - log(v2);
                nu = v3 * exp(-psi);
            }
        }
        x = x + psi;
    }

    return (beta * x);
}


void rand_Dir(double *_xi, int _K, mt19937 &rng, double *_pi) {

    auto *rn_gam = new double[_K];
    double sum_gam = 0.0;
    for (int k = 0; k < _K; k++) {
        rn_gam[k] = rgamma(_xi[k], 1.0, rng);
        sum_gam += rn_gam[k];
    }
    for (int k = 0; k < _K; k++) {
        _pi[k] = rn_gam[k] / sum_gam;
    }

    delete[] rn_gam;
}

double fast_beta(double alpha, double beta, mt19937 &rng) {
    double x = rgamma(alpha, 1, rng);
    double y = rgamma(beta, 1, rng);
    return (x / (x + y));
}


double rnorm(double mu, double sigma, mt19937 &rng) {
    // generate random number from X~N(mu,sigma)
    // f(x) = 1/sqrt(2*pi*si^2) * exp(-0.5*((x-mu)/si)^2)
    // -infty < x < infty, -infty < mu < infty, sigma > 0
    // E(X) = mu, Var(X) = sigma^2
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    double u, v, z, x;

    u = uniform_dist(rng);
    v = uniform_dist(rng);

    z = sqrt(-2 * log(u)) * cos(2 * pi * v);
    x = mu + z * sigma;

    return (x);
}

void _update_exp_dist(int _N, int _N_nei,//dimension
                      double **_nei_dist, double _beta,
                      double **_exp_nei_dist
) {
#pragma omp parallel for
    for (int i = 0; i < _N; ++i) {
        for (int j = 0; j < _N_nei; ++j) {
            if (_nei_dist[i][j] >= 0) {
                _exp_nei_dist[i][j] = exp(-_beta * _nei_dist[i][j]);
            } else {
                break;
            }
        }
    }
}

void _update_mu(int _N, int _K, int _N_nei,//dimension
                double *_m, double _tau_mu_sq,//prior
                double *_Y, int **_nei_list, double *_cell_size,// observed
                double **_nei_exp_dist,
                int *_Z, // latent
                double _gamma, double _sigma_sq, //parameter
                mt19937 &rng,
                double *_mu) {

    double post_mu, post_sigma_sq;
    for (int k = 0; k < _K; ++k) {
        post_mu = 0;
        post_sigma_sq = 0;

#pragma omp parallel for reduction(+:post_mu, post_sigma_sq)
        for (int i = 0; i < _N; ++i) {
            double same_nei_sum = 0;
            double diff_nei_sum = 0;
            int count = 0;

            for (int j = 0; j < _N_nei; ++j) {
                int nei_idx = _nei_list[i][j];
                double nei_dist = _nei_exp_dist[i][j];
                if (nei_idx >= 0) {
                    if (_Z[nei_idx] == k) {
                        same_nei_sum += _cell_size[nei_idx] * nei_dist;
                        count++;
                    } else {
                        diff_nei_sum += _mu[_Z[nei_idx]] * _cell_size[nei_idx] * nei_dist;
                    }

                } else {
                    break;
                }
            }


            if (_Z[i] == k) {
                double factor = _gamma * same_nei_sum + _cell_size[i];
                post_mu += factor * (_Y[i] - _gamma * diff_nei_sum) / _sigma_sq;
                post_sigma_sq += factor * factor / _sigma_sq;
            } else if (count != 0) {
                same_nei_sum = _gamma * same_nei_sum;
                post_mu += same_nei_sum * (_Y[i] - _mu[_Z[i]] * _cell_size[i] - _gamma * diff_nei_sum) / _sigma_sq;
                post_sigma_sq += same_nei_sum * same_nei_sum / _sigma_sq;
            }
        }

        post_sigma_sq = 1.0 / (post_sigma_sq + 1.0 / _tau_mu_sq);
        post_mu = (post_mu + _m[k] / _tau_mu_sq) * post_sigma_sq;

        double normal_sample = rnorm(0, 1, rng);
        _mu[k] = post_mu + sqrt(post_sigma_sq) * normal_sample;
    }
}


void _update_gamma(int _N, int _N_nei,// dimension
                   double _a_gamma, double _b_gamma,// prior
                   double *_Y, int **_neighbor_list, double *_cell_size,//observed and neighbor index
                   double **_nei_exp_dist,
                   int *_Z, //latent
                   double *_mu, double _sigma_sq, //parameters
                   mt19937 &rng,
                   double *_gamma) {

    uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    double log_r = 0;
    double gamma_t = _gamma[0];
    double norm_tune = 3;
    double beta_tune = 5;
    double res = rnorm(gamma_t, gamma_t / norm_tune, rng);
    int indicator = 0;
    if (res <= 0 | res >= 1) {
        res = fast_beta(gamma_t / beta_tune, (1 - gamma_t) / beta_tune, rng);
        indicator = 1;
    }

    int I_beta = int(_gamma[1]);

#pragma omp parallel for reduction(+:log_r)
    for (int i = 0; i < _N; ++i) {
        double nei_sum = 0;


        for (int j = 0; j < _N_nei; ++j) {
            int nei_idx = _neighbor_list[i][j];
            double nei_dist = _nei_exp_dist[i][j];
            if (nei_idx >= 0) {
                nei_sum += _mu[_Z[nei_idx]] * _cell_size[nei_idx] * nei_dist;
            } else {
                break;
            }
        }


        //nominator
        double no_factor = _Y[i] - _mu[_Z[i]] * _cell_size[i] - res * nei_sum;
        log_r += -no_factor * no_factor / (2 * _sigma_sq);
        //denominator
        double de_factor = _Y[i] - _mu[_Z[i]] * _cell_size[i] - gamma_t * nei_sum;
        log_r += de_factor * de_factor / (2 * _sigma_sq);
    }


    //prior
    log_r += _a_gamma * log(res) + _b_gamma * log(1 - res)
             - _a_gamma * log(gamma_t) - _b_gamma * log(1 - gamma_t);


    //proposal
    if (indicator == 0 && I_beta == 0) {
        log_r += -log(res / norm_tune) - (gamma_t - res) * (gamma_t - res) / (2 * res / norm_tune * res / norm_tune)
                 + log(gamma_t / norm_tune) +
                 (gamma_t - res) * (gamma_t - res) / (2 * gamma_t / norm_tune * gamma_t / norm_tune);

    } else if (indicator == 0 && I_beta == 1) {
        log_r += -log(2 * pi) / 2 - log(res / norm_tune) -
                 (gamma_t - res) * (gamma_t - res) / (2 * res / norm_tune * res / norm_tune)
                 - lgamma(1.0 / beta_tune) + lgamma(gamma_t / beta_tune) + lgamma(1 - gamma_t / beta_tune) -
                 gamma_t / beta_tune * log(res) -
                 (1 - gamma_t / beta_tune) * log(1 - res);

    } else if (indicator == 1 && I_beta == 1) {
        log_r += lgamma(1.0 / beta_tune) - lgamma(res / beta_tune) - lgamma(1 - res / beta_tune) +
                 res / beta_tune * log(gamma_t) +
                 (1 - res / beta_tune) * log(1 - gamma_t)
                 - lgamma(1.0 / beta_tune) + lgamma(gamma_t / beta_tune) + lgamma(1 - gamma_t / beta_tune) -
                 gamma_t / beta_tune * log(res) -
                 (1 - gamma_t / beta_tune) * log(1 - res);
    } else {
        log_r += lgamma(1.0 / beta_tune) - lgamma(res / beta_tune) - lgamma(1 - res / beta_tune) +
                 res / beta_tune * log(gamma_t) +
                 (1 - res / beta_tune) * log(1 - gamma_t)
                 + log(gamma_t / norm_tune) +
                 (gamma_t - res) * (gamma_t - res) / (2 * gamma_t / norm_tune * gamma_t / norm_tune);
    }

    _gamma[1] = indicator;

    if (log_r > log(uniform_dist(rng))) {
        _gamma[0] = res;
    }
}

void _update_beta(int _N, int _G, int _N_nei,// dimension
                  double _a, double _b,// prior
                  double **_Y, int **_neighbor_list, double **_neighbor_dist, double *_cell_size,//observed
                  double **_nei_exp_dist_t, double **_nei_exp_dist_r,
                  int *_Z,
                  double **_mu, double **_gamma, double *_sigma_sq, // parameters
                  double *_beta,
                  mt19937 &rng) {

    uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    double log_r = 0;
    double beta_t = _beta[0];
    double norm_tune = 5;
    double gamma_tune = 3;
    double res = rnorm(beta_t, beta_t / norm_tune, rng);
    int indicator = 0;
    if (res < 0) {
        res = rgamma(beta_t / gamma_tune, gamma_tune, rng);
        indicator = 1;
    }

    int I_gamma = int(_beta[1]);

    _update_exp_dist(_N, _N_nei, _neighbor_dist, beta_t, _nei_exp_dist_t);
    _update_exp_dist(_N, _N_nei, _neighbor_dist, res, _nei_exp_dist_r);

#pragma omp parallel for reduction(+:log_r)
    for (int i = 0; i < _N; ++i) {
        for (int g = 0; g < _G; ++g) {
            double numerator_nei_sum = 0;
            double denominator_nei_sum = 0;
            for (int j = 0; j < _N_nei; ++j) {
                int nei_idx = _neighbor_list[i][j];
                double num_nei_dist = _nei_exp_dist_r[i][j];
                double den_nei_dist = _nei_exp_dist_t[i][j];
                if (nei_idx >= 0) {
                    numerator_nei_sum += _mu[g][_Z[nei_idx]] * _cell_size[nei_idx] * num_nei_dist;
                    denominator_nei_sum += _mu[g][_Z[nei_idx]] * _cell_size[nei_idx] * den_nei_dist;
                } else {
                    break;
                }

            }
            double numerator_factor = (_Y[g][i] - _mu[g][_Z[i]] * _cell_size[i] - _gamma[g][0] * numerator_nei_sum);
            double denominator_factor = (_Y[g][i] - _mu[g][_Z[i]] * _cell_size[i] - _gamma[g][0] * denominator_nei_sum);
            log_r += -numerator_factor * numerator_factor / (2 * _sigma_sq[g]) +
                     denominator_factor * denominator_factor / (2 * _sigma_sq[g]);
        }
    }
    //prior
    log_r += (_a - 1) * (log(res) - log(beta_t)) + _b * (-res + beta_t);

    //proposal
    if (indicator == 0 && I_gamma == 0) {
        log_r += -log(res / norm_tune) - (beta_t - res) * (beta_t - res) / (2 * res / norm_tune * res / norm_tune)
                 + log(beta_t / norm_tune) +
                 (beta_t - res) * (beta_t - res) / (2 * beta_t / norm_tune * beta_t / norm_tune);

    } else if (indicator == 0 && I_gamma == 1) {
        log_r += -log(2 * pi) / 2 - log(res / norm_tune) -
                 (beta_t - res) * (beta_t - res) / (2 * res / norm_tune * res / norm_tune)
                 + lgamma(beta_t / gamma_tune) - beta_t / gamma_tune * log(gamma_tune) -
                 (beta_t / gamma_tune - 1) * log(res) + gamma_tune * res;

    } else if (indicator == 1 && I_gamma == 1) {
        log_r += -lgamma(res / gamma_tune) + res / gamma_tune * log(gamma_tune) +
                 (res / gamma_tune - 1) * log(beta_t) - gamma_tune * beta_t
                 + lgamma(beta_t / gamma_tune) - beta_t / gamma_tune * log(gamma_tune) -
                 (beta_t / gamma_tune - 1) * log(res) + gamma_tune * res;
    } else {
        log_r += -lgamma(res / gamma_tune) + res / gamma_tune * log(gamma_tune) +
                 (res / gamma_tune - 1) * log(beta_t) - gamma_tune * beta_t
                 + log(beta_t / norm_tune) +
                 (beta_t - res) * (beta_t - res) / (2 * beta_t / norm_tune * beta_t / norm_tune);
    }

    _beta[1] = indicator;
    if (log_r > log(uniform_dist(rng))) {
        _beta[0] = res;
    }
}

double _update_sigma_sq(int _N, int _N_nei,// dimension
                        double _a, double _b, // prior
                        double *_Y, int **_neighbor_list, double *_cell_size,//observed
                        double **_nei_exp_dist,
                        int *_Z, // observed + latent
                        double *_mu, double _gamma, // parameters
                        mt19937 &rng) {

    double a_post, b_post = 0;
    double sigma_sq;

#pragma omp parallel for reduction(+:b_post)
    for (int i = 0; i < _N; ++i) {
        double nei_sum = 0;

        for (int j = 0; j < _N_nei; ++j) {
            int nei_idx = _neighbor_list[i][j];
            double nei_dist = _nei_exp_dist[i][j];
            if (nei_idx >= 0) {
                nei_sum += _mu[_Z[nei_idx]] * _cell_size[nei_idx] * nei_dist;

            } else {
                break;
            }
        }
        b_post += pow(_Y[i] - _mu[_Z[i]] * _cell_size[i] - _gamma * nei_sum, 2);
    }

    a_post = _a + _N * 0.5;
    b_post = _b + 0.5 * b_post;

    sigma_sq = 1.0 / rgamma(a_post, 1.0 / b_post, rng);
    return sigma_sq;
}

void _update_Z(int _N, int _G, int _K, int _N_nei,// dimension
               double **_Y, int **_neighbor_list, // observed
               double *_nu, double **_mu, double *_alpha, double **_exp_values, double *_sigma_sq,
               double *_pi, //parameters
               mt19937 &rng,
               int *_Z) {


    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    double start, end;
    start = omp_get_wtime();
    for (int i = 0; i < _N; ++i) {

        double prob_max, prob_norm;
        auto *prob = new double[_K]();
#pragma omp parallel for shared(prob)
        for (int k = 0; k < _K; ++k) {
            double self_sum = 0;
            // self sum
#pragma omp parallel for reduction(+:self_sum)
            for (int g = 0; g < _G; ++g) {
                double nei_sum = 0;

#pragma omp simd
                for (int j = 0; j < _N_nei; ++j) {
                    int nei_idx = _neighbor_list[i][j];
                    if (nei_idx >= 0) {
                        nei_sum += _alpha[g] * _exp_values[i][j]
                                   * (_nu[g] + _mu[g][_Z[nei_idx]]);
                    }
                }
                double self_factor = _Y[g][i] - _nu[g] - _mu[g][k] - nei_sum;
                self_sum += self_factor * self_factor /
                            (2.0 * _sigma_sq[g]);

            }

            //neighbor sum
            double nei_self_sum = 0;

#pragma omp parallel for reduction(+:nei_self_sum)
            for (int j = 0; j < _N_nei; ++j) {
                int nei_index = _neighbor_list[i][j];

                if (nei_index >= 0) {

#pragma omp simd
                    for (int g = 0; g < _G; ++g) {
                        double nei_nei_sum = 0;

                        for (int r = 0; r < _N_nei; ++r) {
                            int nei_nei_idx = _neighbor_list[nei_index][r];
                            if ((nei_nei_idx >= 0) && (nei_nei_idx != i)) {
                                nei_nei_sum += _alpha[g] * _exp_values[nei_index][r] *
                                               (_nu[g] + _mu[g][_Z[nei_nei_idx]]);
                            }
                        }
                        double nei_self_factor = _Y[g][nei_index] - _nu[g] - _mu[g][_Z[nei_index]] -
                                                 _alpha[g] * _exp_values[i][j] * (_nu[g] + _mu[g][k])
                                                 - nei_nei_sum;
                        nei_self_sum += nei_self_factor * nei_self_factor /
                                        (2.0 * _sigma_sq[g]);
                    }
                }
            }

            prob[k] = log(_pi[k]) - self_sum - nei_self_sum;
        }

        prob_max = my_max(prob, _K);
        prob_norm = 0;
#pragma omp simd
        for (int k_t = 0; k_t < _K; ++k_t) {
            prob[k_t] = exp(prob[k_t] - prob_max);
            prob_norm += prob[k_t];
        }

#pragma omp simd
        for (int k_t = 0; k_t < _K; ++k_t) {
            prob[k_t] = prob[k_t] / prob_norm;
        }

        // Sample from categorical distribution deterministically
        double u = uniform_dist(rng);
        double cumulative_prob = 0.0;
        for (int k_t = 0; k_t < _K; ++k_t) {
            cumulative_prob += prob[k_t];
            if (u < cumulative_prob) {
                _Z[i] = k_t;
                break;
            }
        }

        if ((i + 1) % 10000 == 0) {
            cout << "Now it's " << i + 1 << "!" << endl;
            end = omp_get_wtime();
            double parallel_time = end - start;
            cout << "Time is " << parallel_time << "!" << endl;
        }

        delete[] prob;
    } //End for-i loop
    end = omp_get_wtime();
    double parallel_time = end - start;
    cout << "All time is " << parallel_time << "!" << endl;


}

void _update_Z_MB(int _G, int _K, int _N_nei,// dimension
                  double **_Y, int **_neighbor_list, double *_cell_size,// observed
                  double **_nei_exp_dist,
                  int _N_MB, int *_length_MB, int **_MB_list, int _tail_length, int *_tail_list, // Markov Blanket
                  double **_mu, double **_gamma, double *_sigma_sq,
                  double *_pi, //parameters
                  mt19937 *&rng,
                  int *_Z) {


    uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    double start;
    start = omp_get_wtime();

    for (int m = 0; m < _N_MB; ++m) {

#pragma omp parallel for
        for (int i = 0; i < _length_MB[m]; ++i) {
            int MB_idx = _MB_list[m][i];
            double prob_max, prob_norm;
            auto *prob = new double[_K]();

//#pragma omp parallel for shared(prob)
            for (int k = 0; k < _K; ++k) {
                double self_sum = 0;
                // self sum
#pragma omp parallel for reduction(+:self_sum)
                for (int g = 0; g < _G; ++g) {
                    double nei_sum = 0;

                    for (int j = 0; j < _N_nei; ++j) {
                        int nei_idx = _neighbor_list[MB_idx][j];
                        double nei_dist = _nei_exp_dist[MB_idx][j];
                        if (nei_idx >= 0) {
                            nei_sum += _mu[g][_Z[nei_idx]] * _cell_size[nei_idx] * nei_dist;
                        } else {
                            break;
                        }
                    }
                    double self_factor = _Y[g][MB_idx] - _mu[g][k] * _cell_size[MB_idx] - _gamma[g][0] * nei_sum;
                    self_sum += self_factor * self_factor / (2.0 * _sigma_sq[g]);
                }

                //neighbor sum
                double nei_self_sum = 0;

//#pragma omp parallel for reduction(+:nei_self_sum)
                for (int j = 0; j < _N_nei; ++j) {
                    int nei_index = _neighbor_list[MB_idx][j];

                    if (nei_index >= 0) {


                        for (int g = 0; g < _G; ++g) {
                            double nei_nei_sum = 0;

                            for (int r = 0; r < _N_nei; ++r) {
                                int nei_nei_idx = _neighbor_list[nei_index][r];
                                double nei_nei_dist = _nei_exp_dist[nei_index][r];

                                if ((nei_nei_idx >= 0) && (nei_nei_idx != MB_idx)) {
                                    nei_nei_sum += _mu[g][_Z[nei_nei_idx]] * _cell_size[nei_nei_idx] * nei_nei_dist;
                                } else if (nei_nei_idx < 0) {
                                    break;
                                }
                            }

                            double nei_self_factor =
                                    _Y[g][nei_index] - _mu[g][_Z[nei_index]] * _cell_size[nei_index] -
                                    _gamma[g][0] *
                                    (_mu[g][k] * _cell_size[MB_idx] * _nei_exp_dist[MB_idx][j] + nei_nei_sum);
                            nei_self_sum += nei_self_factor * nei_self_factor /
                                            (2.0 * _sigma_sq[g]);
                        }
                    } else {
                        break;
                    }
                }

                prob[k] = log(_pi[k]) - self_sum - nei_self_sum;
            }

            prob_max = my_max(prob, _K);
            prob_norm = 0;
#pragma omp simd
            for (int k_t = 0; k_t < _K; ++k_t) {
                prob[k_t] = exp(prob[k_t] - prob_max);
                prob_norm += prob[k_t];
            }

#pragma omp simd
            for (int k_t = 0; k_t < _K; ++k_t) {
                prob[k_t] = prob[k_t] / prob_norm;
            }


            _Z[MB_idx] = rand_cate(prob, rng[omp_get_thread_num()]);
            delete[] prob;
        }

//        cout << "Now it's " << m + 1 << " blanket!" << endl;
//        double end = omp_get_wtime();
//        double parallel_time = end - start;
//        cout << "Time is " << parallel_time << "!" << endl;
    }//End for-i loop

    for (int i = 0; i < _tail_length; ++i) {
        int MB_idx = _tail_list[i];
        double prob_max, prob_norm;
        auto *prob = new double[_K]();

#pragma omp parallel for shared(prob)
        for (int k = 0; k < _K; ++k) {
            double self_sum = 0;
            // self sum
#pragma omp parallel for reduction(+:self_sum)
            for (int g = 0; g < _G; ++g) {
                double nei_sum = 0;


                for (int j = 0; j < _N_nei; ++j) {
                    int nei_idx = _neighbor_list[MB_idx][j];
                    if (nei_idx >= 0) {
                        nei_sum += _mu[g][_Z[nei_idx]] * _cell_size[nei_idx];
                    } else {
                        break;
                    }
                }
                double self_factor = _Y[g][MB_idx] - _mu[g][k] * _cell_size[MB_idx] - _gamma[g][0] * nei_sum;
                self_sum += self_factor * self_factor / (2.0 * _sigma_sq[g]);
            }

            //neighbor sum
            double nei_self_sum = 0;

//#pragma omp parallel for reduction(+:nei_self_sum)
            for (int j = 0; j < _N_nei; ++j) {
                int nei_index = _neighbor_list[MB_idx][j];

                if (nei_index >= 0) {

                    for (int g = 0; g < _G; ++g) {
                        double nei_nei_sum = 0;

                        for (int r = 0; r < _N_nei; ++r) {
                            int nei_nei_idx = _neighbor_list[nei_index][r];
                            if ((nei_nei_idx >= 0) && (nei_nei_idx != MB_idx)) {
                                nei_nei_sum += _mu[g][_Z[nei_nei_idx]] * _cell_size[nei_nei_idx];
                            } else if (nei_nei_idx < 0) {
                                break;
                            }
                        }
                        double nei_self_factor = _Y[g][nei_index] - _mu[g][_Z[nei_index]] * _cell_size[nei_index] -
                                                 _gamma[g][0] * (_mu[g][k] * _cell_size[MB_idx] + nei_nei_sum);
                        nei_self_sum += nei_self_factor * nei_self_factor / (2.0 * _sigma_sq[g]);
                    }
                } else {
                    break;
                }
            }

            prob[k] = log(_pi[k]) - self_sum - nei_self_sum;
        }

        prob_max = my_max(prob, _K);
        prob_norm = 0;
#pragma omp simd
        for (int k_t = 0; k_t < _K; ++k_t) {
            prob[k_t] = exp(prob[k_t] - prob_max);
            prob_norm += prob[k_t];
        }

#pragma omp simd
        for (int k_t = 0; k_t < _K; ++k_t) {
            prob[k_t] = prob[k_t] / prob_norm;
        }

        // Sample from categorical distribution deterministically
//        double u = uniform_dist(rng[omp_get_thread_num()]);
//        double cumulative_prob = 0.0;
//        for (int k_t = 0; k_t < _K; ++k_t) {
//            cumulative_prob += prob[k_t];
//            if (u < cumulative_prob) {
//                _Z[MB_idx] = k_t;
//                break;
//            }
//        }

        _Z[MB_idx] = rand_cate(prob, rng[omp_get_thread_num()]);

        delete[] prob;
    }
//    cout << "Now it's the tail!" << endl;
//    double end = omp_get_wtime();
//    double parallel_time = end - start;
//    cout << "Time is " << parallel_time << "!" << endl;

    double end = omp_get_wtime();
    double parallel_time = end - start;
    cout << "Time for updating Z is " << parallel_time << "!" << endl;
}

void _update_pi(int _N, int _K,
                double _alpha,
                int *_Z,
                mt19937 &rng,
                double *_pi) {
    auto *count_Z = new double[_K]();
    for (int k = 0; k < _K; k++) {
        count_Z[k] = _alpha;
    }

    for (int i = 0; i < _N; ++i) {
        count_Z[_Z[i]]++;
    }

    rand_Dir(count_Z, _K, rng, _pi);
    delete[] count_Z;
}

void _contingency_table(int _N, int _K, int *_Z_last, int *_Z_current, double **_con_table) {
    for (int k = 0; k < _K; ++k) {
        for (int l = 0; l < _K; ++l) {
            _con_table[k][l] = 0;
        }
    }

    for (int i = 0; i < _N; ++i) {
        _con_table[_Z_last[i]][_Z_current[i]]++;
    }
}

#endif //MERIFISH_MCMC_V6_1_H

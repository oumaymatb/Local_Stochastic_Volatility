#ifndef HESTONLOCALVOLATILITY_H
#define HESTONLOCALVOLATILITY_H
#include <iostream>
#include <vector> 
#include <map>
#include<algorithm>
#include "HestonModel.h"

using Vector = std::vector<double>;
using Matrix = std::vector<Vector>;
using Pair = std::pair<double, double>;
using PairVector = std::vector<Pair>;
using namespace std;

#pragma once

class HestonLocalVolatility
{
public:
    HestonLocalVolatility(const HestonModel& model, const size_t& discretization_size, const size_t& number_simulation, const double& maturity, 
                          const Vector& maturities, const Vector& strikes, const Matrix& implied_vols, const double& psi_switch = 1.5, 
                          const double& epsilon_maturity = 0.001, const double& epsilon_strike = 0.001);
    ~HestonLocalVolatility();

    // double local_vol_value(const double& time, const double& spot);

    double price_call(const double& strike);

private:
    HestonModel _model;
    size_t _discretization_size;
    size_t _number_simulation;
    double _maturity;
    double _psi_switch;
    Vector _maturities;
    Vector _strikes;
    Matrix _implied_vols;
    double _epsilon_maturity;
	double _epsilon_strike;
    

    Matrix generate_spots();
    double QEScheme(const double& psi, const double& mean);
    void compute_conditional_expectation(const Vector& last_log_spot, const Vector& bins_time_step, const map<double, double>& spot_to_variance, map<double, double>& spot_to_bin, map<double, double>& bin_to_conditional_expectation);

};

#endif
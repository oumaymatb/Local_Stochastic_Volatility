#include <random>
#include <math.h>
#include <algorithm>
#include <map> 
#include <iostream>
#include <iterator>
#include "HestonLocalVolatility.h"
#include "Chisquare.h"
#include "ImpliedVolatilitySurface.h"
#include "DupireLocalVolatilitySurface.h"

using namespace std;

HestonLocalVolatility::HestonLocalVolatility(const HestonModel& model, const size_t& discretization_size, const size_t& number_simulation, const double& maturity, 
                          const Vector& maturities, const Vector& strikes, const Matrix& implied_vols, const double& psi_switch, const double& epsilon_maturity, const double& epsilon_strike):
                                            _model(model), _discretization_size(discretization_size), _number_simulation(number_simulation), _maturity(maturity), 
                                            _maturities(maturities), _strikes(strikes), _implied_vols(implied_vols), _psi_switch(psi_switch),
                                            _epsilon_maturity(_epsilon_maturity), _epsilon_strike(epsilon_strike)
{}

HestonLocalVolatility::~HestonLocalVolatility()
{}

double HestonLocalVolatility::QEScheme(const double& psi, const double& mean)
{
    random_device rd;  // Will be used to obtain a seed for the random number engine
	mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
	normal_distribution<double> dis(0.0, 1.0);
	uniform_real_distribution<> dis_uni(0.0, 1.0);
	if (psi <= _psi_switch)
			{
				// cout << "I'm up "<<endl;
				double inverse_psi = 1/(double)psi;
				double b_square = 2*inverse_psi - 1 + sqrt(2*inverse_psi*(2*inverse_psi - 1));
				double a = mean/(double)(1 + b_square);
				
				double variance_sim_step = a * pow((sqrt(b_square) + dis(gen)), 2.0);

				return variance_sim_step;
			}
	else
			{
				// cout << "I'm down " <<endl;
				double p = (psi - 1)/(double)(psi + 1);
				// cout << "value of p " << p << endl;
				double beta = 2/(double)(mean*(psi + 1));
				double U_v = dis_uni(gen);
				if (U_v>=0 && U_v<=p)
				{
					double variance_sim_step  = 0.0;
					return variance_sim_step;
				}
				else
				{
					double variance_sim_step  = ( log(1 - p) - log(1 - U_v)  )/(double)beta;
					return variance_sim_step;
				}
            }
}

void HestonLocalVolatility::compute_conditional_expectation(const Vector& last_log_spot, const Vector& bins_time_step, const map<double, double>& spot_to_variance, map<double, double>& spot_to_bin, map<double, double>& bin_to_conditional_expectation)
{
    size_t bin_idx = 0;
    double current_bin = bins_time_step[bin_idx];
    double N_k = 0;
    double sum_variance_bin = 0;
	double value_coniditon_expectation;

	for(auto bin:bins_time_step)
	{
		for(auto spot:last_log_spot)
		{
			if (spot <= bin)
			{
				N_k += 1;
				sum_variance_bin += spot_to_variance.find(spot)->second;
				spot_to_bin.insert(pair<double,double> (spot, bin));
			}
		}
		value_coniditon_expectation = sum_variance_bin/(double)N_k;
		bin_to_conditional_expectation.insert(pair<double,double> (bin, value_coniditon_expectation));
	}
}

Matrix HestonLocalVolatility::generate_spots()
{
    // intialise variances matrix with initial variance
    Matrix variances;

    // initialise log spots matrix with initial log spot
    Matrix log_spots;

    ImpliedVolatilitySurface implied_vol(_maturities, _strikes, _implied_vols, _model.drift());
    DupireLocalVolatilitySurface local_vol(implied_vol, _epsilon_maturity, _epsilon_strike, _model.initial_spot());
    double local_vol_value = local_vol.local_volatility(0.0,_model.initial_spot());
    local_vol_value = local_vol_value/(sqrt(_model.initial_variance()));

    // variables for condition expectation computation
	Vector conditional_expectations;
	map<double, double> spot_to_bin;
	map<double,double> bin_to_conditional_expectation;
	map<double, double> spot_to_variance;

	// initialise everything
	double delta = _maturity/(double)(_discretization_size);
	Vector initialise_variance;
	Vector initialise_log_spot;

	double second, variance_value, mean, s_square, psi, rho_1, c_1, first, r_delta, sigma_variance, log_spot, lsv, sigma_delta, sigma_normal ;

 
	for(size_t nbr=0; nbr<_number_simulation; nbr++)
	{
		second = _model.initial_variance(); //v_i_j
		mean = _model.mean_reversion_level() + (second - _model.mean_reversion_level())*exp(-_model.mean_reversion_speed()*delta);
		s_square = ( second*_model.vol_of_vol()*_model.vol_of_vol()*exp(-_model.mean_reversion_speed()*delta)*(1 - exp(-_model.mean_reversion_speed()*delta)))/_model.mean_reversion_speed();
		s_square += (_model.mean_reversion_level()*_model.vol_of_vol()*_model.vol_of_vol()*pow( (1 - exp(-_model.mean_reversion_speed()*delta)),2.0 ))/(2.0*_model.mean_reversion_speed());
		psi = s_square/(mean*mean);

		random_device rd;  // Will be used to obtain a seed for the random number engine
		mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
		normal_distribution<double> dis(0.0, 1.0);

		variance_value = QEScheme(psi, mean);
		// cout << "variance_value at " <<nbr << " " << variance_value <<endl;
		initialise_variance.push_back(variance_value);

		rho_1 = sqrt(1 - pow(_model.correlation(),2.0));
		c_1 = _model.mean_reversion_speed()*delta - 1;
		first = _model.initial_spot();
		r_delta = _model.drift()*delta;
		sigma_variance = (_model.correlation()*(variance_value + second*c_1 - _model.mean_reversion_speed()*_model.mean_reversion_level()*delta))/(double)(_model.vol_of_vol());
		lsv = local_vol_value;
        sigma_delta = 0.5*lsv*lsv*delta;
        sigma_normal = rho_1*lsv*sqrt(second*delta);
        log_spot = first + r_delta - sigma_delta*second + sigma_variance + sigma_normal*dis(gen);
    	initialise_log_spot.push_back(log_spot);

		spot_to_variance.insert ( pair<double, double>(log_spot, variance_value));
	}
	variances.push_back(initialise_variance);
	log_spots.push_back(initialise_log_spot);

    for(size_t time_step = 1; time_step<_discretization_size-1; time_step++)
    {
        // to store spot and variance values
        Vector variance_time_step;
        Vector log_spot_time_step;
		Vector last_variance = variances.back();
		Vector last_log_spot = log_spots.back();

        // STEP 1:  Generate S_t_i and V_t_i
        for(size_t sim_step =0; sim_step<_number_simulation; sim_step++)
        {	
			// Generate bins
			sort(last_log_spot.begin(), last_log_spot.end());

            size_t bins_size = 5;
            Vector bins_time_step;

            for(size_t bins_step = 0; bins_step<bins_size; bins_step++)
            {
                double value_bin = last_log_spot[0] + (bins_step*(last_log_spot[_number_simulation-1] - last_log_spot[0]))/(double)(bins_size);
				// cout << "value_bin at bin_step for time step " << time_step << " " << bins_step << " " << value_bin << endl;
				bins_time_step.push_back(value_bin);
            }

			//Compute the conditional expectation
			compute_conditional_expectation(last_log_spot, bins_time_step, spot_to_variance, spot_to_bin, bin_to_conditional_expectation);

			// geneate next variance
			second = last_variance[sim_step]; //v_i_j
			mean = _model.mean_reversion_level() + (second - _model.mean_reversion_level())*exp(-_model.mean_reversion_speed()*delta);
			s_square = ( second*_model.vol_of_vol()*_model.vol_of_vol()*exp(-_model.mean_reversion_speed()*delta)*(1 - exp(-_model.mean_reversion_speed()*delta)))/_model.mean_reversion_speed();
			s_square += (_model.mean_reversion_level()*_model.vol_of_vol()*_model.vol_of_vol()*pow( (1 - exp(-_model.mean_reversion_speed()*delta)),2.0 ))/(2.0*_model.mean_reversion_speed());
			psi = s_square/(mean*mean);

			random_device rd;  // Will be used to obtain a seed for the random number engine
            mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
            normal_distribution<double> dis(0.0, 1.0);

			variance_value = QEScheme(psi, mean);
			variance_time_step.push_back(variance_value);

			rho_1 = sqrt(1 - pow(_model.correlation(),2.0));
            c_1 = _model.mean_reversion_speed()*delta - 1;
            first = last_log_spot[sim_step];
            r_delta = _model.drift()*delta;
            sigma_variance = (_model.correlation()*(variance_value + second*c_1 - _model.mean_reversion_speed()*_model.mean_reversion_level()*delta))/(double)(_model.vol_of_vol());
			double bin, condition_expectation_bin, lsv, sigma_delta, sigma_normal, log_spot_sim_step;

			if(isnan(first) or isinf(first))
			{
				continue;
			}
			else
			{
				local_vol_value = local_vol.local_volatility(time_step*delta, first/100.0);
				bin = spot_to_bin[first];
				condition_expectation_bin = bin_to_conditional_expectation[bin];
				lsv = local_vol_value/sqrt(condition_expectation_bin);
				sigma_delta = 0.5*lsv*lsv*delta;
				sigma_normal = rho_1*lsv*sqrt(second*delta);

				log_spot_sim_step = first + r_delta - sigma_delta*second + sigma_variance + sigma_normal*dis(gen);

				log_spot_time_step.push_back(log_spot_sim_step);

				// fill the map ==> assign a variance to each spot
				spot_to_variance.insert ( pair<double, double>(log_spot_sim_step,variance_value) );
			}
        }
		log_spots.push_back(log_spot_time_step);
        variances.push_back(variance_time_step);
    }
	return log_spots;
}

double HestonLocalVolatility::price_call(const double& strike) 
{
    Matrix paths_spots = generate_spots();
    Vector maturity_spot = paths_spots.back();
    double sum = 0.0;
    for (auto spot:maturity_spot)
    {
        sum += max(spot - strike, 0.0);
    }

    return (exp(-_model.drift()*_maturity)*sum)/(double)(_number_simulation);
}

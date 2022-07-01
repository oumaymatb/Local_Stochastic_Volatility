#include "DupireLocalVolatilitySurface.h"
#include <math.h>
#include <iostream>

using namespace std;

DupireLocalVolatilitySurface::DupireLocalVolatilitySurface(const ImpliedVolatilitySurface & implied_vol_surface, const double& eps_mat, const double& eps_strike, const double& init_spot)
	: _implied_volatility_surface(implied_vol_surface), _epsilon_maturity(eps_mat), _epsilon_strike(eps_strike), _initial_spot(init_spot)
{
}

// return sigma_D(T,K)
double DupireLocalVolatilitySurface::local_volatility(const double & maturity, const double & strike) const
{
	// sigma*(T,K)
	double sigma = _implied_volatility_surface.implied_volatility(maturity, strike);
	// cout << "sigma " << sigma << endl;

	double d_sigma_dT = first_order_derivative_impliedvol_maturity(maturity, strike);
	// cout << "d_sigma_dT " << d_sigma_dT << endl;
	double d_sigma_dK = first_order_derivative_impliedvol_strike(maturity, strike);
	// cout << "d_sigma_dK " << d_sigma_dK << endl;
	double d2_sigma_dK2 = second_order_derivative_impliedvol_strike(maturity, strike);
	// cout << "d2_sigma_dK2 " << d2_sigma_dK2 << endl;
	double risk_free_rate = _implied_volatility_surface.risk_free_rate();

	double d1 = (log(_initial_spot / strike) + (risk_free_rate + 0.5 * sigma * sigma)*maturity) / (double) (sigma* sqrt(maturity));
	// cout << "d1 " << d1 << endl;
	double d2 = d1 - sigma * sqrt(maturity);
	// cout << "d2 " << d2 << endl;
	
	double numerator = 1. + (2.*maturity / sigma) * (d_sigma_dT + risk_free_rate * d_sigma_dK);
	// cout << "numerator " << numerator << endl;

	double k_s_sqrtT = strike * sigma * sqrt(maturity);
	// cout << "k_s_sqrtT " << k_s_sqrtT << endl;
	double k_ds_dK_sqrtT = strike * d_sigma_dK * sqrt(maturity);
	// cout << "k_ds_dK_sqrtT " << k_ds_dK_sqrtT << endl;
	double k_d2s_dK2_sqrtT = strike * d2_sigma_dK2 * sqrt(maturity);
	// cout << "k_d2s_dK2_sqrtT " << k_d2s_dK2_sqrtT << endl;

	double denominator = 1. + 2. * d1 * k_ds_dK_sqrtT + d1*d2*k_ds_dK_sqrtT*k_ds_dK_sqrtT + k_d2s_dK2_sqrtT * k_s_sqrtT;
	// cout << "denominator " << denominator << endl;

	double dupire_variance = sigma * sigma * numerator /(double)denominator;
	// cout << "dupire_variance " << dupire_variance << endl;
	double dupire_volatility = sqrt(dupire_variance);
	return dupire_volatility;
}

double DupireLocalVolatilitySurface::first_order_derivative_impliedvol_maturity(const double & maturity, const double & strike) const
{
	double sigma_up = _implied_volatility_surface.implied_volatility(maturity + _epsilon_maturity, strike);
	double sigma_down = _implied_volatility_surface.implied_volatility(maturity - _epsilon_maturity, strike);
	double deriv = 0.5 * (sigma_up - sigma_down) /(double) _epsilon_maturity;
	return deriv;
}

double DupireLocalVolatilitySurface::first_order_derivative_impliedvol_strike(const double & maturity, const double & strike) const
{
	double sigma_up = _implied_volatility_surface.implied_volatility(maturity, strike + _epsilon_strike);
	double sigma_down = _implied_volatility_surface.implied_volatility(maturity, strike - _epsilon_strike);
	double deriv = 0.5 * (sigma_up - sigma_down) /(double)_epsilon_strike;
	return deriv;
}

double DupireLocalVolatilitySurface::second_order_derivative_impliedvol_strike(const double & maturity, const double & strike) const
{
	double sigma_up = _implied_volatility_surface.implied_volatility(maturity, strike + _epsilon_strike);
	double sigma_base = _implied_volatility_surface.implied_volatility(maturity, strike);
	double sigma_down = _implied_volatility_surface.implied_volatility(maturity, strike - _epsilon_strike);
	double deriv = (sigma_up - 2.*sigma_base + sigma_down) /(double) (_epsilon_strike*_epsilon_strike);
	return deriv;
}
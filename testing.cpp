#include <iostream>
#include <typeinfo>
#include "DupireLocalVolatilitySurface.h"
#include "HestonPricer.h"
#include "RootSearch.hpp"
#include "HestonLocalVolatility.h"

using namespace std;

int main()
{

	//FILE* stream;
	//freopen_s(&stream, "Heston_call_prices.csv", "w", stdout);

	/* Risk Free Rate and init spot*/

	double risk_free_rate = 0.05;
	double init_spot = 100.;
	double init_variance = 0.1;
	double vol_of_vol = 0.95;
	double vol_speed = 1.05;
	double vol_level = 0.0855;
	double correl = -0.315;

	/* Implied Volatility object */
	Vector market_maturities = { 0.25, 0.5, 0.75, 1., 2., 3., 4., 5. };
	Vector market_strikes = { 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5};
		Matrix market_implied_vols = {
		{0.39, 0.31, 0.24, 0.22, 0.16, 0.19, 0.23, 0.29, 0.38},
		{0.44, 0.36, 0.27, 0.21, 0.17, 0.21, 0.27, 0.35, 0.40},
		{0.45, 0.30, 0.25, 0.21, 0.18, 0.22, 0.29, 0.37, 0.45},
		{0.48, 0.42, 0.34, 0.28, 0.20, 0.26, 0.31, 0.42, 0.50},
		{0.52, 0.43, 0.34, 0.26, 0.21, 0.27, 0.38, 0.45, 0.55},
		{0.54, 0.46, 0.34, 0.27, 0.23, 0.28, 0.36, 0.49, 0.58},
		{0.57, 0.50, 0.46, 0.35, 0.25, 0.32, 0.45, 0.54, 0.60},
		{0.60, 0.52, 0.41, 0.31, 0.26, 0.34, 0.40, 0.55, 0.62}
	};

	HestonModel Model(init_spot, init_variance, risk_free_rate, vol_speed, vol_level, vol_of_vol, correl);
	double maturity = 2.0;
	size_t trapeze_size = 100;

	size_t discretization_size = 10;
	size_t number_simulation = 100;

	double heston_price_call = 0.0;
	double local_heston_price_call = 0.0;
	vector<double> heston_call_prices;

	for (auto strike:market_strikes)
	{
		MonteCarloHestonPricer HestonModel(Model, strike, 2., 10000, 100);
		// HestonLocalVolatility LocalHestonModel(Model, discretization_size, number_simulation, maturity, market_maturities, market_strikes, market_implied_vols);
		heston_price_call = HestonModel.price();
		// local_heston_price_call = LocalHestonModel.price_call(strike);
		cout << "heston call for strike  " << strike <<" is " << heston_price_call << endl;
		for (size_t i= 0; i<6; i++)
		{
			Vector implied_vols = market_implied_vols[i];
			cout << "BS call for strike " << strike << " is " << BS_call(init_spot, strike, risk_free_rate, implied_vols[5],maturity);
		}

		// cout << "Implied heston volatility " << calc_implied_vol(heston_price_call, init_spot, strike*100., risk_free_rate, 2., 0, 1) << endl;
		// cout << "Implied local heston volatility " << calc_implied_vol(heston_price_call, init_spot, strike*100., risk_free_rate, 2., 0, 1) << endl;
	}

	return 0;
}
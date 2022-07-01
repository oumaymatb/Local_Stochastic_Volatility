#include <iostream>
#include <typeinfo>
#include "DupireLocalVolatilitySurface.h"
#include "HestonLocalModel.h"
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
	

	HestonLocalVolatility heston_local(Model, discretization_size, number_simulation, maturity, market_maturities, market_strikes, market_implied_vols);
	double price_heston_local_call = heston_local.price_call(0.7);
	cout << "let's print sth as a test" << endl;
	cout << "Heston local price call is " << price_heston_local_call << endl;

	/*

	double heston_price_call = 0.0;
	vector<double> heston_call_prices;

	for (auto strike:market_strikes)
	{
		MonteCarloHestonPricer heston_price(Model, strike, 2., 10000, 100);
		heston_price_call = heston_price.price();
		cout << "heston call for strike  " << strike <<" is " << heston_price_call << endl;
		heston_call_prices.push_back(heston_price_call);
		cout << "Implied heston volatility " << calc_implied_vol(heston_price_call, init_spot, strike*100., risk_free_rate, 2., 0, 1) << endl;
	}


	for (auto strike:market_strikes)
	{
		AnalyticalHestonPricer heston_pricer(Model, strike, maturity, trapeze_size);
		price_call = heston_pricer.price();
		std::cout <<  price_call << std::endl;
	}*/




	/*
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

	ImpliedVolatilitySurface implied_volatility_surface(market_maturities, market_strikes, market_implied_vols, risk_free_rate);
	Vector new_mats;
	for (int i = 1; i <= 120; ++i)
		new_mats.push_back((double)i * 0.05);

	Vector new_strikes;
	for (int j = 1; j <= 20; ++j)
		new_strikes.push_back((double)j * 10.);

	//Matrix new_implied_vols;
	for (int i = 0; i < new_mats.size(); ++i)
	{
		for (int j = 0; j < new_strikes.size(); ++j)
		{
			double sigma = implied_volatility_surface.implied_volatility(new_mats[i], new_strikes[j]);
		
			std::cout << sigma << ',';
		}
		std::cout << '\n';
	}*/

	return 0;
}
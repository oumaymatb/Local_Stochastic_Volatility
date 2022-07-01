#include "HestonPricer.h"
#include "HestonModel.h"
#include "Complex.h"
#include "Chisquare.h"
#include "TrapezoidalRule.h"
#include <vector>
#include <math.h>
#include <random>
#include <chrono> 
#define _USE_MATH_DEFINES // for C++
# define M_PIl          3.141592653589793238462643383279502884L /* pi */

using namespace std;

HestonPricer::HestonPricer(const HestonModel& model, const double& strike, const double& maturity):
	_model(model),_strike(strike),_maturity(maturity)
{}

AnalyticalHestonPricer::AnalyticalHestonPricer(
	const HestonModel& model, 
	const double& strike, 
	const double& maturity, 
	const size_t& _gauss_legendre_sample_size):
	HestonPricer(model,strike,maturity), _gauss_legendre_sample_size(_gauss_legendre_sample_size)
{}

AnalyticalHestonPricer* AnalyticalHestonPricer::clone() const
{
	return nullptr;
}

Complex AnalyticalHestonPricer::func_u(const size_t& i) const
{
	double u = (i == 1) ? -1. : 1.;
	return Complex(u);
}

Complex AnalyticalHestonPricer::func_y(const size_t& i, const double& omega) const
{
	double y_real = (i == 1) ? -1. : 0.;
	return Complex(y_real,omega);
}

Complex AnalyticalHestonPricer::func_a(const size_t& i, const double& omega) const
{
	Complex kappa(_model.mean_reversion_speed());
	Complex rho(_model.correlation());
	Complex sigma_v(_model.vol_of_vol());
	Complex y = func_y(i, omega);

	return kappa- rho*sigma_v*y;
}

Complex AnalyticalHestonPricer::func_delta(const size_t& i, const double& omega) const
{
	Complex a = func_a(i, omega);
	Complex sigma_v(_model.vol_of_vol());
	Complex u = func_u(i);
	Complex j(0., 1.);
	Complex omega_c(omega);

	return a*a + sigma_v* sigma_v*(u*j*omega_c + omega_c*omega_c);
}

Complex AnalyticalHestonPricer::func_b(const size_t& i, const double& omega) const
{
	return Complex::square_root(func_delta(i,omega));
}

Complex AnalyticalHestonPricer::func_g(const size_t& i, const double& omega) const
{
	Complex a=func_a(i,omega);
	Complex b = func_b(i, omega);
	return (a - b) / (a + b);
}

Complex AnalyticalHestonPricer::func_C(const size_t& i,
					const double& tau,
					const double& omega) const
{	
	Complex tau_c(tau);
	Complex sigma_v(_model.vol_of_vol());
	Complex a = func_a(i, omega);
	Complex b = func_b(i, omega);
	Complex g = func_g(i, omega);
	Complex k = Complex(_model.mean_reversion_speed());
	Complex theta = Complex(_model.mean_reversion_level());

	Complex res = (Complex(1.) - g*Complex::exponential(Complex(-1.)*b*tau_c));
	res = res/(Complex(1.) - g);
	res = Complex(0.,omega*tau*_model.drift()) + 
			k*theta/(sigma_v*sigma_v)*( (a -b)*tau_c - Complex(2.)*Complex::logarithm(res) );

	return res;
}

Complex AnalyticalHestonPricer::func_D(const size_t& i,
									   const double& tau,
									   const double& omega) const
{
	Complex tau_c(tau);
	Complex sigma_v(_model.vol_of_vol());
	Complex a = func_a(i, omega);
	Complex b = func_b(i, omega);
	Complex g = func_g(i, omega);

	Complex res = (a - b) * (Complex(1.) - Complex::exponential(Complex(-1.) * b * tau_c)) / (sigma_v * sigma_v * (Complex(1.) - g * Complex::exponential(Complex(-1.) * b * tau_c)));
	return res;
}

Complex AnalyticalHestonPricer::func_Phi(const size_t& i,
					 const double& tau, 
					 const double& x,
					 const double& V,
					 const double& omega) const
{	
	Complex C = func_C(i,tau,omega);
	Complex D = func_D(i,tau,omega);
	Complex res = C + D*Complex(V) + Complex(0.,omega*x);

	return Complex::exponential(res);
}

double AnalyticalHestonPricer::func_integrand(const size_t& i, const double& tau, const double& x, const double& V, const double& omega) const
{
	Complex phi = func_Phi(i,tau,x,V,omega);
	Complex K = Complex(0, omega*log(_strike));
	Complex res = (phi*Complex::exponential( Complex(-1.)*K ))/Complex(0.,omega);

	return res.real_part();
}

double AnalyticalHestonPricer::func_P(const size_t& i, const double& tau, const double& x, const double& V) const
{
	TrapezoidalRule glq(_gauss_legendre_sample_size);
	function<double(double)> func =
		[this, i, tau, x, V](double omega) { return func_integrand(i, tau, x, V, omega); };
	double integral = glq.integrate(func, 0.0, 1.0);
	return 0.5*( 1. + integral/M_PIl);
}

double AnalyticalHestonPricer::price() const
{
	double price_call = 0.0;
	double P_1 = func_P(1,_maturity,log(_model.initial_spot()),_model.initial_variance());
	double P_2 = func_P(2,_maturity,log(_model.initial_spot()),_model.initial_variance());
	price_call = _model.initial_spot()*P_1 - _strike*exp(-_model.drift()*_maturity)*P_2;

	return price_call;
}

MonteCarloHestonPricer::MonteCarloHestonPricer(const HestonModel& model,
						   const double& strike,
						   const double& maturity,
						   const size_t& number_simulations,
						   const size_t& discretization_size,
						   const int& is_maximum):
						   HestonPricer(model,strike,maturity),
						   _number_simulations(number_simulations),_discretization_size(discretization_size), _is_maximum(is_maximum)
{
	double dt = _maturity / double(_discretization_size);
	double time_step = 0.0;
	while(time_step <= _maturity)
	{	
		_discretisation_times.push_back(time_step);
		time_step += dt;
	}
}

MonteCarloHestonPricer* MonteCarloHestonPricer::clone() const 
{
	return new MonteCarloHestonPricer(*this);
}

PairVector MonteCarloHestonPricer::path_2Factors(const double& psi_switch) const
{
	PairVector path2F;

	// Initialisation
	double X_0 = log(_model.initial_spot());
	double V_0 = _model.initial_variance();
	Pair init_pair(X_0, V_0);
	path2F.push_back(init_pair);

	double kappa = _model.mean_reversion_speed();
	double theta = _model.mean_reversion_level();
	double sigma_v = _model.vol_of_vol();
	double correl = _model.correlation();

	// Iteration

	for(size_t dis_idx = 1; dis_idx<_discretization_size; dis_idx++)
	{
		double log_spot;
		double variance;
		double delta = _maturity / (double)_discretization_size;
		double second = path2F[dis_idx-1].second;
		double first = path2F[dis_idx-1].first;
		if (_is_maximum == 1) 
		{
			random_device rd;  // Will be used to obtain a seed for the random number engine
			mt19937 gen1(rd()); // Standard mersenne_twister_engine seeded with rd()
			mt19937 gen2(rd());
			normal_distribution<double> dis(0.0, 1.0);

			double n_1 = dis(gen1);
			double n_2 = dis(gen2);

			double z_v = _model.correlation()*n_1 + sqrt(1 - _model.correlation()*_model.correlation())*n_2;

			variance = second + _model.mean_reversion_speed()*(_model.mean_reversion_level()-max(0.0,second))*delta +
						_model.vol_of_vol()*sqrt(max(0.0,second))*sqrt(delta)*z_v;
			log_spot = first - 0.5*max(0.0,second)*delta + sqrt(max(0.0,second))*n_1*sqrt(delta);

			Pair current_pair(log_spot, variance);
			path2F.push_back(current_pair);

		}
		else 
		{
			double mean = _model.mean_reversion_level() + (second - _model.mean_reversion_level())*exp(-_model.mean_reversion_speed()*delta);
			double s_square = ( second*_model.vol_of_vol()*_model.vol_of_vol()*exp(-_model.mean_reversion_speed()*delta)*(1 - exp(-_model.mean_reversion_speed()*delta)))/_model.mean_reversion_speed();
			s_square += (_model.mean_reversion_level()*_model.vol_of_vol()*_model.vol_of_vol()*pow( (1 - exp(-_model.mean_reversion_speed()*delta)),2.0 ))/(2.0*_model.mean_reversion_speed());
			double psi = s_square/(mean*mean);

			random_device rd;  // Will be used to obtain a seed for the random number engine
			mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
			normal_distribution<double> dis(0.0, 1.0);
			uniform_real_distribution<> dis_uni(0.0,1.0);

			double corr_eps = _model.correlation()/(double)_model.vol_of_vol();
			double K_0 = - _model.mean_reversion_speed()*corr_eps*_model.mean_reversion_level()*delta;
			double K_1 = delta*(_model.mean_reversion_speed()*corr_eps - 0.5) - corr_eps;
			double K_2 = corr_eps;
			double K_3 = delta*(1 - _model.correlation()*_model.correlation());
			double A = K_2;

			if (psi <= psi_switch)
			{
				// cout << "I'm up "<<endl;
				double inverse_psi = 1/(double)psi;
				double b_square = 2*inverse_psi - 1 + sqrt(2*inverse_psi*(2*inverse_psi - 1));
				double a = mean/(double)(1 + b_square);
				
				variance = a * pow((sqrt(b_square) + dis(gen)), 2.0);

				// gamma_1 = 1 and gamma_2 = 0 == Euler
				// ==> K_4 = 0
				double M = (exp((A*b_square*a)/(1 - 2*A*a)))/(sqrt(1 - 2*A*a));
				double K_0_star = - log(M) - (K_1 + 0.5*delta);

				log_spot = first + K_0_star + K_1*second + K_2*variance + sqrt(K_3*second)*dis(gen);
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
					variance = 0.0;
				}
				else
				{
					variance = ( log(1 - p) - log(1 - U_v)  )/(double)beta;
				}

				double M = (beta*(1 - p))/(double)(beta - A);
				double K_0_star = - log(M) - (K_1 + 0.5*delta);

				log_spot = first + K_0_star + K_1*second + K_2*variance + sqrt(K_3*second)*dis(gen);
			}

			Pair current_pair(log_spot, variance);
			// cout << "log_spot is " << log_spot << " and variance is " << variance <<endl;
			path2F.push_back(current_pair);
		}
	}

	return path2F;

}

double MonteCarloHestonPricer::price() const
{
	double price = 0.0;
	for (size_t sim_idx = 0; sim_idx<_number_simulations; sim_idx++)
	{
		PairVector path = path_2Factors();
		Pair last_pair = path.back();
		double spot_matu = last_pair.first;
		// cout << "spot at maturity is " << spot_matu <<endl;
		price += max(0.0, spot_matu - _strike);
		// cout << "price cumul is " << price << endl;
	}
	price = exp(-_model.drift()*_maturity) * (price/(double)_number_simulations);
	// cout << "price after is " << price << endl;
	return price;
}

double MonteCarloHestonPricer::PsiInverse(const double& u, const double& p, const double& beta) const
{
	if((u>= 0) && (u<=p))
	{
		return 0.0;
	}
	else
	{
		return log((1 - p)/(1 - u))/(double)beta;
	}
}





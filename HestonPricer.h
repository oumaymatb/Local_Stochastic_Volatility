#ifndef HESTONPRICER_H
#define HESTONPRICER_H

#include "HestonModel.h"
#include "Complex.h"
#include <vector>

using Vector = std::vector<double>;
using Pair = std::pair<double, double>;
using PairVector = std::vector<Pair>;

class HestonPricer
{
public:
	HestonPricer(const HestonModel& model, const double& strike, const double& maturity);
	virtual HestonPricer* clone() const = 0;
	virtual double price() const = 0;
protected:
	HestonModel _model;
	double _strike;
	double _maturity;
};


class AnalyticalHestonPricer : public HestonPricer
{
public:
	AnalyticalHestonPricer(const HestonModel& model,
		const double& strike, const double& maturity,
		const size_t& gauss_legendre_sample_size);

	AnalyticalHestonPricer* clone() const override;

	double price() const override;

protected:
	Complex func_u(const size_t& i) const;
	Complex func_y(const size_t& i, const double& omega) const;
	Complex func_a(const size_t& i, const double& omega) const;
	Complex func_delta(const size_t& i, const double& omega) const;
	Complex func_b(const size_t& i, const double& omega) const;
	Complex func_g(const size_t& i, const double& omega) const;

	Complex func_C(const size_t& i,
					const double& tau,
					const double& omega) const;

	Complex func_D(const size_t& i,
					const double& tau,
					const double& omega) const;

	Complex func_Phi(const size_t& i,
					 const double& tau, 
					 const double& x,
					 const double& V,
					 const double& omega) const;

	double func_integrand(const size_t& i, const double& tau, const double& x, const double& V, const double& omega) const;
	
	double func_P(const size_t& i, const double& tau, const double& x, const double& V) const;

	size_t _gauss_legendre_sample_size;

};

class MonteCarloHestonPricer : public HestonPricer
{
public:
	MonteCarloHestonPricer(const HestonModel& model,
						   const double& strike,
						   const double& maturity,
						   const size_t& number_simulations,
						   const size_t& discretization_size,
						   const int& is_maximum = 0);

	MonteCarloHestonPricer* clone() const override;
	double price() const override;

protected:
	size_t _number_simulations;
	size_t _discretization_size;
	int _is_maximum;

	PairVector path_2Factors(const double& psi_switch = 1.5) const; //if is_maximum = 1 = maximum else QE
	Vector _discretisation_times;
	double PsiInverse(const double& u, const double& p, const double& beta) const;

};

#endif
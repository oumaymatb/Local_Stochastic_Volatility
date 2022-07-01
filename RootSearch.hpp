#include <cmath>

using namespace std;

double toy_f(double x);

double normal_cdf(double x);

double d_j(const int& j, const double& S, const double& K, const double& r, const double& sigma, const double& T);

double BS_call(const double& S, const double& K, const double& r, const double& sigma, const double& T);


double calc_implied_vol(const double& MtM_call, const double& S, const double& K, const double& r, const double& T,
                         double lower_b,  double upper_b, const double& precision=pow(10,-5));

double root_search(const double& f(double),  const double& target,  double lower_b,  double&upper_b, const double& precision=pow(10,-1));
#include <iostream>
#include <cmath>
#include <vector>
#include <array>

using namespace std;

double toy_f(double x){
    return x - 1;
}

double normal_cdf(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

double d_j(const int& j, const double& S, const double& K, const double& r, const double& sigma, const double& T) {
    return (log(S/K) + (r + (pow(-1,j-1))*0.5*sigma*sigma)*T)/(sigma*(pow(T,0.5)));
}

double BS_call(const double& S, const double& K, const double& r, const double& sigma, const double& T) {
    return S * normal_cdf(d_j(1, S, K, r, sigma, T))-K*exp(-r*T) * normal_cdf(d_j(2, S, K, r, sigma, T));
}


double calc_implied_vol(const double& MtM_call, const double& S, const double& K, const double& r, const double& T,
                        double lower_b, double upper_b, const double& precision=pow(10,-1)){
    double mid = 0;
    double BS_price = BS_call(S, K, r, mid, T);
    do{
        mid = (lower_b + upper_b)/2;
        BS_price = BS_call(S, K, r, mid, T);
        if(BS_price - MtM_call> 0){
            upper_b = mid;
        }else{
            lower_b = mid;
        }
    }while(!((fabs(BS_price - MtM_call)) < precision));

    return mid;
}

// Solves f(x) = target, where f is a non-decreasing function
double root_search(double f(double), double target, double lower_b, double upper_b, double precision=pow(10,-1)){
    double mid = 0;
    do{
        mid = (lower_b + upper_b)/2;
        if(f(mid) - target> 0){
            upper_b = mid;
        }else{
            lower_b = mid;
        }
    }while(!(fabs(f(mid)-target)) < precision);

    return mid;
}



int main()
{
    double risk_free_rate = 0.05;
	double init_spot = 100.;
	double init_variance = 0.1;
	double vol_of_vol = 0.95;
	double vol_speed = 1.05;
	double vol_level = 0.0855;
	double correl = -0.315;
    double maturity = 2.;
    double strike = 80.;
    double implied_vol = 0.21;
	cout << "BS call for strike " << strike << " is " << BS_call(init_spot, strike, risk_free_rate, implied_vol,maturity);

}

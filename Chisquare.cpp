#include "Chisquare.h"
#include <random>
#include <chrono> 
#include <math.h>

using namespace std;

Chisquare::Chisquare(const double& degrees_freedom, const double& non_central_param):
_degrees_freedom(degrees_freedom),_non_central_param(non_central_param)
{
    for (size_t i =0; i<degrees_freedom; i++)
    {   
        default_random_engine generator1;
		generator1.seed(std::chrono::system_clock::now().time_since_epoch().count());
        if (i==0){
            normal_distribution<double> distribution(0.0,sqrt(non_central_param));
            _value += pow(distribution(generator1),2.0);}
        else{
            normal_distribution<double> distribution(0.0,1.0);
            _value += pow(distribution(generator1),2.0);
        }
    }
}

Chisquare::~Chisquare()
{

}

double Chisquare::generate() const
{
    return _value;
}
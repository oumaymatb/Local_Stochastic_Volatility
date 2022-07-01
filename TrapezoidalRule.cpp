#include "TrapezoidalRule.h"

TrapezoidalRule::TrapezoidalRule(const size_t& number_sample_points):_number_sample_points(number_sample_points)
{}

TrapezoidalRule::~TrapezoidalRule()
{}

double TrapezoidalRule::integrate(std::function<double(double)> func, const double& a, const double& b) const
{
    double integrale = 0.0;
    double delta_x = (b - a)/(double) _number_sample_points;
    double x = a;
    for(size_t idx = 0; idx<_number_sample_points; idx++)
    {
        double f_1 = func(x);
        x = x + delta_x;
        double f_2 = func(x);
        integrale += ((f_1 + f_2)*delta_x)*0.5;
    }
    return integrale;
}

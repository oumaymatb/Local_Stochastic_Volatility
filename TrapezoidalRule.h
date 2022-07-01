#ifndef TRAPEZOIDALRULE_H
#define TRAPEZOIDALRULE_H

#include <vector>
#include <functional>
#include<iostream>
#include <math.h>

#pragma once

class TrapezoidalRule
{
public:
    TrapezoidalRule(const size_t& number_sample_points);
    ~TrapezoidalRule();
    
    // integrate
    double integrate(std::function<double(double)> func, const double& a, const double& b) const;

private:
    size_t _number_sample_points;
};

#endif
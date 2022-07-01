#ifndef CHISQUARE_H
#define CHISQUARE_H

#pragma once

class Chisquare
{
public:
    Chisquare(const double& degrees_freedom, const double& non_central_param);
    ~Chisquare();
    double generate() const;

private:
    double _degrees_freedom;
    double _non_central_param;
    double _value;
};

#endif
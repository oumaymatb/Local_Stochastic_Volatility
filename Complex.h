#ifndef COMPLEX_H
#define COMPLEX_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cerrno>
#include <cstring>
#include <cfenv>

using Vector = std::vector<double>;
using Matrix = std::vector<Vector>;

class Complex
{
public:
	// Constructors
	// default constructor
	Complex();
	// constructor with parameters
	Complex(const double& real, const double& imaginary = 0.0);
	// Copy constructor
	Complex(const Complex& complex);
	// Destructor
	virtual ~Complex() = default;

	// getters
	double real_part() const;
	double imaginary_part() const;

	// Operators
	Complex& operator=(const Complex&  complex);
	Complex operator+(const Complex&  complex) const;
	Complex operator-(const Complex&  complex) const;
	Complex operator*(const Complex&  complex) const;
	Complex operator/(const Complex&  complex) const;

	// Methods
	static Complex exponential(const Complex& complex);
	static Complex logarithm(const Complex& complex);
	static Complex square_root(const Complex& complex);


private:
	// Check method
	bool isZero() const;
	
	double module() const;
	double argument() const;

	double _real; 
	double _imaginary;

};

#endif 
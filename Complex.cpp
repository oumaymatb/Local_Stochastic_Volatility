#include <iostream>
#include <stdexcept>
#include <limits>
#include "Complex.h"

using namespace std;

//////////////////////////////////  Complex //////////////////////////////////////////////

Complex::Complex():
	_real(0), _imaginary(0)
{}

Complex::Complex(const double& real, const double& imaginary):
	_real(real), _imaginary(imaginary)
{}

Complex::Complex(const Complex& complex):
	_real(complex._real),_imaginary(complex._imaginary)
{}

double Complex::real_part() const
{
	return _real;
}

double Complex::imaginary_part() const
{
	return _imaginary;
}

Complex& Complex::operator=(const Complex&  complex) 
{
	if (this != &complex)
	{
		_real = complex._real;
		_imaginary = complex._imaginary;
	}
	return *this;
}

Complex Complex::operator+(const Complex& complex) const
{
	double sum_real = _real+ complex._real ;
	double sum_imaginary = _imaginary + complex._imaginary; 
	return Complex(sum_real, sum_imaginary);
}

Complex Complex::operator-(const Complex& complex) const
{
	double diff_real = _real - complex._real;
	double diff_imaginary = _imaginary - complex._imaginary;
	return Complex(diff_real, diff_imaginary);
}

Complex Complex::operator*(const Complex& complex) const
{
	double multip_real = _real* complex._real - _imaginary * complex._imaginary;
	double multip_imaginary = _real* complex._imaginary + _imaginary * complex._real;
	return Complex(multip_real, multip_imaginary);
}

Complex Complex::operator/(const Complex& complex) const
{
	if (complex.isZero()) 
	{
		throw(complex);
	}
	else 
	{
		double div_real = (_real * complex._real + _imaginary * complex._imaginary)/(complex.module()* complex.module());
		double div_imaginary = (_real * complex._imaginary - _imaginary * complex._real) /(complex.module() * complex.module());
		return Complex(div_real, div_imaginary);
	}
}

Complex Complex::exponential(const Complex& complex)
{
	double exp_real = std::exp(complex._real)* std::cos(complex._imaginary);
	double exp_imaginary = std::exp(complex._real) * std::sin(complex._imaginary);
	return Complex(exp_real, exp_imaginary);
}

Complex Complex::logarithm(const Complex& complex)
{
	double log_real = log(complex.module());
	double log_imaginary = complex.argument();
	return Complex(log_real, log_imaginary);
}

Complex Complex::square_root(const Complex& complex)
{
	double sqr_real = std::sqrt(complex.module()) * std::cos(complex.argument() / 2);
	double sqr_imaginary = std::sqrt(complex.module()) * std::sin(complex.argument() / 2);
	return Complex(sqr_real, sqr_imaginary);
}

bool Complex::isZero() const
{
	if (_real == 0 && _imaginary == 0) 
	{
		return true;
	}

	return false;
}

double Complex::module() const
{
	return std::sqrt(_real * _real + _imaginary * _imaginary);
}

double Complex::argument() const
{
	return std::atan2(_imaginary,_real);
}

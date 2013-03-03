/*
 * polynomial.h
 * This file is part of libPolynomial
 *
 * Copyright (C) 2012 - prime number
 *
 * libPolynomial is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * libPolynomial is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with libPolynomial. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __POLY_H__
#define __POLY_H__
#include <string>
#include <vector>
#include <stack>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <stdarg.h>

long long _count = 0;

#define THRESHOLD_KARATSUBA		20

template<typename T,typename U=int>
class Polynomial
{
	public:
	Polynomial();
	Polynomial(const Polynomial<T,U> &);
	Polynomial(const std::vector<T> &);
	Polynomial(const T &);
	Polynomial diff()const;
	Polynomial i_integ(const T &)const;
	T d_integ(const T &,const T &)const;
	Polynomial & toMonique()const;
	static Polynomial variable();
	static Polynomial monomial(const int,const T &);
	Polynomial operator+(const Polynomial<T,U> &)const;
	Polynomial operator+(const T &)const;
	Polynomial & operator+=(const Polynomial<T,U> &);
	Polynomial & operator+=(const T &);
	Polynomial operator-(const Polynomial<T,U> &)const;
	Polynomial operator-(const T &)const;
	Polynomial operator-()const;
	Polynomial & operator-=(const Polynomial<T,U> &);
	Polynomial & operator-=(const T &);
	Polynomial operator*(const Polynomial<T,U> &)const;
	Polynomial operator*(const T &)const;
	Polynomial & operator*=(const Polynomial<T,U> &);
	Polynomial operator^(const U &)const;
	Polynomial & operator^=(const U &);
	Polynomial operator/(const Polynomial<T,U> &)const;
	Polynomial operator/(const T &)const;
	Polynomial & operator/=(const Polynomial<T,U> &);
	Polynomial operator%(const Polynomial<T,U> &)const;
	Polynomial & operator%=(const Polynomial<T,U> &);
	Polynomial & operator=(const Polynomial<T,U> &);
	Polynomial & operator=(const T &);
	bool operator==(const Polynomial<T,U> &)const;
	bool operator==(const T &)const;
	bool operator!=(const Polynomial<T,U> &)const;
	bool operator!=(const T &)const;
	T operator()(const T &)const;
	int degree()const;
	std::vector<T> getPoly()const;
	void setPoly(const std::vector<T> &);
	std::string getStr()const;
	bool isFactor(const Polynomial<T,U> &)const;
	private:
	std::vector<T> cefs;
};

template<typename T,typename U>
Polynomial<T,U> operator+(const T & x,const Polynomial<T,U> & y)
{
	return y+x;
}

template<typename T,typename U>
Polynomial<T,U> operator-(const T & x,const Polynomial<T,U> & y)
{
	return -y+x;
}

template<typename T,typename U>
Polynomial<T,U> operator*(const T & x,const Polynomial<T,U> & y)
{
	return y*x;
}

template<typename T,typename U>
Polynomial<T,U> Polynomial<T,U>::variable()
{
	std::vector<T> ary(2);
	ary[0] = 0;
	ary[1] = 1;
	return Polynomial(ary);
}

template<typename T,typename U>
Polynomial<T,U> Polynomial<T,U>::monomial(const int degree,const T & cef)
{
	_count++;
	std::vector<T> ary(degree+1,0);
	ary[degree] = cef;
	return Polynomial(ary);
}

template<typename T,typename U>
Polynomial<T,U>::Polynomial(const std::vector<T> & ary)
{
	_count++;
	cefs = std::vector<T>(ary);
	while(cefs.back() == 0 && cefs.size() > 1)
	{
		cefs.pop_back();
	}
}

template<typename T,typename U>
Polynomial<T,U>::Polynomial(const Polynomial<T,U> & poly)
{
	_count++;
	cefs = std::vector<T>(poly.cefs);
	while(poly.cefs.back() == 0 && poly.cefs.size() > 1)
	{
		cefs.pop_back();
	}
}

template<typename T,typename U>
Polynomial<T,U>::Polynomial(const T & val)
{
	cefs.clear();
	cefs.push_back(val);
}

template<typename T,typename U>
Polynomial<T,U>::Polynomial()
{
	cefs = std::vector<T>(1);
	cefs[0] = 0;
}

template<typename T,typename U>
Polynomial<T,U> Polynomial<T,U>::diff()const
{
	const int size = cefs.size();
	std::vector<T> ary(size-1);
	int i;
	for(i = 1;i < size;i++)
	{
		ary[i-1] = cefs[i]*i;
	}
	return Polynomial(ary);
}

template<typename T,typename U>
Polynomial<T,U> Polynomial<T,U>::i_integ(const T & c)const
{
	const int size = cefs.size();
	std::vector<T> ary(size+1);
	int i;
	for(i = 0;i < size;i++)
	{
		ary[i+1] = cefs[i]/(i+1);
	}
	ary[0] = c;
	return Polynomial(ary);
}

template<typename T,typename U>
T Polynomial<T,U>::d_integ(const T & start,const T & end)const
{
	Polynomial integ = (*this).i_integ(0);
	return integ(end) - integ(start);
}

template<typename T,typename U>
Polynomial<T,U> & Polynomial<T,U>::toMonique()const
{
	return *this / cefs[cefs.size()-1];
}

template<typename T,typename U>
Polynomial<T,U> Polynomial<T,U>::operator+(const Polynomial<T,U> & y)const
{
	Polynomial poly = Polynomial(*this);
	poly += y;
	return poly;
}

template<typename T,typename U>
Polynomial<T,U> Polynomial<T,U>::operator+(const T & y)const
{
	Polynomial poly = Polynomial(*this);
	poly += Polynomial(y);
	return poly;
}

template<typename T,typename U>
Polynomial<T,U> & Polynomial<T,U>::operator+=(const Polynomial<T,U> & y)
{
	const int xsize = cefs.size();
	const int ysize = y.cefs.size();
	if(xsize < ysize)
	{
		cefs.resize(ysize);
	}
	int i;
	for(i = 0;i < ysize;i++)
	{
		cefs[i] += y.cefs[i];
	}
	while(cefs.back() == 0 && cefs.size() > 1)
	{
		cefs.pop_back();
	}
	return *this;
}

template<typename T,typename U>
Polynomial<T,U> & Polynomial<T,U>::operator+=(const T & y)
{
	cefs[0] += y;
	return *this;
}

template<typename T,typename U>
Polynomial<T,U> Polynomial<T,U>::operator-(const Polynomial<T,U> & y)const
{
	Polynomial poly = Polynomial(*this);
	poly -= y;
	return poly;
}

template<typename T,typename U>
Polynomial<T,U> Polynomial<T,U>::operator-(const T & y)const
{
	Polynomial poly = Polynomial(*this);
	poly -= Polynomial(y);
	return poly;
}

template<typename T,typename U>
Polynomial<T,U> Polynomial<T,U>::operator-()const
{
	const int size = cefs.size();
	std::vector<T> ary(size);
	int i;
	for(i = 0;i < size;i++)
	{
		ary[i] = -cefs[i];
	}
	Polynomial poly = Polynomial(ary);
	return poly;
}

template<typename T,typename U>
Polynomial<T,U> & Polynomial<T,U>::operator-=(const T & y)
{
	cefs[0] -= y;
	return *this;
}

template<typename T,typename U>
Polynomial<T,U> & Polynomial<T,U>::operator-=(const Polynomial<T,U> & y)
{
	const int xsize = cefs.size();
	const int ysize = y.cefs.size();
	if(xsize < ysize)
	{
		cefs.resize(ysize);
	}
	int i;
	for(i = 0;i < ysize;i++)
	{
		cefs[i] -= y.cefs[i];
	}
	while(cefs.back() == 0 && cefs.size() > 1)
	{
		cefs.pop_back();
	}
	return *this;
}

template<typename T,typename U>
Polynomial<T,U> Polynomial<T,U>::operator*(const Polynomial<T,U> & y)const
{
	Polynomial poly = Polynomial(*this);
	poly *= y;
	return poly;
}

template<typename T,typename U>
Polynomial<T,U> Polynomial<T,U>::operator*(const T & val)const
{
	Polynomial poly = Polynomial(*this);
	poly *= Polynomial(val);
	return poly;
}

template<typename T,typename U>
Polynomial<T,U> & Polynomial<T,U>::operator*=(const Polynomial<T,U> & y)
{
	const int xsize = cefs.size();
	const int ysize = y.cefs.size();
	std::vector<T> ary(xsize+ysize-1,0);
	if(xsize > THRESHOLD_KARATSUBA && ysize > THRESHOLD_KARATSUBA)
	{
		if(xsize*2 >= ysize && ysize*2 > xsize)
		{
			int half = xsize/2;
			typename std::vector<T>::iterator itx = cefs.begin();
			Polynomial x1p = Polynomial(std::vector<T>(itx,itx+half));
			Polynomial x2p = Polynomial(std::vector<T>(itx+half,cefs.end()));
			Polynomial y1p = Polynomial(std::vector<T>(y.cefs.begin(),y.cefs.begin()+half));
			Polynomial y2p = Polynomial(std::vector<T>(y.cefs.begin()+half,y.cefs.end()));
			std::vector<T> a = (x1p*y1p).getPoly();
			std::vector<T> b = ((x2p-x1p)*(y2p-y1p)).getPoly();
			std::vector<T> c = (x2p*=y2p).getPoly();
			int i;
			const int la = a.size();
			for(i = 0;i < la;i++)
			{
				ary[i] += a[i];
				ary[i+half] += a[i];
			}
			const int lb = b.size();
			for(i = 0;i < lb;i++)
			{
				ary[i+half] -= b[i];
			}
			const int lc = c.size();
			for(i = 0;i < lc;i++)
			{
				ary[i+half] += c[i];
				ary[i+half*2] += c[i];
			}
		}
		else if(ysize*2 <= xsize)
		{
			int half = xsize/2;
			typename std::vector<T>::iterator itx = cefs.begin();
			std::vector<T> x1 = std::vector<T>(itx,itx+half);
			std::vector<T> x2 = std::vector<T>(itx+half,cefs.end());
			Polynomial x1p = Polynomial(x1);
			Polynomial x2p = Polynomial(x2);
			std::vector<T> a = (x1p*y).getPoly();
			std::vector<T> b = (x2p*y).getPoly();
			int i;
			const int la = a.size();
			for(i = 0;i < la;i++)
			{
				ary[i] += a[i];
			}
			const int lb = b.size();
			for(i = 0;i < lb;i++)
			{
				ary[i+half] += b[i];
			}
		}
		else
		{
			int half = ysize/2;
			std::vector<T> y1 = std::vector<T>(y.cefs.begin(),y.cefs.begin()+half);
			std::vector<T> y2 = std::vector<T>(y.cefs.begin()+half,y.cefs.end());
			Polynomial y1p = Polynomial(y1);
			Polynomial y2p = Polynomial(y2);
			std::vector<T> a = ((*this)*y1p).getPoly();
			a.resize(half*2,0);
			std::vector<T> b = ((*this)*y2p).getPoly();
			b.resize(ysize,0);
			int i;
			const int la = a.size();
			for(i = 0;i < la;i++)
			{
				ary[i] += a[i];
			}
			const int lb = b.size();
			for(i = 0;i < lb;i++)
			{
				ary[i+half] += b[i];
			}
		}
	}
	else
	{
		int i,j;
		for(i = 0;i < xsize;i++)
		{
			for(j = 0;j < ysize;j++)
			{
				ary[i+j] += cefs[i]*y.cefs[j];
			}
		}
	}
	cefs = std::vector<T>(ary);
	return *this;
}

template<typename T,typename U>
Polynomial<T,U> Polynomial<T,U>::operator^(const U & index)const
{
	Polynomial poly = Polynomial(*this);
	poly ^= index;
	return poly;
}

template<typename T,typename U>
Polynomial<T,U> & Polynomial<T,U>::operator^=(const U & index)
{
	if(index <= 0)
	{
		(*this) = Polynomial(1);
		return *this;
	}
	else if(index == 1)
	{
		return *this;
	}
	else
	{
		Polynomial poly = Polynomial(*this);
		std::vector<bool> ary;
		int k = 0;
		U i = index;
		while(i > 0)
		{
			if((i % 2) != 0)
			{
				ary.push_back(true);
			}
			else
			{
				ary.push_back(false);
			}
			k++;
			i /= 2;
		}
		int j;
		for(j = k-2;j >= 0;j--)
		{
			*this *= *this;
			if(ary[j])
			{
				*this *= poly;
			}
		}
		return *this;
	}
}

template<typename T,typename U>
Polynomial<T,U> Polynomial<T,U>::operator/(const Polynomial<T,U> & y)const
{
	Polynomial poly = Polynomial(*this);
	poly /= y;
	return poly;
}

template<typename T,typename U>
Polynomial<T,U> Polynomial<T,U>::operator/(const T & val)const
{
	Polynomial poly = Polynomial(*this);
	poly /= Polynomial(val);
	return poly;
}

template<typename T,typename U>
Polynomial<T,U> & Polynomial<T,U>::operator/=(const Polynomial<T,U> & y)
{
	const int xsize = cefs.size();
	const int ysize = y.cefs.size();
	if(xsize < ysize)
	{
		*this = Polynomial(T(0));
		return *this;
	}
	const int deg = ysize-1;
	const int loop = xsize-deg;
	std::vector<T> ary(loop);
	Polynomial poly_ = Polynomial(*this);
	int i;
	for(i = loop-1;i >= 0;i--)
	{
		ary[i] = (poly_.cefs[i+deg])/(y.cefs[deg]);
		poly_ -= Polynomial::monomial(i,1)*y*ary[i];
	}
	cefs = std::vector<T>(ary);
	return *this;
}

template<typename T,typename U>
Polynomial<T,U> Polynomial<T,U>::operator%(const Polynomial<T,U> & y)const
{
	Polynomial poly = Polynomial(*this);
	poly %= y;
	return poly;
}

template<typename T,typename U>
Polynomial<T,U> & Polynomial<T,U>::operator%=(const Polynomial<T,U> & y)
{
	const int xsize = cefs.size();
	const int ysize = y.cefs.size();
	if(xsize < ysize)
	{
		return *this;
	}
	const int deg = ysize-1;
	const int loop = xsize-deg;
	std::vector<T> ary(loop);
	Polynomial poly_ = Polynomial(*this);
	int i;
	for(i = loop-1;i >= 0;i--)
	{
		ary[i] = (poly_.cefs[i+deg])/(y.cefs[deg]);
		poly_ -= Polynomial::monomial(i,1)*y*ary[i];
	}
	*this = poly_;
	return *this;
}

template<typename T,typename U>
Polynomial<T,U> & Polynomial<T,U>::operator=(const Polynomial<T,U> & y)
{
	cefs.clear();
	int i;
	const int ysize = y.cefs.size();
	for(i = 0;i < ysize;i++)
	{
		cefs.push_back(y.cefs[i]);
	}
	while(cefs.back() == 0 && cefs.size() > 1)
	{
		cefs.pop_back();
	}
	return *this;
}

template<typename T,typename U>
Polynomial<T,U> & Polynomial<T,U>::operator=(const T & val)
{
	cefs.clear();
	cefs.push_back(val);
	return *this;
}

template<typename T,typename U>
bool Polynomial<T,U>::operator==(const Polynomial<T,U> & y)const
{
	return cefs == y.cefs;
}

template<typename T,typename U>
bool Polynomial<T,U>::operator!=(const Polynomial<T,U> & y)const
{
	return cefs != y.cefs;
}

template<typename T,typename U>
bool Polynomial<T,U>::operator==(const T & y)const
{
	return cefs == std::vector<T>(1,y);
}

template<typename T,typename U>
bool Polynomial<T,U>::operator!=(const T & y)const
{
	return cefs != std::vector<T>(1,y);
}

template<typename T,typename U>
T Polynomial<T,U>::operator()(const T & x)const
{
	T val = 0;
	T expx = 1;
	int i;
	const int size = cefs.size();
	for(i = 0;i < size;i++)
	{
		val += cefs[i]*expx;
		expx *= x;
	}
	return val;
}

template<typename T,typename U>
int Polynomial<T,U>::degree()const
{
	return cefs.size()-1;
}

template<typename T,typename U>
std::vector<T> Polynomial<T,U>::getPoly()const
{
	return cefs;
}

template<typename T,typename U>
void Polynomial<T,U>::setPoly(const std::vector<T> & ary)
{
	cefs = ary;
	while(cefs.back() == 0 && cefs.size() > 1)
	{
		cefs.pop_back();
	}
}

template<typename T,typename U>
bool Polynomial<T,U>::isFactor(const Polynomial<T,U> & y)const
{
	return *this == ((*this)/y)*y && (*this)%y == 0;
}

template<typename T,typename U>
std::string Polynomial<T,U>::getStr()const
{
	bool flg = false;
	std::stringstream ss;
	int i = cefs.size()-1;
	if(i > 1)
	{
		if(cefs[i] != 0)
		{
			flg = true;
			if(cefs[i] == 1)
				ss << "x^" << i;
			else if(cefs[i] == -1)
				ss << "-x^" << i;
			else
				ss << cefs[i] << "x^" << i;
		}
	}
	for(i--;i > 1;i--)
	{
		if(cefs[i] != 0)
		{
			flg = true;
			if(cefs[i] == 1)
				ss << "+x^" << i;
			else if(cefs[i] == -1)
				ss << "-x^" << i;
			else
				ss << std::showpos << cefs[i] << "x^" << std::noshowpos << i;
		}
	}
	if(cefs.size() > 2)
	{
		if(cefs[1] != 0)
		{
			flg = true;
			if(cefs[1] == 1)
				ss << "+x";
			else if(cefs[1] == -1)
				ss << "-x";
			else
				ss << std::showpos << cefs[1] << "x";
		}
	}
	else if(cefs.size() == 2)
	{
		if(cefs[1] != 0)
		{
			flg = true;
			if(cefs[1] == 1)
				ss << "x";
			else if(cefs[1] == -1)
				ss << "-x";
			else
				ss << std::noshowpos << cefs[1] << "x";
		}
	}
	if(flg == false)
		ss << cefs[0];
	else if(cefs[0] != 0)
		ss << std::showpos << cefs[0];
	return ss.str();
}

#endif

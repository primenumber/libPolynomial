#ifndef __POLY_H__
#define __POLY_H__
#include <string>
#include <vector>
#include <stack>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <stdarg.h>

using namespace std;

template<class T,class U=int>
class Polynomial
{
	public:
	Polynomial();
	Polynomial(const Polynomial<T,U> &);
	Polynomial(const vector<T> &);
	Polynomial(const T);
	Polynomial diff()const;
	Polynomial i_integ(T)const;
	T d_integ(T,T)const;
	Polynomial toMonique()const;
	static Polynomial variable();
	Polynomial operator+(const Polynomial<T,U>)const;
	Polynomial operator+(const T)const;
	Polynomial & operator+=(const Polynomial<T,U>);
	Polynomial operator-(const Polynomial<T,U>)const;
	Polynomial operator-(const T)const;
	Polynomial operator-()const;
	Polynomial & operator-=(const Polynomial<T,U>);
	Polynomial operator*(const Polynomial<T,U>)const;
	Polynomial operator*(const T)const;
	Polynomial & operator*=(const Polynomial<T,U>);
	Polynomial operator^(const U)const;
	Polynomial operator/(const Polynomial<T,U>)const;
	Polynomial operator/(const T)const;
	Polynomial & operator/=(const Polynomial<T,U>);
	Polynomial operator%(const Polynomial<T,U>)const;
	Polynomial & operator%=(const Polynomial<T,U>);
	Polynomial & operator=(const Polynomial<T,U>);
	Polynomial & operator=(const T);
	bool operator==(const Polynomial<T,U>)const;
	bool operator!=(const Polynomial<T,U>)const;
	bool operator==(const T)const;
	bool operator!=(const T)const;
	T operator()(T)const;
	unsigned int degree()const;
	vector<T> & getPoly()const;
	void setPoly(const vector<T> &);
	string getStr()const;
	Polynomial map_cefs(T (*)(T,int,T),T)const;
	Polynomial map_cefs(T (*)(T,int))const;
	Polynomial map_cefs(T (*)(T))const;
	bool isFactor(const Polynomial<T,U>)const;
	private:
	vector<T> cefs;
};

template<class T,class U>
Polynomial<T,U> operator+(const T x,const Polynomial<T,U> y)
{
	return y+x;
}

template<class T,class U>
Polynomial<T,U> operator-(const T x,const Polynomial<T,U> y)
{
	return -y+x;
}

template<class T,class U>
Polynomial<T,U> operator*(const T x,const Polynomial<T,U> y)
{
	return y*x;
}

template<class T,class U>
Polynomial<T,U> Polynomial<T,U>::variable()
{
	vector<T> ary(2);
	ary[0] = 0;
	ary[1] = 1;
	return Polynomial(ary);
}

template<class T,class U>
Polynomial<T,U>::Polynomial(const vector<T> & ary)
{
	cefs.clear();
	int i;
	const int size = ary.size();
	for(i = 0;i < size;i++)
	{
		cefs.push_back(ary[i]);
	}
	while(cefs.back() == 0 && cefs.size() > 1)
	{
		cefs.pop_back();
	}
}

template<class T,class U>
Polynomial<T,U>::Polynomial(const Polynomial<T,U> &poly)
{
	cefs.clear();
	int i;
	const int size = poly.cefs.size();
	for(i = 0;i < size;i++)
	{
		cefs.push_back(poly.cefs[i]);
	}
	while(cefs.back() == 0 && cefs.size() > 1)
	{
		cefs.pop_back();
	}
}

template<class T,class U>
Polynomial<T,U>::Polynomial(const T val)
{
	cefs.clear();
	cefs.push_back(val);
}

template<class T,class U>
Polynomial<T,U>::Polynomial()
{
	cefs = vector<T>(1);
	cefs[0] = 0;
}

template<class T,class U>
Polynomial<T,U> Polynomial<T,U>::diff()const
{
	const int size = cefs.size();
	vector<T> ary(size-1);
	int i;
	for(i = 1;i < size;i++)
	{
		ary[i-1] = cefs[i]*i;
	}
	return Polynomial(ary);
}

template<class T,class U>
Polynomial<T,U> Polynomial<T,U>::i_integ(T c)const
{
	const int size = cefs.size();
	vector<T> ary(size+1);
	int i;
	for(i = 0;i < size;i++)
	{
		ary[i+1] = cefs[i]/(i+1);
	}
	ary[0] = c;
	return Polynomial(ary);
}

template<class T,class U>
T Polynomial<T,U>::d_integ(T start,T end)const
{
	Polynomial integ = (*this).i_integ(0);
	return integ(end) - integ(start);
}

template<class T,class U>
Polynomial<T,U> Polynomial<T,U>::toMonique()const
{
	return *this / cefs[cefs.size()-1];
}

template<class T,class U>
Polynomial<T,U> Polynomial<T,U>::operator+(const Polynomial<T,U> y)const
{
	const int xsize = cefs.size();
	const int ysize = y.cefs.size();
	vector<T> ary(max(xsize,ysize));
	int i;
	for(i = 0;i < xsize;i++)
	{
		ary[i] = cefs[i];
	}
	for(i = 0;i < ysize;i++)
	{
		ary[i] += y.cefs[i];
	}
	while(ary.back() == 0 && ary.size() > 1)
	{
		ary.pop_back();
	}
	return Polynomial(ary);
}

template<class T,class U>
Polynomial<T,U> Polynomial<T,U>::operator+(const T y)const
{
	const int size = cefs.size();
	vector<T> ary(size);
	int i;
	for(i = 0;i < size;i++)
	{
		ary[i] = cefs[i];
	}
	ary[0] += y;
	return Polynomial(ary);
}

template<class T,class U>
Polynomial<T,U> & Polynomial<T,U>::operator+=(const Polynomial<T,U> y)
{
	*this = *this + y;
	return *this;
}

template<class T,class U>
Polynomial<T,U> Polynomial<T,U>::operator-(const Polynomial<T,U> y)const
{
	const int xsize = cefs.size();
	const int ysize = y.cefs.size();
	vector<T> ary(max(xsize,ysize));
	int i;
	for(i = 0;i < xsize;i++)
	{
		ary[i] = cefs[i];
	}
	for(i = 0;i < ysize;i++)
	{
		ary[i] -= y.cefs[i];
	}
	while(ary.back() == 0 && ary.size() > 1)
	{
		ary.pop_back();
	}
	return Polynomial(ary);
}

template<class T,class U>
Polynomial<T,U> Polynomial<T,U>::operator-(const T y)const
{
	const int size = cefs.size();
	vector<T> ary(size);
	int i;
	for(i = 0;i < size;i++)
	{
		ary[i] = cefs[i];
	}
	ary[0] -= y;
	return Polynomial(ary);
}

template<class T,class U>
Polynomial<T,U> Polynomial<T,U>::operator-()const
{
	const int size = cefs.size();
	vector<T> ary(size);
	int i;
	for(i = 0;i < size;i++)
	{
		ary[i] = -cefs[i];
	}
	return Polynomial(ary);
}

template<class T,class U>
Polynomial<T,U> & Polynomial<T,U>::operator-=(const Polynomial<T,U> y)
{
	*this = *this - y;
	return *this;
}

template<class T,class U>
Polynomial<T,U> Polynomial<T,U>::operator*(const Polynomial<T,U> y)const
{
	const int xsize = cefs.size();
	const int ysize = y.cefs.size();
	vector<T> ary(xsize+ysize-1);
	int i,j;
	for(i = 0;i < xsize;i++)
	{
		for(j = 0;j < ysize;j++)
		{
			ary[i+j] += cefs[i]*y.cefs[j];
		}
	}
	return Polynomial(ary);
}

template<class T,class U>
Polynomial<T,U> Polynomial<T,U>::operator*(const T val)const
{
	const int size = cefs.size();
	vector<T> ary(size);
	int i;
	for(i = 0;i < size;i++)
	{
		ary[i] = cefs[i]*val;
	}
	return Polynomial(ary);
}

template<class T,class U>
Polynomial<T,U> & Polynomial<T,U>::operator*=(const Polynomial<T,U> y)
{
	*this = *this * y;
	return *this;
}

template<class T,class U>
Polynomial<T,U> Polynomial<T,U>::operator^(const U index)const
{
	if(index <= 0)
		return Polynomial();
	else if(index == 1)
		return Polynomial(*this);
	else
	{
		Polynomial half = Polynomial(*this)^(index/2);
		half *= half;
		if((index % 2) == 1)
		{
			half *= Polynomial(*this);
		}
		return half;
	}
}

template<class T,class U>
Polynomial<T,U> Polynomial<T,U>::operator/(const Polynomial<T,U> y)const
{
	const int xsize = cefs.size();
	const int ysize = y.cefs.size();
	if(xsize < ysize)
		return Polynomial(0);
	const int deg = ysize-1;
	const int loop = xsize-deg;
	const Polynomial<T,U> x = Polynomial<T,U>::variable();
	vector<T> ary(loop);
	Polynomial poly = *this;
	int i;
	for(i = loop-1;i >= 0;i--)
	{
		ary[i] = poly.cefs[i+deg]/y.cefs[deg];
		poly -= (x^i)*y*ary[i];
	}
	return Polynomial(ary);
}

template<class T,class U>
Polynomial<T,U> Polynomial<T,U>::operator/(const T y)const
{
	const int size = cefs.size();
	vector<T> ary(size);
	int i;
	for(i = 0;i < size;i++)
	{
		ary[i] = cefs[i]/y;
	}
	return Polynomial(ary);
}


template<class T,class U>
Polynomial<T,U> & Polynomial<T,U>::operator/=(const Polynomial<T,U> y)
{
	*this = *this / y;
	return *this;
}
template<class T,class U>
Polynomial<T,U> Polynomial<T,U>::operator%(const Polynomial<T,U> y)const
{
	return Polynomial(*this - (*this / y) * y);
}

template<class T,class U>
Polynomial<T,U> & Polynomial<T,U>::operator%=(const Polynomial<T,U> y)
{
	*this = *this % y;
	return *this;
}

template<class T,class U>
Polynomial<T,U> & Polynomial<T,U>::operator=(const Polynomial<T,U> y)
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

template<class T,class U>
Polynomial<T,U> & Polynomial<T,U>::operator=(const T val)
{
	cefs.clear();
	cefs.push_back(val);
	return *this;
}

template<class T,class U>
bool Polynomial<T,U>::operator==(const Polynomial<T,U> y)const
{
	return cefs == y.cefs;
}

template<class T,class U>
bool Polynomial<T,U>::operator!=(const Polynomial<T,U> y)const
{
	return cefs != y.cefs;
}

template<class T,class U>
bool Polynomial<T,U>::operator==(const T y)const
{
	return cefs == vector<T>(1,y);
}

template<class T,class U>
bool Polynomial<T,U>::operator!=(const T y)const
{
	return cefs != vector<T>(1,y);
}

template<class T,class U>
T Polynomial<T,U>::operator()(const T x)const
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

template<class T,class U>
unsigned int Polynomial<T,U>::degree()const
{
	return cefs.size()-1;
}

template<class T,class U>
vector<T> & Polynomial<T,U>::getPoly()const
{
	return cefs;
}

template<class T,class U>
void Polynomial<T,U>::setPoly(const vector<T> & ary)
{
	cefs = ary;
	while(cefs.back() == 0 && cefs.size() > 1)
	{
		cefs.pop_back();
	}
}

template<class T,class U>
Polynomial<T,U> Polynomial<T,U>::map_cefs(T (*func)(T,int,T),T data)const
{
	const int size = cefs.size();
	vector<T> ary(size);
	int i;
	for(i = 0;i < size;i++)
	{
		ary[i] = (*func)(cefs[i],i,data);
	}
	return Polynomial(ary);
}

template<class T,class U>
Polynomial<T,U> Polynomial<T,U>::map_cefs(T (*func)(T,int))const
{
	const int size = cefs.size();
	vector<T> ary(size);
	int i;
	for(i = 0;i < size;i++)
	{
		ary[i] = (*func)(cefs[i],i);
	}
	return Polynomial(ary);
}

template<class T,class U>
Polynomial<T,U> Polynomial<T,U>::map_cefs(T (*func)(T))const
{
	const int size = cefs.size();
	vector<T> ary(size);
	int i;
	for(i = 0;i < size;i++)
	{
		ary[i] = (*func)(cefs[i]);
	}
	return Polynomial(ary);
}

template<class T,class U>
bool Polynomial<T,U>::isFactor(const Polynomial<T,U> y)const
{
	return *this == ((*this)/y)*y && (*this)%y == 0;
}

template<class T,class U>
string Polynomial<T,U>::getStr()const
{
	bool flg = false;
	stringstream ss;
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
				ss << showpos << cefs[i] << "x^" << noshowpos << i;
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
				ss << showpos << cefs[1] << "x";
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
				ss << noshowpos << cefs[1] << "x";
		}
	}
	if(flg == false)
		ss << cefs[0];
	else if(cefs[0] != 0)
		ss << showpos << cefs[0];
	return ss.str();
}

#endif

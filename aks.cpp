#include <iostream>
#include <stdio.h>
#include <gmp.h>
#include <gmpxx.h>
#include "polynomial.h"

using namespace std;

mpz_class mod_mpz(const mpz_class & x,const mpz_class & data)
{
	mpz_class r;
	mpz_mod(r.get_mpz_t(),x.get_mpz_t(),data.get_mpz_t());
	return r;
}

class ModuloInt
{
	public:
	friend std::ostream& operator<<(std::ostream &, const ModuloInt &);
	static mpz_class defaultMod;
	ModuloInt();
	ModuloInt(const ModuloInt &);
	ModuloInt(const mpz_class &,const mpz_class &);
	ModuloInt(const int,const mpz_class &);
	ModuloInt operator+(const ModuloInt &)const;
	ModuloInt & operator+=(const ModuloInt &);
	ModuloInt operator-(const ModuloInt &)const;
	ModuloInt & operator-=(const ModuloInt &);
	ModuloInt operator*(const ModuloInt &)const;
	ModuloInt & operator*=(const ModuloInt &);
	ModuloInt operator/(const ModuloInt &)const;
	ModuloInt & operator/=(const ModuloInt &);
	ModuloInt & operator=(const ModuloInt &);
	ModuloInt & operator=(const mpz_class &);
	ModuloInt & operator=(const int);
	bool operator==(const ModuloInt &)const;
	bool operator==(const int)const;
	bool operator!=(const ModuloInt &)const;
	bool operator!=(const int)const;
	private:
	mpz_class val;
	mpz_class _mod;
};

std::ostream& operator<<(std::ostream & os, const ModuloInt & N)
{
	os << N.val.get_str();
	return os;
}

ModuloInt::ModuloInt()
{
	val = 0;
	_mod = defaultMod;
}

ModuloInt::ModuloInt(const ModuloInt & N)
{
	val = N.val;
	_mod = N._mod;
}

ModuloInt::ModuloInt(const mpz_class & N,const mpz_class & mod = defaultMod)
{
	val = mod_mpz(N,mod);
	_mod = mod;
}

ModuloInt::ModuloInt(const int N,const mpz_class & mod = defaultMod)
{
	val = mod_mpz(N,mod);
	_mod = mod;
}

ModuloInt ModuloInt::operator+(const ModuloInt & y)const
{
	return ModuloInt(mod_mpz((*this).val + y.val,_mod),_mod);
}

ModuloInt & ModuloInt::operator+=(const ModuloInt & y)
{
	val = mod_mpz(val + y.val,_mod);
	return *this;
}

ModuloInt ModuloInt::operator-(const ModuloInt & y)const
{
	return ModuloInt(mod_mpz((*this).val - y.val,_mod),_mod);
}

ModuloInt & ModuloInt::operator-=(const ModuloInt & y)
{
	val = mod_mpz(val - y.val,_mod);
	return *this;
}

ModuloInt ModuloInt::operator*(const ModuloInt & y)const
{
	return ModuloInt(mod_mpz((*this).val * y.val,_mod),_mod);
}

ModuloInt & ModuloInt::operator*=(const ModuloInt & y)
{
	val = mod_mpz(val * y.val,_mod);
	return *this;
}

ModuloInt ModuloInt::operator/(const ModuloInt & y)const
{
	return ModuloInt(mod_mpz((*this).val / y.val,_mod),_mod);
}

ModuloInt & ModuloInt::operator/=(const ModuloInt & y)
{
	val = mod_mpz(val / y.val,_mod);
	return *this;
}

ModuloInt & ModuloInt::operator=(const ModuloInt & y)
{
	val = y.val;
	_mod = y._mod;
	return *this;
}

ModuloInt & ModuloInt::operator=(const mpz_class & y)
{
	val = y;
	return *this;
}

ModuloInt & ModuloInt::operator=(const int y)
{
	val = y;
	return *this;
}

bool ModuloInt::operator==(const ModuloInt & y)const
{
	return val == y.val;
}

bool ModuloInt::operator==(const int y)const
{
	return val == y;
}

bool ModuloInt::operator!=(const ModuloInt & y)const
{
	return val != y.val;
}

bool ModuloInt::operator!=(const int y)const
{
	return val != y;
}

mpz_class ModuloInt::defaultMod = 2;

mpf_class log(const mpf_class x,int depth = 0)
{
	if(x > 1.25)
	{
		return 2.0*log(sqrt(x));
	}
	else
	{
		mpf_class val = 0.0;
		mpf_class xexp = x-1.0;
		int i;
		for(i = 1;i < depth+10;i++)
		{
			val += xexp/i;
			xexp *= 1.0-x;
		}
		return val;
	}
}

mpz_class findr(mpz_class n)
{
	mpz_class e_min = mpz_class((log(mpf_class(n))*log(mpf_class(n)))*4.0+1.0);
	mpz_class i=e_min;
	while(1)
	{
		mpz_class gcd;
		mpz_gcd(gcd.get_mpz_t(),i.get_mpz_t(),n.get_mpz_t());
		if(1 < gcd && gcd < n)
			return 0;
		mpz_class j = 1;
		mpz_class nexp = n%i;
		if(gcd != n)
		{
			while(1)
			{
				if(nexp == 1)
				{
					if(j >= e_min)
						return i;
					else
						break;
				}
				else
				{
					nexp *= n;
					nexp %= i;
					j += 1;
				}
			}
		}
		i += 1;
	}
}

Polynomial<ModuloInt,mpz_class> polymodexpr(
	const Polynomial<ModuloInt,mpz_class> & poly,
	const mpz_class & r)
{
	vector<ModuloInt> ary = poly.getPoly();
	if(ary.size() <= r.get_si())
		return Polynomial<ModuloInt,mpz_class>(poly);
	for(int i = ary.size()-1;i >= r;i--)
	{
		ary[i-r.get_si()] += ary[i];
		ary.pop_back();
	}
	return Polynomial<ModuloInt,mpz_class>(ary);
}

Polynomial<ModuloInt,mpz_class> powmodpoly(
	const Polynomial<ModuloInt,mpz_class> & poly,
	const mpz_class & index,
	const mpz_class & r)
{
	if(index == 1)
	{
		return poly;
	}
	else
	{
		Polynomial<ModuloInt,mpz_class> half = powmodpoly(poly,index/2,r);
		Polynomial<ModuloInt,mpz_class> val = polymodexpr(half^2,r);
		if((index % 2) == 1)
		{
			val *= poly;
			val = polymodexpr(val,r);
		}
		return val;
	}
}

mpz_class phi(mpz_class n)
{
	mpz_class i;
	mpz_class val = n;
	for(i = 2;i < sqrt(n);i += 1)
	{
		if((n % i) == 0)
		{
			val = (val * (i-1)) / i;
			while((n % i) == 0)
			{
				n /= i;
			}
		}
	}
	if(n == val)
		val = n-1;
	return val;
}

int main()
{
	mpz_class n;
	string s;
	do
	{
		cout << "primary testing:input number(n >= 2)" << endl;
		cin >> s;
		n.set_str(s,10);
	}
	while(n < 2);
	cout << "start" << endl;
	ModuloInt::defaultMod = n;
	if(mpz_perfect_power_p(n.get_mpz_t()))//1st step
	{
		cout << s << " isn't prime." << endl;
		return 0;
	}
	mpz_class r = findr(n);//2nd/3rd step
	if(r == 0)
	{
		cout << s << " isn't prime." << endl;
		return 0;
	}
	mpz_class i;
	for(i = 2;i < r;i++)
	{
		mpz_class gcd;
		mpz_gcd(gcd.get_mpz_t(),i.get_mpz_t(),n.get_mpz_t());
		if(1 < gcd && gcd < n)
		{
			cout << s << " isn't prime." << endl;
			return 0;
		}
	}
	if(n < r)//4th step
	{
		cout << s << " is prime!" << endl;
		return 0;
	}
	mpz_class amax = mpz_class(2*sqrt(phi(r))*log(n));
	cout << "final step:take some time,r = " << r.get_str() << ",amax = " << amax.get_str() << endl;
	Polynomial<ModuloInt,mpz_class> x = Polynomial<ModuloInt,mpz_class>::variable();
	ModuloInt one = ModuloInt(1);
	Polynomial<ModuloInt,mpz_class> xexpr = (Polynomial<ModuloInt,mpz_class>::monomial(r.get_si(),1))-one;
	Polynomial<ModuloInt,mpz_class> xexpnmod = x^(n%r);
	for(i = 1;i < amax;i += 1)
	{
		ModuloInt m_I = ModuloInt(i);
		if(powmodpoly(x+m_I,n,r) != (xexpnmod+m_I))
		{
			cout << s << " isn't prime." << endl;
			return 0;
		}
	}
	cout << s << " is prime!" << _count << endl;
	return 0;
}

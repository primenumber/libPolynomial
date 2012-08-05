#include <iostream>
#include <stdio.h>
#include <gmp.h>
#include <gmpxx.h>
#include "polynomial.h"

mpz_class mod_mpz(mpz_class x,int index,mpz_class data)
{
	mpz_class r;
	mpz_mod(r.get_mpz_t(),x.get_mpz_t(),data.get_mpz_t());
	return r;
}

mpf_class log(const mpf_class x)
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
		for(i = 1;i < 10;i++)
		{
			val += xexp/i;
			xexp *= 1.0-x;
		}
		return val;
	}
}

mpz_class findr(mpz_class n)
{
	mpz_class e_min = 4*(log(mpf_class(n.get_str()))*log(mpf_class(n.get_str())))+1;
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

Polynomial<mpz_class,mpz_class> powmodpoly(
	Polynomial<mpz_class,mpz_class> poly,
	mpz_class index,
	Polynomial<mpz_class,mpz_class> modpoly,
	mpz_class mod)
{
	if(index == 1)
	{
		return poly.map_cefs(mod_mpz,mod);
	}
	else
	{
		Polynomial<mpz_class,mpz_class> half = powmodpoly(poly,index/2,modpoly,mod);
		Polynomial<mpz_class,mpz_class> val = (half^2)%modpoly;
		val = val.map_cefs(mod_mpz,mod);
		if((index % 2) == 1)
		{
			val *= poly;
			val %= modpoly;
			val = val.map_cefs(mod_mpz,mod);
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
	cin >> s;
	n.set_str(s,10);
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
	mpz_class amax = 2*sqrt(phi(r))*log(n);
	Polynomial<mpz_class,mpz_class> x = Polynomial<mpz_class,mpz_class>::variable();
	mpz_class one = 1;
	Polynomial<mpz_class,mpz_class> xexpr = (x^r)-one;
	Polynomial<mpz_class,mpz_class> xexpnmod = powmodpoly(x,n,xexpr,n);
	for(i = 1;i < amax;i += 1)
	{
		if(powmodpoly(x+i,n,xexpr,n) != (xexpnmod+i)%xexpr)
		{
			cout << s << " isn't prime." << endl;
			return 0;
		}
	}
	cout << s << " is prime!" << endl;
	return 0;
}

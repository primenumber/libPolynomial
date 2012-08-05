#include <iostream>
#include "polynomial.h"

int main()
{
	const Polynomial<int> x = Polynomial<int>::variable();
	const Polynomial<double> w = Polynomial<double>::variable();
	Polynomial<int> f = 3*(x^3)+5*(x^2)+8*x+4;
	Polynomial<double> fd = 3.0*(w^3)+5.0*(w^2)+8.0*w+4.0;
	Polynomial<int> g = 3*(x^3)+7*(x^2)+4*x+6;
	Polynomial<int> h = (x^2)+x+2;								//Output:
	cout << f.getStr() << endl;									//3x^3+5x^2+8x+4
	cout << g.getStr() << endl;									//3x^3+7x^2+4x+6
	cout << h.getStr() << endl;									//x^2+x+2
	cout << (f+g).getStr() << endl;								//6x^3+12x^2+12x+10
	cout << (f-g).getStr() << endl;								//-2x^2+4x-2
	cout << (f*g).getStr() << endl;								//9x^6+36x^5+71x^4+106x^3+90x^2+64x+24
	cout << (h^2).getStr() << endl;								//x^4+2x^3+5x^2+4x+4
	cout << (f/h).getStr() << endl;								//3x+2
	cout << (f%g).getStr() << endl;								//-2x^2+4x-2
	cout << f.diff().getStr() << endl;							//9x^2+10x+8
	cout << fd.i_integ(0).getStr() << endl;					//0.75x^4+1.66667x^3+4x^2+4x
	cout << fd.d_integ(-2,2) << endl;							//42.6667
	cout << f(5) << endl;											//544
	cout << (f == g) << endl;										//0
	cout << (f == (3*(x^3)+5*(x^2)+8*x+4)) << endl;			//1
	return 0;
}

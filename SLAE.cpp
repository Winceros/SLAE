// SLAE.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <iostream>
#include <vector>
using namespace std;

#define h 1.0/n
#define PRINT_X_PROGONKA
#define PRINT_MATRIX
const int n = 40;

vector<vector<double>> getAG() {
	vector<double> a;
	vector<double> g;
	a.push_back(0);
	g.push_back(0);

	for(int i=1 ; i<=n ; i++)
		a.push_back(1+i*h);

	g = a;

	vector<vector<double>> res;
	res.push_back(a);
	res.push_back(g);

	return res;
}

vector<double> get_f() {
	vector<double> f;
	f.push_back(0);
	for(int i=1 ; i<=n ; i++)
		f.push_back( h*h*(10*(i*h)*(i*h) + 7*(i*h) - (i*h)*(i*h)*(i*h)*(i*h) - 1) );
	return f;
}

vector<vector<double>> getABC() {
	vector<double> A;
	vector<double> B;
	vector<double> C;
	
	vector<vector<double>> v = getAG();
	vector<double> a = v[0];
	vector<double> g = v[1];

	A.push_back(0);
	B.push_back(0);
	C.push_back(0);

	A.push_back(0);
	B.push_back( a[1]+a[2]+h*h*g[1] );
	C.push_back( -a[2] );
	for(int i = 2 ; i<=n-1 ; i++){
		A.push_back( -a[i] );
		B.push_back( a[i] + a[i+1] + h*h*g[i] );
		C.push_back( -a[i+1] );
	}
	A.push_back( a[n-1] + a[n] + h*h*g[n-1] );
	B.push_back( a[n-1]);
	C.push_back(0);
	
	vector<vector<double>> res;
	res.push_back(A);
	res.push_back(B);
	res.push_back(C);
	return res;
}

vector<vector<double>> getAlphaBeta(){
	vector<double> alpha;
	vector<double> beta;

	vector<vector<double>> v = getABC();
	vector<double> a = v[0];
	vector<double> b = v[1];
	vector<double> c = v[2];
	vector<double> f = get_f();

#ifdef PRINT_MATRIX
	FILE *f_mat;
	f_mat = fopen("matrix.txt","w");
	fprintf(f_mat,"%9f\t%9f\t = %9f\n",b[1],c[1],f[1]);
	for(int i=2; i<=n-1; i++)
		fprintf(f_mat,"%9f\t%9f\t%9f\t = %9f\n",a[i],b[i],c[i],f[i]);	
	fprintf(f_mat,"%9f\t%9f\t = %9f\n",a[n],b[n],f[n]);
	fclose(f_mat);
#endif	
	
	alpha.push_back(0);
	beta.push_back(0);

	alpha.push_back(0);
	alpha.push_back( -c[1]/b[1] );
	beta.push_back(0);
	beta.push_back( f[1]/b[1] );
	for(int i = 2 ; i <= n-1 ; i++){
		alpha.push_back( -c[i]/(a[i]*alpha[i]+b[i]) );
		beta.push_back( (f[i]-a[i]*beta[i])/(a[i]*alpha[i]+b[i]) );
	}

	vector<vector<double>> res;
	res.push_back(alpha);
	res.push_back(beta);
	return res;
}

vector<double> getXbyProgonkaMethod() {
	vector<double> x(n+1);
	
	vector<vector<double>> v = getAlphaBeta();
	vector<double> alpha = v[0];
	vector<double> beta = v[1];
	vector<double> f = get_f();

	v = getABC();
	vector<double> a = v[0];
	vector<double> b = v[1];
	
	x[n] = ( (f[n] - a[n]*b[n])/(b[n] + a[n]*alpha[n]) );

	for( int i=n ; i>=2 ; i-- )
		x[i-1] = alpha[i]*x[i] + beta[i];

#ifdef PRINT_X_PROGONKA
	FILE *f_x;
	f_x = fopen("x_progonka.txt","w");
	for(vector<double>::iterator iter = x.begin()+1 ; iter!=x.end(); iter++)
		fprintf(f_x,"%9f\n",*iter);
	fclose(f_x);
#endif

	return x;
}

int _tmain(int argc, _TCHAR* argv[])
{
	vector<double> x = getXbyProgonkaMethod();
	return 0;
}


// SLAE.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <iostream>
#include <vector>
using namespace std;

#define h 1.0/n
#define PRINT_X_PROGONKA
#define PRINT_MATRIX
#define PRINT_X_ZEIDEL
#define PRINT_X_FASTEST_DESCENT
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
	
	x[n] = ( (f[n] - a[n]*beta[n])/(b[n] + a[n]*alpha[n]) );

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

double max(vector<double> &v) {
	double max = v[1];
	for(vector<double>::iterator iter = v.begin()+1 ; iter!=v.end(); iter++)
		if(*iter > max)
			max = *iter;
	return max;
}

class threeDiagMat {
	vector<double> a;
	vector<double> b;
	vector<double> c;
	int size;
public:
	threeDiagMat(vector<vector<double>> mat) {
		a = mat[0];
		b = mat[1];
		c = mat[2];
		size = b.size()-1;
	}

	double operator()(int i,int j) {
		if(i==1) {
			if(j==1)
				return b[1];
			else if(j==2)
				return c[1];
			else
				return 0;
		}
		else if(i==n) {
			if(j==n-1)
				return a[n];
			else if(j==n)
				return b[n];
			else
				return 0;
		}
		else {
			if(j==i-1)
				return a[i];
			else if(j==i)
				return b[i];
			else if(j==i+1)
				return c[i];
			else
				return 0;
		}
	}

	vector<double> operator *(vector<double> v) {
		vector<double> res(size+1,0);
		for(int i=1 ; i<=size ; i++)
			for(int i2=1 ; i2<=size ; i2++)
				res[i] += (*this)(i,i2)*v[i2];
		return res;
	}
};

vector<double> calcNeviazka(threeDiagMat &mat, vector<double> &x, vector<double> &f) {
	vector<double> r(n+1);
	r[1] = (mat(1,1)*x[1] + mat(1,2)*x[2]) - f[1];
	for(int i = 2 ; i<=n-1 ; i++) {
		double Ay_i = mat(i,i-1)*x[i-1] + mat(i,i)*x[i] + mat(i,i+1)*x[i+1];
		r[i] = Ay_i - f[i];
	}	
	r[n] = (mat(n,n-1)*x[n-1] + mat(n,n)*x[n]) - f[n];
	return r;
}

vector<double> getXbyZeidel(double eps) {
	threeDiagMat mat(getABC());
	vector<double> f = get_f();

#ifdef PRINT_X_ZEIDEL
	FILE *f_x;
	f_x = fopen("x_zeidel.txt","w");
#endif

	vector<double> x(n+1,0);
	x = f;
	vector<double> r = calcNeviazka(mat,x,f);
	double max_r = max(r);
	int count = 0;
	while( max_r > eps ) {
		for(int i = 1 ; i<=n ; i++) {
			double term1 = 0;
			for(int j = 1 ; j <= i-1 ; j++)
				term1 += mat(i,j)/mat(i,i)*x[j];
			double term2 = 0;
			for(int j = i+1 ; j <= n ; j++)
				term2 += mat(i,j)/mat(i,i)*x[j];
			x[i] = -term1 - term2 + f[i]/mat(i,i);
		}
		r = calcNeviazka(mat,x,f);
		max_r = max(r);
		count++;
#ifdef PRINT_X_ZEIDEL
		fprintf(f_x,"%f\n",max_r);
#endif
	}

#ifdef PRINT_X_ZEIDEL
	for(vector<double>::iterator iter = x.begin()+1 ; iter!=x.end(); iter++)
		fprintf(f_x,"%9f\n",*iter);	
	fprintf(f_x,"%i\n",count);
	fprintf(f_x,"%f\n",max_r);
	fclose(f_x);
#endif

	return x;
}

double operator *(vector<double> &v1, vector<double> &v2) {
	int len = v1.size();
	double res = 0;
	for(int i=1 ; i<len ; i++)
		res += v1[i] * v2[i]; 
	return res;
}

vector<double> operator -(vector<double> &v1, vector<double> &v2) {
	int size = v1.size();
	vector<double> res(size,0);
	for(int i=1 ; i < size ; i++)
		res[i] = v1[i] - v2[i]; 
	return res;
}

vector<double> operator *(double p, vector<double> &v1) {
	int size = v1.size();
	vector<double> res(size,0);
	for(int i=1 ; i<size ; i++)
		res[i] = p * v1[i]; 
	return res;
}

vector<double> getXbyFastestDescent(double eps) {
	threeDiagMat mat(getABC());
	vector<double> f = get_f();

#ifdef PRINT_X_FASTEST_DESCENT
	FILE *f_x;
	f_x = fopen("x_fastest_descent.txt","w");
#endif

	vector<double> x(n+1,0);
	vector<double> r = calcNeviazka(mat,x,f);
	double max_r = max(r);
	int count = 0;
	while( max_r > eps ) {
		x = x - ( (r*r)/((mat*r)*r) ) * r;
		r = calcNeviazka(mat,x,f);
		max_r = max(r);
		count++;
#ifdef PRINT_X_FASTEST_DESCENT
		fprintf(f_x,"%f\n",max_r);
#endif
	}

#ifdef PRINT_X_FASTEST_DESCENT
	for(vector<double>::iterator iter = x.begin()+1 ; iter!=x.end(); iter++)
		fprintf(f_x,"%9f\n",*iter);	
	fprintf(f_x,"%i\n",count);	
	fprintf(f_x,"%f\n",max_r);
	fclose(f_x);
#endif

	return x;
}

int _tmain(int argc, _TCHAR* argv[])
{
	getXbyProgonkaMethod();
	getXbyZeidel(h*h*h);
	getXbyFastestDescent(h*h*h);
	return 0;
}


#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>

#include <math.h>

#define EPS 1e-4

using namespace std;

typedef vector<double> vd;
typedef vector<vd> vdd;

const double max_m() 
{
    return 255; //tp
}

double ecm(const vd orig, const vd recup) 
{
    double e=0.0;
    
    if(orig.size()!=recup.size())
    {
        cerr << "Error-ECM: distintas longitudes de vectores" << endl;
        exit(1);
    }
    unsigned int n = orig.size();

    for(unsigned int i=0; i<n; i++)
    {
        double t = orig[i]-recup[i];
        e+=(t*t);
    }
    e=e/n;
    return e;
}


double psnr(const vd orig, const vd recup)
{
    return 10 * log10((max_m()*max_m())/ecm(orig, recup));
}
vector<double> frecuencias;

vector<double> muestreo;

vector<double> lecturas;

vector<vector<double> > MSombrero;

vector<vector<double> > M;

vector<double> gauss(vector<vector<double> > &mat, vector<double> y)
{
	int n = mat.size();
	assert(mat.size()==n&&mat[0].size()==n&&y.size()==n);
	vector<vector<double> > sistema = mat;
	forn(i,n)
		sistema[i].push_back(y[i]);
	vector<int> perm(n);
	forn(i,n)
		perm[i] = i;
	for(int i = 0;i < n-1; i++)
	{
		if(abs(sistema[i][i]-0)<EPS)
		{
			for(int j=i+1;j<n;j++)
			{
				if(abs(sistema[j][i]-0)>EPS)
				{
					swap(sistema[i],sistema[j]);
					swap(perm[i],perm[j]);
				}
			}
		}
		for(int j = i+1;j<n;j++)
		{
			double db = sistema[j][i]/sistema[i][i];
			for(int t=0;t<n+1;t++)
				sistema[j][t] -= db*sistema[i][t];
		}
	}
	vector<double> aux(n,0);
	for(int i=n-1;i>=0;i--)
	{
		int a = sistema[i][n];
		int b = 0;
		for(int j=i+1;j<n;j++)
			b += sistema[i][j]*aux[j];
		aux[i] = (a-b)/sistema[i][i];
	}
	vector<double> sol(n,0);
	for(int i=0;i<n;i++)
		sol[perm[i]] = aux[i];
	return sol;
}

void generarMatrizDCT(int n, double _max)
{
    frecuencias = vector<double> (n);
    muestreo = vector<double> (n);
    MSombrero=vector<vector<double> >(n);
    M=vector<vector<double> >(n);
    
    for(int i = 0; i< n; i++)
    {
        MSombrero[i]=vector<double>(n);
        M[i]=vector<double>(n);
    }

    double k = (M_PI/(double)n);
    for(int i = 0; i < n; i++)
    {
        frecuencias[i]=i;
        muestreo[i]=k*(1.0+((i/2)));
    }
    double sq1n=sqrt(1.0/n);
    double sq2n=sqrt(2.0/n);
#pragma omp parallel for
    for(int i = 0; i< n ; i++)
    {
        #pragma omp parallel for
        for(int j = 0; j< n ; j++)
        {
            double v;
            v=cos(frecuencias[i]*muestreo[j]);
            if(i==0) v = v*sq1n;
            else     v = v*sq2n;
            MSombrero[i][j]=v;
            M[i][j]=floor((max_m()*v + 1)/2);
        }
    }

    double calc=0.0;
    for(int i=0;i<n;i++) 
        for(int j=0; j<n;j++)
            calc+=MSombrero[i][j];
    double constante=6.6274170;
    if(abs(calc-constante)>EPS)
    {
        cerr << "DIFIEREN!" << endl;
        fprintf(stderr, "%.7lf vs %.7lf\n", calc,constante);
    }

}

vector<double> obtenerY(int n)
{
    vector<double> y(n);
#pragma omp parallel for
    for(int i = 0; i<n; i++)
    {
        double q=0.0;
        for(int j =0 ; j<n;j++)
        {
            q+=MSombrero[i][j]*lecturas[j];
        }
        y[i]=q;
        
    }
    return y;
}

double ruido1(double x)
{
    return x*sin(x)*M_PI*0.3;
}

void modificar(vector<double> &y, int n)
{
    #pragma omp parallel for
    for(int i = 0 ; i < n ; i++)
    {
        y[i]=ruido1(y[i]);
    }
    return;
}


int main(int argc, char* argv[])
{

    int n;

    cin >> n;
    lecturas = vector<double> (n);
    double _max = 0.0;
    for(int i = 0; i < n; i++)
    {
        cin >> lecturas[i];        
        _max=(_max>abs(lecturas[i])?_max:abs(lecturas[i]));
    }
    generarMatrizDCT(n, _max);
    vector<double> y = obtenerY(n);
    modificar(y,n);
    
    return 0;
}

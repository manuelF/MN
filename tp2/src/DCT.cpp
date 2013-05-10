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

#define EPS 1e-8

#define double long double

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

vector<vector<double> > mbmt(vector<vector<double> > &b) //Dado B calculo MBM^t como pide en el apendice A1
{
    int n = M.size();
    assert((int)M[0].size() == n && (int)b.size() == n &&(int) b[0].size() == n);
    vector<vector<double> > aux(n,vector<double>(n,0)); // la creo con todos ceros la matriz
    for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
    {
        for(int t=0;t<n;t++)
            aux[i][j] += M[i][t]*b[t][j];
    }
    vector<vector<double> > sol(n,vector<double>(n,0));
    for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
    {
        for(int t=0;t<n;t++)
            sol[i][j] = aux[i][t]*M[j][t]; // (j,t) y no (t,j) porque es la transpuesta
    }
    return sol;
}

vector<double> gauss(vector<vector<double> > &mat, vector<double> y) // mat * ret = y
{
	int n = mat.size();
	assert((int)mat.size()==n&&(int)mat[0].size()==n&&(int)y.size()==n);
	vector<vector<double> > sistema = mat;
	for(int i=0;i<n;i++)
		sistema[i].push_back(y[i]); // sistema es [mat|y]
	vector<int> perm(n); // el vector de permutacion
	//for(int i=0;i<n;i++)
		//perm[i] = i; // comienzo con la permutacion identidad
	for(int j = 0;j < n-1; j++)
	{
		int mx = j;
		for(int t = j+1; t<n ;t++)
		if(abs(sistema[mx][j])< abs(sistema[t][j]))
			mx = t;
		swap(sistema[mx],sistema[j]);
		for(int i = j+1;i<n;i++)
		{
			double db = sistema[i][j]/sistema[j][j];
			for(int t=j;t<n+1;t++)
				sistema[i][t] -= db*sistema[j][t];
		}
	}
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
			cout << sistema[i][j] << " ";
			cout<< endl;
	}
	vector<double> x(n,0);// los valores de x antes de aplicar la permutacion
	for(int i=n-1;i>=0;i--)
	{
		int y_i = sistema[i][n]; // el que seria el numero en y en [mat|y]
		int a = 0;
		for(int j=i+1;j<n;j++)
			a += sistema[i][j]*x[j];// 
		x[i] = (y_i-a)/sistema[i][i];
	}
	return x;
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

vector<double> f(vector<double> y, int imp) // imp es la implementacion
{
	if(imp == 0)
	{
		// Me quedo con los que son mayores al 20% del mayor
		int max = 0;
		for(int i=0;i<(int)y.size();i++)
		if(y[i]>y[max])
			max = y[i];
		for(int i=0;i<(int)y.size();i++)
		if(5*y[i]<max)
			y[i] = 0;
		return y;
	}
	if(imp == 1)
	{
		// Me quedo con el 10% mas grande y los otros los descarto
		vector<pair<double,int> > aux(y.size());
		for(int i=0;i<(int)y.size();i++)
			aux[i] = make_pair(y[i],i);
		sort(aux.begin(),aux.end());
		reverse(aux.begin(),aux.end()); // ordeno de mayor a menor
		for(int i=(int)y.size()/10;i<(int)y.size();i++)
			aux[i].first = 0; // los que se pasan del 10% los convierto en cero
		for(int i=0;i<(int)y.size();i++)
			y[aux[i].second] = aux[i].first;
		return y;
	}
	// podemos hacer mas implementaciones de ser necesario
	return y;
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
    y = f(y,0);
//  for(int i=0;i<(int)y.size();i++)
    	//cout << y[i] << endl;
    vector<double> x = gauss(M,y);
	//for(int i=0;i<(int)x.size();i++)
		//cout << x[i] << endl;
	/*vector<vector<double> > vec(3,vector<double>(3,0));
	vec[0][0] = 1.;
	vec[1][2] = 2.;
	vec[2][1] = 2.;
	vec[2][2] = 2.;
	vector<double> aux(3);
	aux[0] = 1;
	aux[1] = 4;
	aux[2] = 12;
	aux = gauss(vec,aux);
	//cout << aux[0] <<" "<<aux[1]<<" " <<aux[2]<<endl;
	/*
	 * 1 0 0  1   1
	 * 0 0 2  2   8
	 * 0 2 2  4   12
	 */ 
    return 0;
}

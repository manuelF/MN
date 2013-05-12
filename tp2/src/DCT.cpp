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

#include <omp.h>
#include <math.h>

#include <random>

#define EPS 1e-8

#define double long double

using namespace std;


const double max_m()
{
    return 255; //tp
}

double ecm(const vector<double>& orig, const vector<double>& recup)
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
    e=e/(double)n;
    return e;
}


double psnr(const vector<double>&  orig, const vector<double> &recup)
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

vector<double> gauss(const vector<vector<double> > &mat, const vector<double> &y) // mat * ret = y
{    
	int n = mat.size();
	assert((int)mat.size()==n&&(int)mat[0].size()==n&&(int)y.size()==n);
	vector<vector<double> > sistema(mat);
	for(int i=0;i<n;i++)
		sistema[i].push_back(y[i]); // sistema es [mat|y]
	for(int j = 0;j < n-1; j++)
	{
		int mx = j;
		for(int t = j+1; t<n ;t++)
		if(abs(sistema[mx][j])< abs(sistema[t][j]))
			mx = t;
		swap(sistema[mx],sistema[j]);
        double m = sistema[j][j];
        
		for(int i = j+1;i<n;i++)
		{
			double db = sistema[i][j]/m;
			for(int t=j;t<n+1;t++)
				sistema[i][t] -= db*sistema[j][t];
		}
	}
	vector<double> x(n,0);// los valores de x antes de aplicar la permutacion
	for(int i=n-1;i>=0;i--)
	{
		double y_i = sistema[i][n]; // el que seria el numero en y en [mat|y]
		double a = 0;
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
        muestreo[i]=k*(i+((1./2.)));
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
}

vector<double> obtenerY(const vector<double>& l, int n)
{
    vector<double> y(n);
    
//#pragma omp parallel for
    for(int i = 0; i<n; i++)
    {
        double q=0.0;
        for(int j =0 ; j<n;j++)
        {
            q+=M[i][j]*l[j];
        }
        y[i]=q;

    }
    return y;
}


void f(vector<double> &y, int imp) // imp es la implementacion
{
	if(imp == 0)
	{
		// Me quedo con los que son mayores al 20% del mayor
		int max = y[0];
		for(int i=0;i<(int)y.size();i++)
            if(abs(y[i])>abs(max))
                max = y[i];
		for(int i=0;i<(int)y.size();i++)
            if(5*abs(y[i])<abs(max))
                y[i] = 0;		
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
	}
    if(imp == 2)
	{
        std::default_random_engine generator;
        std::normal_distribution<double> distribution(0.0,10.0);
		for(int i=0;i<(int)y.size();i++)
        {
            double randval= distribution(generator);
            y[i]=y[i]+randval;
        }
	}
    if(imp == 3)
	{
        for(int i =0; i<(int)y.size();i++)
        {
            y[i]=(double)y[i]+sin(i);
        }
	}
	// podemos hacer mas implementaciones de ser necesario
	return;
}

void dump(char* name, const vector<double>& src)
{
    int n = src.size();
    FILE* fileptr = fopen(name,"w");
    for(int i=0; i< n; i++)
        fprintf(fileptr,"%.6Lf\n",src[i]);
    fclose(fileptr);
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
    
    dump("orig",lecturas);

    generarMatrizDCT(n, _max);
    vector<double> l (lecturas);
    f(l,3);
    dump("mod",l);
    vector<double> y = obtenerY(l,n);
    dump("recovered",y);
    f(y,0);
    //dump("mod",y);
    vector<double> x = gauss(M,y);
    //dump("recovered",x);
    
    cerr<< "PNSR: " << psnr(lecturas,x) << endl;
    return 0;
}

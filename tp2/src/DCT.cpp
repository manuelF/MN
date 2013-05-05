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

#include <math.h>

#define EPS 1e-4

using namespace std;

typedef vector<double> vd;
typedef vector<vd> vdd;

double max_m()
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

void generarMatrizDCT(int n)
{
    frecuencias = vector<double> (n);
    muestreo = vector<double> (n);
    MSombrero=vector<vector<double> >(n);
    
    for(int i = 0; i< n; i++)
        MSombrero[i]=vector<double>(n);

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


int main(int argc, char* argv[])
{

    int n;

    cin >> n;
    lecturas = vector<double> (n);
    for(int i = 0; i < n; i++)
    {
        cin >> lecturas[i];        
    }
    generarMatrizDCT(n);
    vector<double> y = obtenerY(n);

    return 0;
}

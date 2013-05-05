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

vector<vector<double> > M;

void generarMatrizDCT(int n)
{
    frecuencias = vector<double> (n);
    muestreo = vector<double> (n);
    M=vector<vector<double> >(n);
    
    for(int i = 0; i< n; i++)
        M[i]=vector<double>(n);

    double k = (M_PI/(double)n);
    for(int i = 0; i < n; i++)
    {
        cin >> lecturas[i];        
        frecuencias[i]=i;
        muestreo[i]=k*(1.0+((i/2)));
    }
#pragma omp parallel for
    for(int i = 0; i< n ; i++)
    {
        #pragma omp parallel for
        for(int j = 0; j< n ; j++)
        {
            M[i][j]=cos(frecuencias[i]*muestreo[j]);
        }
    }

    double calc=0.0;
    for(int i=0;i<n;i++) 
        for(int j=0; j<n;j++)
            calc+=M[i][j];
    double constante=256.0;
    cerr << calc<<endl;
    if(abs(calc-constante)>EPS)
        cerr << "DIFIEREN!" << endl;

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

    return 0;
}

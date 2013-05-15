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


#define GAUSSIAN_NOISE 0
#define SIN_NOISE 1
#define ADDITIVE_NOISE 2

#define ZERO_FILTER 0
#define EXPONENTIAL_FILTER 1
#define AVERAGER_FILTER 2

//typedef long double double;

using namespace std;


const double max_m()
{
    return 255.0; //tp
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
    double _max = *max_element(orig.begin(), orig.end());
    
    return 10 * log10((_max*_max)/ecm(orig, recup));
}
vector<double> frecuencias;
vector<double> muestreo;
vector<double> lecturas;
vector<double> ruido;
vector<vector<double> > MSombrero;
vector<vector<double> > M;

vector<vector<double> > mbmt(vector<vector<double> > &b) //Dado B calculo MBM^t como pide en el apendice A1
{
    int n = M.size();
    assert((int)M[0].size() == n && (int)b.size() == n &&(int) b[0].size() == n);
    vector<vector<double> > aux(n,vector<double>(n,0)); // la creo con todos ceros la matriz
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            for(int t=0;t<n;t++)
                aux[i][j] += M[i][t]*b[t][j];
        }
    }
    vector<vector<double> > sol(n,vector<double>(n,0));
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            for(int t=0;t<n;t++)
            {
                sol[i][j] = aux[i][t]*M[j][t]; // (j,t) y no (t,j) porque es la transpuesta
            }
        }
    }
    return sol;
}

vector<double> antitransformar(const vector<double> &y) // mat * ret = y
{   
    const vector<vector<double> > &mat = M;
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

void generarMatrizDCT(int n)
{
    frecuencias = vector<double> (n);
    muestreo = vector<double> (n);
    MSombrero=vector<vector<double> >(n, vector<double>(n));
    M=vector<vector<double> >(n, vector<double>(n));


    double k = (M_PI/(double)n);
    for(int i = 0; i < n; i++)
    {
        frecuencias[i]=i;
        muestreo[i]=k*(i+((1./2.)));
    }
    double sq1n=sqrt(1.0/n);
    double sq2n=sqrt(2.0/n);
    double _max = max_m();
    double v = sq1n;
#pragma omp parallel for
    for(int i = 0; i< n ; i++)
    {
        double freq = frecuencias[i];
        #pragma omp parallel for
        for(int j = 0; j< n ; j++)
        {
            if(i>0) v=sq2n;
            double theta=v*cos(freq*muestreo[j]);
            MSombrero[i][j]=theta;
            M[i][j]=floor((_max*theta + 1)/2);
        }
    }
}

vector<double> transformar(const vector<double>& l)
{
    int n = (int) l.size();
    vector<double> y(n);

#pragma omp parallel for
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


void generarRuido(vector<double> &y, int imp)
{
    if(imp == GAUSSIAN_NOISE)
    {
        ruido = vector<double>(y.size());
        std::default_random_engine generator;
        std::normal_distribution<double> distribution(0.0,10.0);
		for(int i=0;i<(int)y.size();i++)
        {
            double randval= distribution(generator);
            ruido[i]=randval;
            y[i]=y[i]+ruido[i];
        }
	}
    if(imp ==  SIN_NOISE)
	{
        ruido = vector<double>(y.size());
        #pragma omp parallel for
        for(int i =0; i<(int)y.size();i++)
        {
            ruido[i]=50*sin(i);
            y[i]=(double)y[i]+ruido[i];
        }
	}
	// podemos hacer mas implementaciones de ser necesario
	return;

}

void filtrarRuido(vector <double> &y, int imp)
{
    const int n = (int) y.size();
    const int contorno = 20;
    const double diffmax = 2.0;
    const int onepercent = n/100;
    const int startfilter = 5*onepercent;
    vector<double> replace(y);
	if(imp == ZERO_FILTER)
	{
	    double mx = abs(y[startfilter]);
	    for(int i=startfilter;i<n;i++) mx = max(mx, abs(y[i]));
        for(int i=startfilter;i<n;i++)
        {
            if(abs(y[i])*diffmax>mx)
            {
                for(int j=max(startfilter,i-contorno);j<min(n,i+contorno);j++)
                    replace[j] = 0;
            }        
        }
	}
	if(imp == EXPONENTIAL_FILTER)
	{
		double mx = abs(y[startfilter]);
	    for(int i=startfilter;i<n;i++) mx = max(mx, abs(y[i]));
        for(int i=startfilter;i<n;i++)
        {
            if(abs(y[i])*diffmax>mx)
            {
                 for(int j=max(0,i-contorno);j<min(n,i+contorno);j++)
                 {
                        double decay = (1.0/pow(1.2,contorno-abs(j-i)+1));
                        replace[j] = y[j]*decay;
                 }
        
            }
        }       
	}
    if(imp == AVERAGER_FILTER)
    {
       for(int i=startfilter+2;i<n-2;i++)
        {
            replace[i]= .15*y[i-2]+
                        .20 *y[i-1]+
                        .30 *y[i]  +
                        .20 *y[i+1]+
                        .15*y[i+2];
        }
    }
    for (int i = 0; i< (int) y.size(); i++)
        y[i]=replace[i];

}
void dump(string name, const vector<double>& src)
{
    int n = src.size();
    FILE* fileptr = fopen(name.c_str(),"w");
    for(int i=0; i< n; i++)
        fprintf(fileptr,"%.6lf\n",src[i]);
    fclose(fileptr);
}

void procesar1D()
{
    int n;
    cin >> n;
    lecturas = vector<double> (n);
    double _max=0.0;
    for(int i = 0; i < n; i++)
    {
        cin >> lecturas[i];
        _max=(_max>abs(lecturas[i])?_max:abs(lecturas[i]));
    }


    generarMatrizDCT(n);
    vector<double> l (lecturas);
    //generarRuido(l,GAUSSIAN_NOISE);
    generarRuido(l,SIN_NOISE);

    vector<double> q = transformar(lecturas); //original

    vector<double> y = transformar(l); //transformada

    
    dump("orig",q);
    
    filtrarRuido(y,EXPONENTIAL_FILTER);
    filtrarRuido(y,AVERAGER_FILTER );

//    dump("mod",y);
    vector<double> x = antitransformar(y);
//    dump("recovered",x);
    cerr<< "PNSR: " << psnr(lecturas,x) << endl;
}



void procesar2D()
{
    string magic;
    cin >> magic;
    if(magic!="P5")
    {
        cerr << "ERROR: PGM tipo " << magic << " no implementado, solo P5"<< endl;
        exit(1);
    }
    int x, y;
    cin >> x >> y;
    vector<vector<double>> _img(y,vector<double>(x));
    
    for(int j = 0; j < y; j++)
    {
        for(int i = 0; i < x; i++)
        {
            unsigned char read;
            if(scanf("%c",&read)!=1)
            {
    //              cerr << "ERROR: PGM no se puede leer correctamente" << endl;
            }
            _img[j][i]=(double) read;
        }
    }


    return;

}

int main(int argc, char* argv[])
{
    switch(argc)
    {
        case 1:
            procesar1D();
            break;
        case 2:
            if(atoi(argv[1])==1)
            {
            
                procesar2D();
            }
            else 
            {   
                if(atoi(argv[1])==0)
                {
                    procesar1D();
                }
            }

            break;
        default:
            break;

    }
   return 0;

}

#include<algorithm>
#include<cmath>
#include<cstdio>
#include<iostream>
#include<map>
#include<queue>
#include<set>
#include<sstream>
#include<string>
#include<vector>
#include "TFloat.h"
#include <math.h>

using namespace std;

TFloat::TFloat()
{
	_valor = 0.0;
	set_precision(52);
	recortar();
}

TFloat::TFloat(size_t t)
{
	_valor = 0.0;
	set_precision(t);
	recortar();
}

TFloat::TFloat(int i, size_t t)
{
	_valor = (double) i;
	set_precision(t);
	recortar();
}

TFloat::TFloat(float f, size_t t)
{
	_valor = (double) f;
	set_precision(t);
	recortar();
}

TFloat::TFloat(double f, size_t t)
{
	_valor = f;
	set_precision(t);
	recortar();
}

double TFloat::dbl() const
{
	return _valor;
}

void  TFloat::operator =(const TFloat& f)
{
	_valor = f._valor;
	_t = f._t;
}

void  TFloat::operator =(const double& f)
{
	_valor = f;
	recortar();
}

bool  TFloat::operator ==(const TFloat& f) const
{
	return (_valor == f._valor);
}

TFloat TFloat::operator + (const TFloat& f) const
{
	double result = _valor + f._valor;
    return TFloat(result, _t);
}

TFloat TFloat::operator - (const TFloat& f) const
{
	double result = _valor - f._valor;
    return TFloat(result, _t);
}

TFloat TFloat::operator * (const TFloat& f) const
{
	double result = _valor * f._valor;
    return TFloat(result, _t);
}

TFloat TFloat::operator / (const TFloat& f) const
{
	double result = _valor / f._valor;
    return TFloat(result, _t);
}

TFloat TFloat::exponencial() const
{
	double result = exp(_valor);
    return TFloat(result, _t);
}

TFloat TFloat::operator+ (const double& f) const
{
	return this->operator+(TFloat(f,_t));
}

TFloat TFloat::operator- (const double& f) const
{
	return this->operator-(TFloat(f,_t));
}

TFloat TFloat::operator* (const double& f) const
{
	return this->operator*(TFloat(f,_t));
}

TFloat TFloat::operator/(const double& f) const
{
	return this->operator/(TFloat(f,_t));
}


void TFloat::recortar()
{
	bitset<64> * bits;
	bits = (bitset<64>*) &_valor;

	for(size_t i = 0 ; i < (52 - _t); i++ )
	{
		(*bits)[i] = 0;
	}
}

vector<TFloat> valores;
int n, t;
const double epsilon =  1.e-9;
int maximoIteraciones = 100;
TFloat pot(TFloat base, TFloat exp)
{
    double result = pow(base.dbl(),exp.dbl());
    return TFloat(result, t);
}

TFloat log(TFloat d)
{
    double result = log(d.dbl());
    return TFloat(result, t);
}

TFloat M(TFloat s)
{
    TFloat res = TFloat(0.,t);
    for(int i=0;i<n;i++)
        res = res + pot(valores[i],s);
    res = res/TFloat(n,t);
    return res;
}

TFloat MSombrero(TFloat s)
{
    TFloat res = TFloat(0.,t);
    for(int i=0;i<n;i++)
        res = res + pot(valores[i],s)*log(valores[i]);
    res = res/TFloat(n,t);
    return res;
}

TFloat R(TFloat s)
{
    return MSombrero(s)/M(s);
}

TFloat Lambda(TFloat beta)
{
    vector<TFloat> logvals(n);
    vector<TFloat> xbeta(n);
    TFloat sumalogvals(0.0,t);
    TFloat sumaxbeta(0.0, t);
    TFloat sumalogxbeta(0.0, t);
    for(int i=0; i<n; i++)
    {
        logvals[i]=log(valores[i]);
        xbeta[i]=pot(valores[i],beta);
    }
    for(int i=0; i<n; i++)
    {
        sumalogvals= sumalogvals + logvals[i];
        sumaxbeta= sumaxbeta + xbeta[i];
        sumalogxbeta= sumalogxbeta+ logvals[i]*xbeta[i];
    }
    TFloat lambda = pot(beta*(((sumalogxbeta)/(sumaxbeta))-((sumalogvals)/(TFloat(n,t)))),TFloat(-1,t));
    return lambda;

}

TFloat Sigma(TFloat beta, TFloat lambda)
{
    TFloat sigma;
    TFloat xbetasuma(0.0,t);

    for(int i=0; i<n; i++)
    {
        xbetasuma=xbetasuma + pot(valores[i],beta);
    }
    sigma = pot((xbetasuma/(TFloat(n,t)*lambda)),(TFloat(1.0,t)/beta));

    return sigma;
}

TFloat f(TFloat beta)
{
    return M(beta*2)/(M(beta)*M(beta)) - beta * (R(beta)-R(0))  - 1.;
}


TFloat regulaFalsi(TFloat beta, TFloat beta2, int& iteraciones)
{
    iteraciones=0;
    TFloat beta3;
    while(iteraciones<maximoIteraciones && abs(beta.dbl()-beta2.dbl())>epsilon)
    {
        iteraciones++;
        beta3 = beta2 - f(beta2)*(beta2-beta)/(f(beta2)-f(beta));
        if(f(beta).dbl()*f(beta3).dbl()<0)
            beta2 = beta3;
        else
        {
            beta = beta2;
            beta2 = beta3;
        }
    }
    return beta;
}

void uso()
{
    cout << "\"./Newton\" precision = 52, iteraciones maximas = 10, Beta entre 0.0 y 2.0"<<endl;
    cout << "\"./Newton <t>\" precision = t, iteraciones maximas = 10, Beta entre 0.0 y 2.0"<<endl;
    cout << "\"./Newton <t> <n> <b1> <b2>\" precision = t, iteraciones maximas= n, Beta entre b1 y b2"<<endl;
}

int main(int argc, char* argv[])
{
    TFloat beta = TFloat(1.,52);
    TFloat beta2 = TFloat(8.,52);
    switch(argc)
    {
        case 1:
            t=52;
            break;
        case 2:
            t=atoi(argv[1]);
            break;
        case 5:
            t=atoi(argv[1]);
            maximoIteraciones=atoi(argv[2]);
            beta=TFloat(atof(argv[3]),t);
            beta2=TFloat(atof(argv[4]),t);
            break;
        default:
            uso();
            return 1;
            break;
    }
    cin >> n;
    valores.resize(n);
    double db;
    for(int i=0;i<n;i++)
    {
        cin >> db;
        valores[i] = TFloat(db,t);
    }
    int iteraciones = 0;
    TFloat outBeta;
    TFloat outLambda;
    TFloat outSigma;

    outBeta = regulaFalsi(beta, beta2, iteraciones);
    outLambda = Lambda(outBeta);
    outSigma = Sigma(outBeta, outLambda);

    cout << outBeta.dbl() << " " << outLambda.dbl() << " " << outSigma.dbl() << endl;
    return 0;
}

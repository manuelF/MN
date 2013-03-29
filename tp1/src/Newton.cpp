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
int maximoIteraciones = 10;
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

TFloat MPrima(TFloat s)
{
    /// Resulta ser que Mprima es igual a Msombrero
    return MSombrero(s);
}

TFloat MSombreroPrima(TFloat s)
{
    TFloat res = TFloat(0.,t);
    for(int i=0;i<n;i++)
        res = res + pot(valores[i],s)*log(valores[i])*log(valores[i]);
    res = res/TFloat(n,t);
    return res;
}

TFloat R(TFloat s)
{
    return MSombrero(s)/M(s);
}

TFloat RPrima(TFloat s)
{
    return (MSombreroPrima(s)*M(s)-MSombrero(s)*MPrima(s))/(M(s)*M(s));
}

TFloat f(TFloat beta)
{
    return M(beta*2)/(M(beta)*M(beta)) - beta * (R(beta)-R(0))  - 1.;
}


TFloat fprima(TFloat beta)
{

    /** Mprima(beta*2)*2 = Derivada de M(beta*2)
     *  M(beta)*Mprima(beta)*2 = Derivada de M(beta)*M(beta) o M(beta)^2
     *  Derivada de M(beta*2)/M(beta)^2 =
     *      (Mprima(beta*2)*M(beta)*M(beta)*2 - M(beta*2)*M(beta)*Mprima(beta)*2)/(M(beta)*M(beta)*M(beta)*M(beta))
     *  Derivada de beta * (R(beta) - R(0)) =
     *      beta * (Rprima(beta)) + R(beta)-R(0) //// Rprima(0) = 0 porque es una constante?
     */
  

return  (MPrima(beta*2)*M(beta)*M(beta)*2 - M(beta*2)*M(beta)*MPrima(beta)*2)/(M(beta)*M(beta)*M(beta)*M(beta))
            - (beta * RPrima(beta) + R(beta) - R(0));
}


//TFloat f(TFloat x) {return x*x - 2.0;}
//TFloat fprima(TFloat x){ return TFloat(2.)*x;}


TFloat newton(TFloat beta, TFloat beta2, int& iteraciones)
{
    iteraciones=0;
    
    while(iteraciones<maximoIteraciones && abs(beta.dbl()-beta2.dbl())>epsilon)
    {
        iteraciones++;
        beta2 = beta;
        beta = beta - f(beta)/fprima(beta);
    }
    return beta;
}

TFloat secante(TFloat p0, TFloat p1, int& iteraciones)
{
    TFloat q1=f(p1);
    TFloat q0=f(p0);
    TFloat pnew=p0;
   
    iteraciones=0;
    
    while(iteraciones<maximoIteraciones && (abs(pnew.dbl()-p1.dbl())>epsilon))
    {
        iteraciones++;
        cout << "p0: " << p0.dbl()<< " p1: " << p1.dbl() << " q0: " << q0.dbl() << " q1: " << q1.dbl() << endl;
        pnew = p1 - ((q1*(p1-p0))*(TFloat(1.,52)/(q1-q0)));          

        p0=p1; p1=pnew;
        q0=q1; q1=f(pnew);

    }
    return pnew;
 
}

void uso()
{
    cout << "\"./Newton\" precision = 52, iteraciones maximas = 10, Beta entre 0.0 y 1.0"<<endl;
    cout << "\"./Newton <t>\" precision = t, iteraciones maximas = 10, Beta entre 0.0 y 1.0"<<endl;
    cout << "\"./Newton <t> <n> <b1> <b2>\" precision = t, iteraciones maximas= n, Beta entre b1 y b2"<<endl;
}

int main(int argc, char* argv[])
{
    {
    }
    TFloat beta = TFloat(2.,52);
    TFloat beta2 = TFloat(0.,52);
    switch(argc)
    {
        case 1:
            t=52;
            break;
        case 2:
            t=atoi(argv[1]);
            break;
        case 3:
            uso();
            return 1;
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
    TFloat out;

    out = secante(beta, beta2, iteraciones);
    cout << "Secante: " << out.dbl() << " "<< iteraciones << endl;
    out = newton(beta, beta2, iteraciones);
    cout << "Newton: " << out.dbl() << " "<< iteraciones << endl;
    return 0;
}

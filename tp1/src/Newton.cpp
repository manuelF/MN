#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <vector>
#include <map>
#include "TFloat.h"
#include <math.h>

#ifdef __gnu_linux__
    #include "timetools.h"
#endif
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

bool TFloat::operator<(const TFloat &f) const
{
	return this->dbl()<f.dbl();
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
const double epsilon =  1.e-4;
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

map<TFloat,TFloat> mapaMSombreroPrima;

TFloat MSombreroPrima(TFloat s)
{
    if(mapaMSombreroPrima.find(s)!=mapaMSombreroPrima.end())
	return mapaMSombreroPrima[s];
    TFloat res = TFloat(0.,t);
    for(int i=0;i<n;i++)
        res = res + pot(valores[i],s)*log(valores[i])*log(valores[i]);
    res = res/TFloat(n,t);
    mapaMSombreroPrima[s] = res;
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
    TFloat Mbeta = M(beta);
    TFloat res = M(beta*2)/(Mbeta*Mbeta) - beta * (R(beta)-R(0))  - 1.;    
    return res;
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

    TFloat Mbeta=M(beta);
    TFloat Mbeta2=Mbeta*Mbeta;
        //Aca podemos simplificar varios terminos sacando factores comunes
return  (TFloat(2,t))*(MPrima(beta*2)*Mbeta - M(beta*2)*MPrima(beta))/(Mbeta*Mbeta2)
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

TFloat illinois(TFloat beta, TFloat beta2, int& iteraciones)
{
    iteraciones=0;
    TFloat beta3;
    TFloat fbeta=f(beta);
    TFloat fbeta2=f(beta2);
    for(int i=0;i<10;i++)
    if(fbeta.dbl()>0)
    {
	beta = beta/2.;
	fbeta = f(beta);
    }
    for(int i=0;i<10;i++)
    if(fbeta2.dbl()<0)
    {
    	beta2 = beta2*2;
	fbeta2 = f(beta2);
    }
    if(fbeta.dbl()*fbeta2.dbl()>0)
    {
	cerr << setprecision(15);
        cerr<< "==================fruta follows===================\n";
        cerr<< beta.dbl() << ": " << fbeta.dbl() << "   "<<beta2.dbl() << ": " << fbeta2.dbl() << endl;
    }
    int side = -1;
    while(iteraciones<maximoIteraciones && fabs(beta.dbl()-beta2.dbl())>epsilon)
    {
        iteraciones++;
//r = (fs*t - ft*s) / (fs - ft);
	beta3 = (fbeta*beta2 - fbeta2*beta)/(fbeta-fbeta2);
        //beta3 = beta2 - fbeta2*(beta2-beta)/(fbeta2-fbeta);
        TFloat fbeta3=f(beta3); 
        if(fbeta2.dbl()*fbeta3.dbl()>0)
        {
            beta2 = beta3;
            fbeta2= fbeta3;
            if(side == -1) fbeta=fbeta/2;
	    side = -1;
        }
        else if(fbeta.dbl()*fbeta3.dbl()>0)
        {
            beta = beta3;
            fbeta=fbeta3;
            if(side == 1) fbeta2=fbeta2/2;
	    side = 1;
        }
    }
    cerr << abs(beta2.dbl()-beta.dbl()) <<" " << iteraciones << endl;
    return beta3;
}

void uso()
{
    cout << "\"./Newton\" precision = 52, iteraciones maximas = 10, Beta entre 0.0 y 10.0"<<endl;
    cout << "\"./Newton <t>\" precision = t, iteraciones maximas = 10, Beta entre 0.0 y 10.0"<<endl;
    cout << "\"./Newton <t> <n> <b1> <b2> <m> \" precision = t, iteraciones maximas= n, Beta entre b1 y b2, <m> metodo (0 Newton - 1 RegulaFalsi)"<<endl;
}

int main(int argc, char* argv[])
{
    TFloat beta = TFloat(10.,52);
    TFloat beta2 = TFloat(1.,52);
    int metodo=0; //0 Newton 1 RegulaFalsi
    switch(argc)
    {
        case 1:
            t=52;
            break;
        case 2:
            t=atoi(argv[1]);
            break;
        case 6:
            metodo = atoi(argv[5]);
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
    mapaMSombreroPrima.clear();
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

    cout << setprecision(15);
    #ifdef __gnu_linux__
        timespec startTime;
        timespec endTime;

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &startTime);
    #endif
    if (metodo==0)
    {
        outBeta = newton(beta, beta2, iteraciones);
    }
    else if (metodo==1)
    {
        outBeta = illinois(beta, beta2, iteraciones);
    }
    else if(metodo==2)
    {
    
        for(int i =40;i<100;i++)
        {
            double val = 0.1 * i;
            cout << val << " " << f(TFloat(val,t)).dbl() << endl;
        }
    }

    outLambda = Lambda(outBeta);
    outSigma = Sigma(outBeta, outLambda);
    #ifdef __gnu_linux__
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &endTime);
        const timespec delta = diff(startTime, endTime);
    #endif
    /*
    outBeta = puntofijo(beta, iteraciones);
    outLambda = Lambda(outBeta);
    outSigma = Sigma(outBeta, outLambda);

    cout << "RF: " <<iteraciones <<" iteraciones (de un maximo de "<< maximoIteraciones << ") "<< endl;
    cout << "    - BetaMoño: " << outBeta.dbl() << endl;
    cout << "    - LambdaMoño: " << outLambda.dbl() << endl;
    cout << "    - SigmaMoño: " << outSigma.dbl() << endl;
    */

    /*
    cout << "Newton: " <<iteraciones << " iteraciones (de un maximo de "<< maximoIteraciones << ") "<< endl;
    cout << "    - BetaMoño: " << outBeta.dbl() << endl;
    cout << "    - LambdaMoño: " << outLambda.dbl() << endl;
    cout << "    - SigmaMoño: " << outSigma.dbl() << endl;
    */
    cout << outBeta.dbl() << " " << outLambda.dbl() << " " << outSigma.dbl() << " " << iteraciones;
    #ifdef __gnu_linux__
        cout << " " << (delta.tv_sec * 1000*1000*1000 + delta.tv_nsec)/(1000*1000);
    #endif
    cout << endl;
    return 0;
}

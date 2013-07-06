#include<iostream>
#include<algorithm>
#include<vector>
#include<cstring>
#include<cstdio>
#include<cmath>
#include<cassert>

using namespace std;

vector<vector<double> > input, X, Mx, Xt, testImages, matrizAuxiliar, av; //av = autovectores
vector<int> labels, testLabels;
vector<double> average;

void transpose(vector<vector<double> > &mat)
{
    /** Transposicion de matrices **/
    int n = mat.size();
    int m = mat[0].size();
    matrizAuxiliar.clear();
    matrizAuxiliar.resize(m,vector<double>(n));
	for(int i=0;i<m;i++)
    for(int j=0;j<n;j++)
	{
	    matrizAuxiliar[i][j] = mat[j][i];
    }
    mat = matrizAuxiliar;
	return;
}

vector<vector<double> > mult(vector<vector<double> > &A, vector<vector<double> > &B)
{
    /** Multiplicacion de matrices standard en o(n^3) **/
	int n = A.size();
	int m = B[0].size();
	int t = A[0].size(); // tiene que ser igual a B.size();
    vector<vector<double> > B2 = B;
    transpose(B2);
	vector<vector<double> > res(n,vector<double>(m,0));
	for(int i=0;i<n;i++)
	{
    	for(int j=0;j<m;j++)
        {
            double sum =0.0;
	        for(int k=0;k<t;k++)
            {
		        sum += A[i][k]*B2[j][k];
            }
            res[i][j]=sum;
        }
    }
	return res;
}

void generateX()
{
	int n = input.size();
	int m = input[0].size();
	average.clear();
	average.resize(m,0);
	for(int i=0;i<n;i++)
    	for(int j=0;j<m;j++)
	    	average[j] += input[i][j];
	for(int j=0;j<m;j++)
		average[j] /= (double)n;
    /** Calculo el promedio de cada pixel **/
	X.clear();
	X.resize(n,vector<double>(m));
	for(int i=0;i<n;i++)
    	for(int j=0;j<m;j++)
	    	X[i][j] = (input[i][j]-average[j])/sqrt(n-1);
            /** A cada pixel le asigno el pixel en su imagen menos el promedio sobre la raiz de la cantidad de imagenes menos uno**/
	return;
}

void generateMx()
{
	Xt = X;
	transpose(Xt);
	Mx = mult(Xt,X); /** Genero Mx matriz de covarianza como X^t por X **/
}

double norm(vector<double> &vec) /** Calculo norma 2 de vec **/
{
	double res = 0;
	for(int i=0;i<(int)vec.size();i++)
		res += vec[i]*vec[i];
	return sqrt(res);
}

vector<vector<double> > Id(int n) /** Identidad de n x n **/
{
    vector<vector<double> > A(n,vector<double>(n,0));
    for(int i = 0; i< n; i++)
        A[i][i]=1;
    return A;
}

vector<double> u,v;
vector<vector<double> > Q,R,aux,aux2;

void householder(vector<vector<double> > &A)
{
    /** Factorizacion QR de Householder en o(n^3) **/
	int n = A.size();
	aux.clear();
	aux2.clear();
	aux.resize(n,vector<double>(n));
	aux2.resize(n,vector<double>(n));
	R = A;
	Q = Id(n);
	for(int i=0;i<n-1;i++)
	{
	    u.clear();
	    for(int j=0;j<i;j++)
            u.push_back(0);
	    for(int j=i;j<n;j++)
            u.push_back(R[j][i]);
        double alpha = norm(u);
        if(abs(abs(alpha)-abs(u[i])) < 1e-6) /** Si debajo de la diagonal son todos ceros no itero **/
            continue;
        if(alpha*u[i]>0) /** Cambio el signo si es necesario **/
            alpha *= -1.;
        u[i] += alpha;
        alpha = norm(u);
        v.clear();
        for(int j=0;j<n;j++)
            v.push_back(u[j]/alpha);
        /** Inicio calculo R **/
        for(int j=0;j<n;j++)
            u[j] = 0;
        for(int t=0;t<n;t++)
        for(int j=0;j<n;j++)
            u[j] += v[t]*R[t][j];
        for(int j=0;j<n;j++)
        for(int t=0;t<n;t++)
            aux[j][t] = 2*v[j]*u[t];
        for(int j=0;j<n;j++)
        for(int t=0;t<n;t++)
            R[j][t] -= aux[j][t];
        /** Fin calculo R **/
        /** Inicio calculo Q **/
        for(int j=0;j<n;j++)
            u[j] = 0;
        for(int t=0;t<n;t++)
        for(int j=0;j<n;j++)
            u[j] += Q[j][t]*v[t];
        for(int j=0;j<n;j++)
        for(int t=0;t<n;t++)
            aux2[j][t] = 2*u[j]*v[t];
        for(int j=0;j<n;j++)
        for(int t=0;t<n;t++)
            Q[j][t] -= aux2[j][t];
        /** Fin calculo Q **/
	}
	return;
}

const double delta = 15;

int iteraciones;

bool sigoIterando(vector<vector<double> > &A)
{
	double res = 0.;
	int n = A.size();
	for(int i=0;i<n;i++)
	    for(int j=0;j<i;j++)
		    res += fabs(A[i][j]);
    printf("Iteracion %d: %lf\n",++iteraciones, res);
	return res>delta;
	/** Itero hasta que los elementos debajo de la diagonal sumen menos de 15 **/
}

void eig(vector<vector<double> > &A, vector<vector<double> > &auVec)
{
	int n = A.size(); /// A es cuadrada
	auVec = Id(n);
	iteraciones = 0;
	while(sigoIterando(A)) /** Condicion de parada **/
	{
		householder(A); /** Calculo QR con Householder **/
		A = mult(R,Q); /** Multiplico RQ para obtener la nueva A que es la matriz de covarianza **/
		auVec = mult(auVec,Q); /** Multiplico todas las Q para obtener los autovectores **/
	}
	/** Ordeno los autovectore segun la magnitud de los autovalores **/
	vector<pair<int,int> > aux(n);
	for(int i=0;i<n;i++)
		aux[i] = make_pair(A[i][i],i);
	sort(aux.begin(),aux.end()); /** Ordeno los autovalores **/
	reverse(aux.begin(),aux.end());
	matrizAuxiliar = auVec;
	for(int i=0;i<n;i++)
	for(int j=0;j<n;j++)
	{
		auVec[i][j] = matrizAuxiliar[aux[i].second][j];
		/** Asigno a la i-esima columna el autovector correspondiente al i-esimo autovalor **/
	}
	/** Fin ordenamiento de autovectores **/
	return;
}

vector<vector<double> > tc;

vector<double> calctc(vector<double> imagen, int k) /** Calculo la transfomacion caracteristica de imagen con parametro k **/
{
	vector<double> res(k,0);
	int n = imagen.size();
	for(int i=0;i<k;i++)
	for(int j=0;j<n;j++)
	{
		res[i]+=av[i][j]*imagen[j];
	}
	return res;
}

void fillTC(int k) /** LLeno tc con las transformaciones caracteristicas **/
{
	int n = input.size();
	tc.resize(n,vector<double>(k,0));
	for(int i=0;i<n;i++)
	{
		tc[i] = calctc(input[i],k);
	}
	return;
}

double dist(vector<double> &v1, vector<double> &v2, int norm)
{
	double res = 0;
	if(norm==0) /** Norma infinito de v1 - v2 | Distancia infinito de v1 a v2 **/
	{
	    for(int i=0;i<(int)v1.size();i++)
            res = max(res,v1[i]-v2[i]);
	}
	if(norm==1) /** Norma 1 de v1 - v2 | Distancias 1 de v1 a v2 **/
	{
	    for(int i=0;i<(int)v1.size();i++)
            res += abs(v1[i]-v2[i]);
	}
    if(norm==2) /** Norma 2 de v1 - v2 | Distancia 2 de v1 a v2 **/
    {
        for(int i=0;i<(int)v1.size();i++)
            res += (v1[i]-v2[i])*(v1[i]-v2[i]);
        res = sqrt(res);
    }
	return res;
}

int parse(char* c)
{
    int res = 0;
    int pos = 0;
    while(c[pos]!=0)
    {
        res *= 10;
        res += c[pos]-'0';
        pos++;
    }
    return res;
}

void usage()
{
    cout << "Uso: ./OCR <k> <imp> <norma>" << endl;
    cout << "donde k = cantidad de columnas a tomar de las transformaciones "<< endl;
    cout << "      imp = 0 (usando nearest neighbours, 1 usando distancia al promedio " << endl;
    cout << "      norma = 0 (norma infinito), 1 (norma 1), 2 (norma 2) " << endl;
    return;
}

int main(int argc, char* argv[])
{
    if(argc!=4){ usage(); exit(1);}
    int k = parse(argv[1]), imp = parse(argv[2]), norm = parse(argv[3]);
    /** k es el parametro k del enunciado, imp es la implementacion y norm es la norma que usamos para medir distancias **/
	FILE* v = fopen("../datos/trainingImages.txt","r");
	int n, t;
    int g;
    const int training_count = 5000; /** Usamos 5000 imagenes de entrenamiento **/
    const int test_count = 1000; /** Usamos 1000 imagenes de test **/
	fscanf(v,"%d %d",&n,&t);
	input.clear();
	input.resize(training_count,vector<double>(t));
	testImages.resize(test_count,vector<double>(t));
	/** Comienzo lectura imagenes de entrenamiento y test **/
	for(int i=0;i<training_count;i++)
    {
        for(int j=0;j<t;j++)
        {
            g = fscanf(v,"%lf",&input[i][j]);
        }
    }
    for(int i=0;i<test_count;i++)
    {
		for(int j=0;j<t;j++)
		{
			g = fscanf(v,"%lf",&testImages[i][j]);
		}
	}
	fclose(v);
	/** Fin lectura imagenes de entrenamiento y test **/
	/** Comienzo lectura labels de entrenamiento y test **/
    v = fopen("../datos/trainingLabels.txt","r");
    fscanf(v,"%d",&n);
    labels.resize(training_count);
    testLabels.resize(test_count);
    for(int i=0;i<training_count;i++)
    {
        g = fscanf(v,"%d",&labels[i]);
    }
    for(int i=0;i<test_count;i++)
    {
		g = fscanf(v,"%d",&testLabels[i]);
	}
	fclose(v);
	/** Fin lectura labels de entrenamiento y test **/
	v = fopen("V.txt","r");
	if(v==NULL) /** Si V no existe la genero **/
	{
        generateX(); /** Genero la matriz X **/
        generateMx(); /** Genero Mx la matriz de covarianza **/
        eig(Mx,av); /** Calculo los autovectores de la matriz de covarianza **/
        fclose(v);
        v = fopen("V.txt","w"); /** Escribo la matriz V en un archivo **/
        fprintf(v,"%d %d\n",av.size(),av[0].size());
        for(int i=0;i<av.size();i++)
        {
            for(int j=0;j<av[0].size();j++)
                fprintf(v,"%.6lf ",av[i][j]);
            fprintf(v,"\n");
        }
	}
	else /** Si V ya fue generada previamente la levanto del archivo **/
	{
	    int N,M;
	    fscanf(v,"%d %d",&N,&M);
	    av.clear();
	    av.resize(N,vector<double>(M));
	    for(int i=0;i<N;i++)
	    for(int j=0;j<M;j++)
            fscanf(v,"%lf",&av[i][j]);

	}
    fillTC(k); /** Genero las transformaciones caracteristicas de cada imagen de entrenamiento **/
    vector<double> vec;
    vector<pair<double,int> > distancias;
    vector<int> cant(10);
    int bien = 0;
    int mal = 0;
    if(imp==0)
    {
        /** La implementacion 0 toma las 100 imagenes mas cercanas y toma la mas repetida de esas 100 **/
        for(int i=0;i<test_count;i++) /** Iteramos sobre la imagenes de test **/
        {
            vec = calctc(testImages[i],k); /** Calculamos la transformacion caracteristica de la imagen dada **/
            distancias.resize(training_count);
            for(int j=0;j<training_count;j++)
            {
                distancias[j] = make_pair(dist(tc[j],vec,norm),j);
                /** Guardamos en un par la distancia a cada imagen de entrenamiento con la imagen **/
            }
            sort(distancias.begin(),distancias.end()); /** Ordenamos por distancia **/
            for(int j=0;j<10;j++)
                cant[j] = 0;
            for(int j=0;j<100;j++)
            {
                cant[labels[distancias[j].second]]++;
                /** Contamos cuantas imagenes hay de las 100 mas cercanas con cada label **/
            }
            int cualEs = 0;
            for(int j=0;j<10;j++)
            if(cant[j]>cant[cualEs])
                cualEs = j;
            if(testLabels[i]==cualEs)
                bien++;
            else
                mal++;
        }
    }
    else if(imp == 1)
    {
        /** La implementacion 1 toma el promedio de las transformaciones de cada digito y compara distancia a cada promedio **/
        vector<double> promedio[10];
        int cuantos[10];
        memset(cuantos,0,sizeof(cuantos));
        /** Comienzo calculo promedios **/
        for(int i=0;i<10;i++)
        {
            promedio[i].clear();
            promedio[i].resize(k,0);
        }
        for(int i=0;i<training_count;i++)
        {
            cuantos[labels[i]]++;
            for(int j=0;j<k;j++)
                promedio[labels[i]][j] += tc[i][j];
        }
        for(int i=0;i<10;i++)
        for(int j=0;j<k;j++)
        if(cuantos[i]!=0)
            promedio[i][j] /= (double)cuantos[i];
        /** Fin calculo promedios **/
        for(int i=0;i<test_count;i++)
        {
            vec = calctc(testImages[i],k);
            distancias.resize(10);
            for(int j=0;j<10;j++)
            {
                distancias[j] = make_pair(dist(promedio[j],vec,norm),j);
                /** Comparamos distancia al promedio de cada digito **/
            }
            sort(distancias.begin(),distancias.end());
            int cualEs = distancias[0].second;
            /** Ordenamos y nos quedamos con el mas cercano **/
            if(testLabels[i]==cualEs)
                bien++;
            else
                mal++;
        }
    }
	cout << bien <<" "<<mal<<" " << (double)bien/(double)(bien+mal)<< endl;
	/** Imprimimos cantidad de hits, cantidad de misses y proporcion de hits **/
    return 0;
}

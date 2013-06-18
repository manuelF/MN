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
	X.clear();
	X.resize(n,vector<double>(m));
	for(int i=0;i<n;i++)
    	for(int j=0;j<m;j++)
	    	X[i][j] = (input[i][j]-average[j])/sqrt(n-1);
	return;
}

void generateMx()
{
	Xt = X;
	transpose(Xt);
	Mx = mult(Xt,X);
}

double norm(vector<double> &vec)
{
	double res = 0;
	for(int i=0;i<(int)vec.size();i++)
		res += vec[i]*vec[i];
	return sqrt(res);
}

vector<vector<double> > Id(int n)
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
        if(abs(abs(alpha)-abs(u[i])) < 1e-6)
            continue;
        if(alpha*u[i]>0)
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

const double delta = 1e-5;

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
}

bool sim(const vector<vector<double> > & A)
{
    int n = A.size();
    int m = A[0].size();
    assert(n==m);
    for(int i =0; i<n; i++)
    {
        for(int j = 0; j<i; j++)
        {
            assert(A[i][j]==A[j][i]); /// Ojo! son doubles, comparar con epsilon
        }
    }
    return true;
}

void eig(vector<vector<double> > &A, vector<vector<double> > &auVec)
{
	int n = A.size(); /// A es cuadrada
	auVec = Id(n);
	iteraciones = 0;
	while(sigoIterando(A))
	{
		householder(A);
		A = mult(R,Q);
		auVec = mult(auVec,Q);
	}
	vector<pair<int,int> > aux(n);
	for(int i=0;i<n;i++)
		aux[i] = make_pair(A[i][i],i);
	sort(aux.begin(),aux.end());
	reverse(aux.begin(),aux.end());
	matrizAuxiliar = auVec;
	for(int i=0;i<n;i++)
	for(int j=0;j<n;j++)
	{
		auVec[i][j] = matrizAuxiliar[aux[i].second][j];
	}
	return;
}

vector<vector<double> > tc;

vector<double> calctc(vector<double> imagen, int k)
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

void fillTC(int k)
{
	int n = input.size();
	tc.resize(n,vector<double>(k,0));
	for(int i=0;i<n;i++)
	{
		tc[i] = calctc(input[i],k);
	}
	return;
}

double dist(vector<double> &v1, vector<double> &v2)
{
	double res = 0;
	for(int i=0;i<(int)v1.size();i++)
		res += (v1[i]-v2[i])*(v1[i]-v2[i]);
	return sqrt(res);
}

int main()
{
	int k = 28*28;
	FILE* v = fopen("../datos/trainingImages.txt","r");
	int n, t;
    int g;
    const int training_count = 5000;
    const int test_count = 1000;
	fscanf(v,"%d %d",&n,&t);
	input.clear();
	input.resize(training_count,vector<double>(t));
	testImages.resize(test_count,vector<double>(t));
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
	v = fopen("V.txt","r");
	if(v==NULL)
	{
        generateX();
        generateMx();
        eig(Mx,av);
        fclose(v);
        v = fopen("V.txt","w");
        fprintf(v,"%d %d\n",av.size(),av[0].size());
        for(int i=0;i<av.size();i++)
        {
            for(int j=0;j<av[0].size();j++)
                fprintf(v,"%.6lf ",av[i][j]);
            fprintf(v,"\n");
        }
	}
	else
	{
	    int N,M;
	    fscanf(v,"%d %d",&N,&M);
	    av.clear();
	    av.resize(N,vector<double>(M));
	    for(int i=0;i<N;i++)
	    for(int j=0;j<M;j++)
            fscanf(v,"%lf",&av[i][j]);

	}
    fillTC(k);
    vector<double> vec;
    vector<pair<double,int> > distancias;
    vector<int> cant(10);
    int bien = 0;
    int mal = 0;
    for(int i=0;i<test_count;i++)
    {
		vec = calctc(testImages[i],k);
		distancias.resize(training_count);
		for(int j=0;j<training_count;j++)
		{
			distancias[j] = make_pair(dist(tc[j],vec),j);
		}
		sort(distancias.begin(),distancias.end());
		for(int j=0;j<10;j++)
			cant[j] = 0;
		for(int j=0;j<100;j++)
		{
			cant[labels[distancias[j].second]]++;
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
	cout << bien <<" "<<mal<<" " << (double)bien/(double)(bien+mal)<< endl;
    return 0;
}

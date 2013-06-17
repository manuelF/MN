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
	cerr << n <<" "<< m <<" "<<t<<endl;
	vector<vector<double> > res(n,vector<double>(m,0));
	for(int i=0;i<n;i++)
	{
	    cerr << i<< endl;
    	for(int j=0;j<m;j++)
	        for(int k=0;k<t;k++)
		        res[i][j] += A[i][k]*B2[j][k];
    }
    cerr <<"Termino mult"<< endl;
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
/*	int n = Mx.size();
	for(int i=0;i<n;i++)
	for(int j=0;j<n;j++)
		Mx[i][j] /= (n-1);*/
}

double norm(vector<double> &vec)
{
	double res = 0;
	for(int i=0;i<(int)vec.size();i++)
		res += vec[i]*vec[i];
	return sqrt(res);
}

vector<vector<double> > aux;

void householder(vector<vector<double> > &A, vector<vector<double> > &Q, vector<vector<double> > &R)
{
	int n = A.size();
	Q.clear();
	Q.resize(n,vector<double>(n,0));
	for(int i=0;i<n;i++)
		Q[i][i] = 1;
	aux.clear();
	aux.resize(n,vector<double>(n));
	R = A;
	for(int i=0;i<n;i++)
	{
		vector<double> v;
		for(int j=i;j<n;j++)
			v.push_back(R[j][i]);
		double alpha = norm(v);
		v[0] = alpha;
		for(int j=1;j<(int)v.size();j++)
			v[j] = 0;
		for(int j=0;j<(int)v.size();j++)
			v[j] = R[j+i][i] - v[j];
		alpha = norm(v);
		for(int j=0;j<(int)v.size();j++)
			v[j] = v[j]/alpha;
		for(int a=0;a<n;a++)
        {
            for(int b=0;b<n;b++)
            {
                if(a<i||b<i)
                {
                    if(a==b)
                        aux[a][b] = 1.;
                    else
                        aux[a][b] = 0.;
                }
                else if(a==b)
                    aux[a][b] = 1.-2.*v[a-i]*v[b-i];
                else
                    aux[a][b] = -2.*v[a-i]*v[b-i];
            }
        }
		Q = mult(Q,aux);
		transpose(aux);
		R = mult(aux,R);
	}
	return;
}

const double delta = 1e-5;

bool sigoIterando(vector<vector<double> > &A)
{
	double res = 0.;
	int n = A.size();
	for(int i=0;i<n;i++)
	    for(int j=0;j<i;j++)
		    res += abs(A[i][j]);
	return res>delta;
}

vector<vector<double> > Id(int n)
{
    vector<vector<double> > A(n,vector<double>(n,0));
    for(int i = 0; i< n; i++)
        A[i][i]=1;
    return A;
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
	vector<vector<double> > Q,R;
	int n = A.size(); /// A es cuadrada
	auVec = Id(n);
	int iter = 0;
	while(sigoIterando(A))
	{
		householder(A,Q,R);
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
	FILE* v = freopen("../datos/trainingImages.txt","r",stdin);
	int n, t;
    int g;
    const int training_count = 1000;
    const int test_count = 1000;
	cin >> n >> t;
	input.clear();
	input.resize(training_count,vector<double>(t));
	testImages.resize(test_count,vector<double>(t));
	for(int i=0;i<training_count;i++)
    {
        for(int j=0;j<t;j++)
        {
            g = scanf("%lf",&input[i][j]);
        }
    }
    for(int i=0;i<test_count;i++)
    {
		for(int j=0;j<t;j++)
		{
			g = scanf("%lf",&testImages[i][j]);
		}
	}
    v = freopen("../datos/trainingLabels.txt","r",stdin);
    cin >> n;
    labels.resize(training_count);
    testLabels.resize(test_count);
    for(int i=0;i<training_count;i++)
    {
        g = scanf("%d",&labels[i]);
    }
    for(int i=0;i<test_count;i++)
    {
		g = scanf("%d",&testLabels[i]);
	}
    generateX();
    cerr << 1 << endl;
    generateMx();
    cerr<<2 << endl;
    eig(Mx,av);
    cerr << 3 << endl;
    fillTC(k);
    cerr << 4 << endl;
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
			distancias[i] = make_pair(dist(tc[j],vec),j);
		}
		sort(distancias.begin(),distancias.end());
		for(int j=0;j<10;j++)
			cant[j] = 0;
		for(int j=0;j<100;j++)
			cant[labels[distancias[j].second]]++;
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

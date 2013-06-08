#include<iostream>
#include<algorithm>
#include<vector>
#include<cstring>
#include<cstdio>
#include<cmath>

using namespace std;

vector<vector<double> > input, X;
vector<int> labels;
vector<double> average;


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
		average[j] /= (double)m;
	X.clear();
	X.resize(n,vector<double>(m));
	for(int i=0;i<n;i++)
	for(int j=0;j<m;j++)
		X[i][j] = (input[i][j]-average[j])/sqrt(m-1);
	return;
}

double norm(vector<double> &vec)
{
	double res = 0;
	for(int i=0;i<(int)vec.size();i++)
		res += vec[i]*vec[i];
	return sqrt(res);
}

vector<vector<double> > mult(vector<vector<double> > &A, vector<vector<double> > &B)
{
	int n = A.size();
	int m = B[0].size();
	int t = A[0].size(); // tiene que ser igual a B.size();
	vector<vector<double> > res(n,vector<double>(m,0));
	for(int i=0;i<n;i++)
	for(int j=0;j<m;j++)
	for(int k=0;k<t;k++)
		res[i][j] += A[i][k]*B[k][j];
	return res;
}

void transpose(vector<vector<double> > &mat)
{
	for(int i=0;i<(int)mat.size();i++)
	for(int j=0;j<i;j++)
	{
		swap(mat[i][j],mat[j][i]);
	}
	return;
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
		Q = mult(Q,aux);
		transpose(aux);
		R = mult(aux,R);
	}
	return;
}

#define delta 1e-5

bool sigoIterando(vector<vector<double> > &A)
{
	double res = 0.;
	int n = A.size();
	for(int i=0;i<n;i++)
	for(int j=0;j<i;j++)
		res += abs(A[i][j]);
	return res>delta;
}

vector<vector<double> > autoVectores(vector<vector<double> > A) /// Ojo! A va por copia porque lo modifico.
{
	vector<vector<double> > Q,R;
	while(sigoIterando(A))
	{
		householder(A,Q,R);
		A = mult(R,Q);
	}
	return A;
}

int main()
{
	freopen("../datos/trainingImages.txt","r",stdin);
	int n, t;
	cin >> n >> t;
	input.clear();
	input.resize(n,vector<double>(t));
	for(int i=0;i<n;i++)
    for(int j=0;j<t;j++)
        scanf("%lf",&input[i][j]);
    freopen("../datos/trainingLabels.txt","r",stdin);
    cin >> n;
    labels.resize(n);
    cout << n << endl;
    for(int i=0;i<n;i++)
    {
        if(i%1000==0)
            cout << i<<endl;
        scanf("%d",&labels[i]);
    }
    generateX();
    return 0;
}

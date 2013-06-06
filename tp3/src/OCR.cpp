#include<iostream>
#include<algorithm>
#include<vector>
#include<cstring>
#include<cmath>

using namespace std;

vector<vector<double> > input, X;
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
		X[i][j] = (input[i][j]-average[j])/sqrt(n-1);
	return;
}

int main()
{
	
}

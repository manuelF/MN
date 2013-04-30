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

int main(int argc, char* argv[])
{
    return 0;
}

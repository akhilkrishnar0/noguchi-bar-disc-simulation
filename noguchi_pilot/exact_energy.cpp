#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;
struct P{double x[3],v[3],m,e,c;};
int main(int argc,char**argv){for(int z=1;z<argc;z++){ifstream f(argv[z],ios::binary|ios::ate);int n=f.tellg()/sizeof(P);f.seekg(0);vector<P>p(n);f.read((char*)p.data(),n*sizeof(P));double U=0,K=0;
#pragma omp parallel for reduction(+:U,K) schedule(dynamic,64)
for(int i=0;i<n;i++){for(int a=0;a<3;a++)K+=.5*p[i].m*p[i].v[a]*p[i].v[a];for(int j=0;j<i;j++){double rr=(p[i].e*p[i].e+p[j].e*p[j].e)/2;for(int a=0;a<3;a++)rr+=pow(p[i].x[a]-p[j].x[a],2);U-=p[i].m*p[j].m/sqrt(rr);}}cout<<setprecision(16)<<argv[z]<<','<<K<<','<<U<<','<<K+U<<endl;}}

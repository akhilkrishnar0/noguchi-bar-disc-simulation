#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;struct P{double x[3],v[3],m,e,c;};
vector<P>load(char*s){ifstream f(s,ios::binary|ios::ate);size_t n=f.tellg()/sizeof(P);f.seekg(0);vector<P>p(n);f.read((char*)p.data(),n*sizeof(P));return p;}
int main(int argc,char**argv){auto p=load(argv[1]),q=load(argv[2]);double K=0,U=0;
#pragma omp parallel for reduction(+:K,U) schedule(dynamic,32)
for(int i=0;i<(int)p.size();i++){for(int k=0;k<3;k++)K+=.5*p[i].m*p[i].v[k]*p[i].v[k];for(int j=0;j<(int)q.size();j++)if(j!=i){double rr=(p[i].e*p[i].e+q[j].e*q[j].e)/2;for(int k=0;k<3;k++)rr+=pow(p[i].x[k]-q[j].x[k],2);U-=p[i].m*q[j].m/sqrt(rr);}}
cout<<setprecision(16)<<K<<','<<U<<','<<K+U<<endl;}

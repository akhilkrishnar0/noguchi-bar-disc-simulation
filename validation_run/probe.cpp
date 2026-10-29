// Direct summation on azimuthal rings: no tree approximation in circular-force diagnostic.
#include <vector>
#include <fstream>
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std; struct P{double x[3],v[3],m,e,c;};
int main(int argc,char**argv){ifstream f(argv[1],ios::binary|ios::ate);int n=f.tellg()/sizeof(P);f.seekg(0);vector<P>p(n);f.read((char*)p.data(),n*sizeof(P));double com[3]={};for(auto&q:p)for(int k=0;k<3;k++)com[k]+=q.m*q.x[k];for(auto&q:p)for(int k=0;k<3;k++)q.x[k]-=com[k];const int nr=48,np=48;vector<double>a(nr),scatter(nr);
#pragma omp parallel for schedule(dynamic)
for(int i=0;i<nr;i++){double R=.025+i*.02,sum=0,sq=0;for(int l=0;l<np;l++){double ph=(l+.5)*2*M_PI/np,c=cos(ph),s=sin(ph),ar=0;for(auto&q:p){double dx=q.x[0]-R*c,dy=q.x[1]-R*s,dz=q.x[2];double r2=dx*dx+dy*dy+dz*dz+(.02*.02+q.e*q.e)/2;ar+=q.m*(dx*c+dy*s)/(r2*sqrt(r2));}sum+=ar;sq+=ar*ar;}a[i]=sum/np;scatter[i]=sqrt(max(0.,sq/np-a[i]*a[i]));}
cout<<setprecision(15)<<"R,mean_ar,azimuthal_std_ar\n";for(int i=0;i<nr;i++)cout<<.025+i*.02<<','<<a[i]<<','<<scatter[i]<<'\n';}

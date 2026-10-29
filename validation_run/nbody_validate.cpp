// 3D collisionless Barnes-Hut octree, monopoles per species, symmetric Plummer
// softening e_ij^2=(e_i^2+e_j^2)/2, kick-drift-kick leapfrog. G=1.
#include <array>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <numeric>
using namespace std;
struct P{double x[3],v[3],m,e,c;};
struct Node{double center[3],half;int start,end,child[8]; double mass[3],cm[3][3];};
vector<P> p,source;bool frozen=false;vector<Node> nodes;vector<int> order;vector<array<double,3>> acc;vector<double> phi;
double theta=.8;const double ep[3]={.02,.04,.004};
int build(int lo,int hi,array<double,3> c,double half,int depth){
 int id=nodes.size();Node n{};for(int k=0;k<3;k++)n.center[k]=c[k];n.half=half;n.start=lo;n.end=hi;fill(n.child,n.child+8,-1);
 for(int j=lo;j<hi;j++){auto &q=source[order[j]];int s=int(q.c);n.mass[s]+=q.m;for(int k=0;k<3;k++)n.cm[s][k]+=q.m*q.x[k];}
 for(int s=0;s<3;s++)if(n.mass[s])for(int k=0;k<3;k++)n.cm[s][k]/=n.mass[s];nodes.push_back(n);
 if(hi-lo<=8||depth>=40)return id;
 auto oct=[&](int i){return (source[i].x[0]>c[0])+2*(source[i].x[1]>c[1])+4*(source[i].x[2]>c[2]);};
 sort(order.begin()+lo,order.begin()+hi,[&](int a,int b){return oct(a)<oct(b);});
 int start=lo;for(int o=0;o<8;o++){int end=start;while(end<hi&&oct(order[end])==o)end++;if(end>start){auto cc=c;for(int k=0;k<3;k++)cc[k]+=((o>>k)&1?1:-1)*half/2;int ch=build(start,end,cc,half/2,depth+1);nodes[id].child[o]=ch;}start=end;}
 return id;
}
inline void add(int i,const double *x,double m,double e,double *a,double &pot){
 double d[3],r2=(p[i].e*p[i].e+e*e)/2;for(int k=0;k<3;k++){d[k]=x[k]-p[i].x[k];r2+=d[k]*d[k];}
 double inv=1/sqrt(r2),fac=m*inv*inv*inv;for(int k=0;k<3;k++)a[k]+=fac*d[k];pot-=m*inv;
}
void walk(int i,int id,double *a,double &pot){auto &n=nodes[id];
 bool leaf=true;for(int k=0;k<8;k++)if(n.child[k]>=0){leaf=false;break;}
 if(leaf){for(int j=n.start;j<n.end;j++){int q=order[j];if(q!=i)add(i,source[q].x,source[q].m,source[q].e,a,pot);}return;}
 bool contains=true;for(int k=0;k<3;k++)if(abs((frozen?source[i].x[k]:p[i].x[k])-n.center[k])>n.half)contains=false;
 bool accept=!contains;
 for(int s=0;s<3&&accept;s++)if(n.mass[s]){double d2=0;for(int k=0;k<3;k++)d2+=pow(p[i].x[k]-n.cm[s][k],2);if(4*n.half*n.half>=theta*theta*d2)accept=false;}
 if(accept){for(int s=0;s<3;s++)if(n.mass[s])add(i,n.cm[s],n.mass[s],ep[s],a,pot);}
 else for(int k=0;k<8;k++)if(n.child[k]>=0)walk(i,n.child[k],a,pot);
}
void forces(){if(!frozen||nodes.empty()){if(!frozen)source=p;nodes.clear();array<double,3> c{};double lo[3]={1e99,1e99,1e99},hi[3]={-1e99,-1e99,-1e99};for(auto&q:source)for(int k=0;k<3;k++){lo[k]=min(lo[k],q.x[k]);hi[k]=max(hi[k],q.x[k]);}double half=0;for(int k=0;k<3;k++){c[k]=(lo[k]+hi[k])/2;half=max(half,(hi[k]-lo[k])/2);}build(0,p.size(),c,half*1.00001,0);}
 #pragma omp parallel for schedule(dynamic,32)
 for(int i=0;i<(int)p.size();i++){double a[3]={},pot=0;walk(i,0,a,pot);for(int k=0;k<3;k++)acc[i][k]=a[k];phi[i]=pot;}
}
void snapshot(string dir,int step){ofstream f(dir+"/snapshot_"+to_string(step)+".bin",ios::binary);f.write((char*)p.data(),p.size()*sizeof(P));ofstream a(dir+"/acceleration_"+to_string(step)+".bin",ios::binary);a.write((char*)acc.data(),acc.size()*sizeof(array<double,3>));}
void diagnostic(ofstream &f,int step,double dt){double K=0,U=0,L[3]={},com[3]={},mom[3]={};for(int i=0;i<(int)p.size();i++){auto&q=p[i];for(int k=0;k<3;k++){K+=.5*q.m*q.v[k]*q.v[k];com[k]+=q.m*q.x[k];mom[k]+=q.m*q.v[k];L[k]+=q.m*(q.x[(k+1)%3]*q.v[(k+2)%3]-q.x[(k+2)%3]*q.v[(k+1)%3]);}U+=(frozen?1.:.5)*q.m*phi[i];}
 f<<step<<','<<step*dt<<','<<K<<','<<U<<','<<K+U;for(double z:L)f<<','<<z;for(double z:com)f<<','<<z;for(double z:mom)f<<','<<z;f<<'\n';f.flush();}
int main(int argc,char**argv){static_assert(sizeof(P)==72);if(argc!=7){cerr<<"Usage: nbody INPUT OUTPUT_DIR STEPS DT THETA FROZEN_0_OR_1\n";return 2;}string input=argv[1],dir=argv[2];int steps=stoi(argv[3]);double dt=stod(argv[4]);theta=stod(argv[5]);frozen=stoi(argv[6]);ifstream in(input,ios::binary|ios::ate);if(!in || in.tellg()<=0 || (size_t)in.tellg()%sizeof(P)){cerr<<"Invalid particle file\n";return 2;}size_t n=in.tellg()/sizeof(P);in.seekg(0);p.resize(n);in.read((char*)p.data(),n*sizeof(P));source=p;acc.resize(n);phi.resize(n);order.resize(n);iota(order.begin(),order.end(),0);nodes.reserve(n*2);
 auto t=chrono::steady_clock::now();forces();
 ofstream check(dir+"/force_validation.csv");check<<"particle,component,relative_error,absolute_error,direct_norm\n";
 for(int j=0;j<192;j++){int i=(j*313)%n;double a[3]={},pot=0;for(int k=0;k<(int)n;k++)if(k!=i)add(i,source[k].x,source[k].m,source[k].e,a,pot);double err=0,norm=0;for(int k=0;k<3;k++){err+=pow(a[k]-acc[i][k],2);norm+=a[k]*a[k];}check<<i<<','<<p[i].c<<','<<sqrt(err/norm)<<','<<sqrt(err)<<','<<sqrt(norm)<<'\n';}
 ofstream f(dir+"/diagnostics.csv");f<<setprecision(15)<<"step,t,K,U,E,Lx,Ly,Lz,cx,cy,cz,px,py,pz\n";diagnostic(f,0,dt);snapshot(dir,0);
 for(int st=1;st<=steps;st++){for(int i=0;i<(int)n;i++)for(int k=0;k<3;k++){p[i].v[k]+=.5*dt*acc[i][k];p[i].x[k]+=dt*p[i].v[k];}forces();for(int i=0;i<(int)n;i++)for(int k=0;k<3;k++)p[i].v[k]+=.5*dt*acc[i][k];
 if(st%10==0||st==steps){diagnostic(f,st,dt);cout<<st<<" t="<<st*dt<<" wall="<<chrono::duration<double>(chrono::steady_clock::now()-t).count()<<endl;}
 if(st%20==0||st==steps)snapshot(dir,st);
 }
 ofstream finalcheck(dir+"/force_final.csv");finalcheck<<"particle,component,relative_error\n";
 for(int j=0;j<192;j++){int i=(j*313)%n;double a[3]={},pot=0;for(int k=0;k<(int)n;k++)if(k!=i)add(i,source[k].x,source[k].m,source[k].e,a,pot);double err=0,norm=0;for(int k=0;k<3;k++){err+=pow(a[k]-acc[i][k],2);norm+=a[k]*a[k];}finalcheck<<i<<','<<p[i].c<<','<<sqrt(err/norm)<<'\n';}

}

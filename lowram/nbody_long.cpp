// Exploratory low-particle disc + live halo + live bulge experiment (G=1).
// Independent implementation inspired by Noguchi (2001), NOT a reproduction.
// Build: g++ -O3 -std=c++17 -fopenmp nbody_long.cpp -o nbody_long
#include <array>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <string>
#include <filesystem>
using namespace std;
struct P{double x[3],v[3],m,e,c;};
struct Node{double center[3],half;int start,end,ch[8];double mass[3],com[3][3];};
bool exclude_bulge=false;vector<P> p;vector<Node> tree;vector<int> order;vector<array<double,3>> a;vector<double> pot;double theta=.5;
constexpr double e_comp[3]={.05,.08,.04};
int construct(int lo,int hi,array<double,3> mid,double half,int depth){
 int id=tree.size();Node n{};for(int k=0;k<3;k++)n.center[k]=mid[k];n.half=half;n.start=lo;n.end=hi;fill(n.ch,n.ch+8,-1);
 for(int j=lo;j<hi;j++){const P &q=p[order[j]];int s=int(q.c);if(exclude_bulge&&s==2)continue;n.mass[s]+=q.m;for(int k=0;k<3;k++)n.com[s][k]+=q.m*q.x[k];}
 for(int s=0;s<3;s++)if(n.mass[s])for(int k=0;k<3;k++)n.com[s][k]/=n.mass[s];tree.push_back(n);
 if(hi-lo<=8||depth>=35||half<1e-10)return id;
 auto oct=[&](int i){return int(p[i].x[0]>mid[0])+2*int(p[i].x[1]>mid[1])+4*int(p[i].x[2]>mid[2]);};
 sort(order.begin()+lo,order.begin()+hi,[&](int i,int j){return oct(i)<oct(j);});
 int start=lo;for(int o=0;o<8;o++){int stop=start;while(stop<hi&&oct(order[stop])==o)stop++;if(stop>start){auto mm=mid;for(int k=0;k<3;k++)mm[k]+=((o>>k)&1?1:-1)*half/2;int ch=construct(start,stop,mm,half/2,depth+1);tree[id].ch[o]=ch;}start=stop;}
 return id;
}
inline void force_one(int i,const double *x,double m,double es,double *out,double &ph){
 double rr=(p[i].e*p[i].e+es*es)/2,d[3];for(int k=0;k<3;k++){d[k]=x[k]-p[i].x[k];rr+=d[k]*d[k];}
 double inv=1/sqrt(rr);double fac=m*inv*inv*inv;for(int k=0;k<3;k++)out[k]+=fac*d[k];ph-=m*inv;
}
void walk(int i,int id,double *out,double &ph){const Node &n=tree[id];bool leaf=true;for(int j=0;j<8;j++)if(n.ch[j]>=0)leaf=false;
 if(leaf){for(int j=n.start;j<n.end;j++){int k=order[j];if(k!=i&&!(exclude_bulge&&p[k].c==2))force_one(i,p[k].x,p[k].m,p[k].e,out,ph);}return;}
 bool contains=true;for(int k=0;k<3;k++)if(abs(p[i].x[k]-n.center[k])>n.half)contains=false;
 bool accept=!contains;for(int s=0;s<3&&accept;s++)if(n.mass[s]){double d2=0;for(int k=0;k<3;k++){double d=p[i].x[k]-n.com[s][k];d2+=d*d;}if(4*n.half*n.half>=theta*theta*d2)accept=false;}
 if(accept){for(int s=0;s<3;s++)if(n.mass[s])force_one(i,n.com[s],n.mass[s],e_comp[s],out,ph);}else for(int j=0;j<8;j++)if(n.ch[j]>=0)walk(i,n.ch[j],out,ph);
}
// Controlled external rotating quadrupole. It injects energy and angular momentum.
void bar(int i,double t,double *acc){const double strength=.024, scale=.22, pattern_speed=1.8, ramp=4.;
 double f=(t<=0?0:t>=ramp?1:0.5*(1-cos(M_PI*t/ramp)));
 double c=cos(pattern_speed*t),s=sin(pattern_speed*t),x=c*p[i].x[0]+s*p[i].x[1],y=-s*p[i].x[0]+c*p[i].x[1];
 double ss=x*x+y*y+scale*scale,q=x*x-y*y,den=pow(ss,-2.5);
 double gx=-strength*f*(2*x*ss-3*x*q)*den,gy=-strength*f*(-2*y*ss-3*y*q)*den;
 acc[0]-=c*gx-s*gy;acc[1]-=s*gx+c*gy;
}
void forces(bool forced,double t){tree.clear();double lo[3]={1e99,1e99,1e99},hi[3]={-1e99,-1e99,-1e99};array<double,3> mid{};
 for(const auto &q:p)for(int k=0;k<3;k++){lo[k]=min(lo[k],q.x[k]);hi[k]=max(hi[k],q.x[k]);}
 double half=0;for(int k=0;k<3;k++){mid[k]=(lo[k]+hi[k])/2;half=max(half,(hi[k]-lo[k])/2);}construct(0,int(p.size()),mid,half*1.000001+1e-9,0);
 #pragma omp parallel for schedule(dynamic,16)
 for(int i=0;i<(int)p.size();i++){double aa[3]={},pp=0;walk(i,0,aa,pp);if(forced)bar(i,t,aa);for(int k=0;k<3;k++)a[i][k]=aa[k];pot[i]=pp;}
}
void dump(const string &name){ofstream f(name,ios::binary);f.write((char*)p.data(),p.size()*sizeof(P));}
void diagnostics(ofstream &f,double t,bool forced){double K=0,U=0,L=0,sum=0,zr=0,rr=0,c2=0,s2=0;int nd=0;
 for(int i=0;i<(int)p.size();i++){const auto &q=p[i];K+=.5*q.m*(q.v[0]*q.v[0]+q.v[1]*q.v[1]+q.v[2]*q.v[2]);U+=.5*q.m*pot[i];L+=q.m*(q.x[0]*q.v[1]-q.x[1]*q.v[0]);
 if(q.c==0){double r=hypot(q.x[0],q.x[1]);zr+=q.x[2]*q.x[2];rr+=r;nd++;if(r>.05&&r<.6){c2+=(q.x[0]*q.x[0]-q.x[1]*q.x[1])/(r*r);s2+=2*q.x[0]*q.x[1]/(r*r);sum++;}}}
 f<<setprecision(14)<<t<<','<<K+U<<','<<L<<','<<sqrt(zr/nd)<<','<<rr/nd<<','<<hypot(c2,s2)/max(1.,sum)<<','<<forced<<'\n';f.flush();
}
void evolve(int steps,double dt,const string &mode,const string &dir,int cadence){
 filesystem::create_directories(dir);bool forced=mode=="forced";forces(forced,0);
 ofstream diag(dir+"/diagnostics.csv");diag<<"t,E_self,Lz,zrms,meanR,A2,forced\n";
 diagnostics(diag,0,forced);dump(dir+"/snap_0000.bin");
 for(int step=1;step<=steps;step++){
  double t=(step-1)*dt;
  for(int i=0;i<(int)p.size();i++){
   bool freeze=(mode=="halo"?p[i].c!=1:mode=="spheroids"?p[i].c==0:false);
   if(freeze)continue;
   for(int k=0;k<3;k++){p[i].v[k]+=.5*dt*a[i][k];p[i].x[k]+=dt*p[i].v[k];}
  }
  forces(forced,t+dt);
  for(int i=0;i<(int)p.size();i++){
   bool freeze=(mode=="halo"?p[i].c!=1:mode=="spheroids"?p[i].c==0:false);
   if(freeze)continue;for(int k=0;k<3;k++)p[i].v[k]+=.5*dt*a[i][k];
  }
  if(step%cadence==0||step==steps){diagnostics(diag,step*dt,forced);dump(dir+"/snap_"+to_string(step/cadence)+".bin");}
 }
}
void vertical_calibrate(){ // local vertical tensor-virial matching after spheroid relaxation
 forces(false,0);const int nb=12;double numerator[nb]={},denom[nb]={};int counts[nb]={};
 for(int i=0;i<(int)p.size();i++)if(p[i].c==0){double r=hypot(p[i].x[0],p[i].x[1]);int b=min(nb-1,int(r/1.05*nb));if(b<0)b=0; numerator[b]+=-p[i].m*p[i].x[2]*a[i][2];denom[b]+=p[i].m*p[i].v[2]*p[i].v[2];counts[b]++;}
 double numall=accumulate(numerator,numerator+nb,0.0),denall=accumulate(denom,denom+nb,0.0);cerr<<"Vertical W="<<numall<<" K="<<denall<<" pre ratio="<<denall/numall<<"\n";
 for(int i=0;i<(int)p.size();i++)if(p[i].c==0){int b=clamp(int(hypot(p[i].x[0],p[i].x[1])/1.05*nb),0,nb-1);
  double ratio=(counts[b]>=30&&denom[b]>0&&numerator[b]>0)?numerator[b]/denom[b]:numall/denall;
  p[i].v[2]*=sqrt(clamp(ratio,.005,100.));
 }
 forces(false,0);double kz=0,wz=0;for(int i=0;i<(int)p.size();i++)if(p[i].c==0){kz+=p[i].m*p[i].v[2]*p[i].v[2];wz-=p[i].m*p[i].x[2]*a[i][2];}
 cerr<<"Disc vertical virial post-calibration: "<<kz/wz<<" (necessary, not sufficient for equilibrium)\n";
}
int main(int argc,char **argv){static_assert(sizeof(P)==72);
 if(argc!=9){cerr<<"Usage: nbody_long INPUT OUTPUT MODE [halo|spheroids|live|forced] STEPS DT THETA OUTPUT_EVERY_STEPS CALIBRATE_0_OR_1\n";return 2;}
 string input=argv[1],out=argv[2],mode=argv[3];int steps=stoi(argv[4]),cadence=stoi(argv[7]);double dt=stod(argv[5]);theta=stod(argv[6]);bool calib=stoi(argv[8]);
 if(dt<=0||steps<0||cadence<=0||theta<=0||theta>=2||!(mode=="halo"||mode=="spheroids"||mode=="live"||mode=="forced")){cerr<<"Invalid arguments\n";return 2;}
 ifstream in(input,ios::binary|ios::ate);if(!in||in.tellg()<=0||(size_t)in.tellg()%sizeof(P))return 3;
 int n=(size_t)in.tellg()/sizeof(P);in.seekg(0);p.resize(n);in.read((char*)p.data(),n*sizeof(P));order.resize(n);iota(order.begin(),order.end(),0);a.resize(n);pot.resize(n);tree.reserve(n*2);
 exclude_bulge=(mode=="halo");if(calib)vertical_calibrate();evolve(steps,dt,mode,out,cadence);dump(out+"/final.bin");
 cout<<"completed "<<mode<<" n="<<n<<" t="<<steps*dt<<" outputs="<<steps/cadence+1<<"\n";
}
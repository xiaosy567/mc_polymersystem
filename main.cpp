
#include "cmc.h"
#include <ctime>

/*#define SQUARE(X) X*X
template<class T> inline T Energy_LJ1(T EPSILON, T SIGMA, T d_2) {//pay attention d_2=d*d; 
  return 4.0*(fabs(EPSILON)*HEXIAL(SQUARE(SIGMA)/d_2)-EPSILON*TRIPLE(SQUARE(SIGMA)/d_2)); 
}*/
int main(void) {
  /*cout<<"the angle is: "<<dihedral(0,0,1,
  	             0,0,0,
  	             1,0,0,
  	             1,-1,0)/_PI<<"*PI"<<endl;
  cout<<"the angle is: "<<atan2(1,2)<<endl;
  string shit=string(" i am a boy  _");
  cout<<"*"<<shit<<"*"<<endl;
  cout<<"*"<<Split_By_Underline_Header(BFilter(shit))<<"*"<<endl;
  cout<<"*"<<Split_By_Underline(BFilter(shit))<<"*"<<endl;*/
  /*double* p1; p1=new double[3];
  p1[0]=1;p1[1]=1;p1[2]=0;
  double* p2; p2=new double[3];
  p2[0]=0;p2[1]=1;p2[2]=1;
  double* p3; p3=new double[3];
  p3=p_X_p(p1,p2);
  cout<<p3[0]<<" "<<p3[1]<<" "<<p3[2]<<endl;*/
  //cout<<"*"<<Split(BFilter(shit))<<"*"<<endl;
  //CMolecule tmol;
  /*double x=2.5;
  double x2=x*x;
  double x6=x*x*x*x*x*x;
  double x12=x*x*x*x*x*x*x*x*x*x*x*x;
  cout<<Energy_LJ(1.0,1.0,1.0,x12,x6)<<endl;
  cout<<Energy_LJ1(1.0,1.0,x2)<<endl;*/
  for(int i=128; i<256; i*=2) {
    cmc mymcsystem;
    //mymcsystem;
    //mymcsystem.load("test.pdb",1.0,false);
    //tmol.writepdbinfo("result.pdb",false);
    mymcsystem.init_conformation_a(i,1,1.0,"test.pdb");
    mymcsystem.initialization();
    double start=clock();
    mymcsystem.run_with_stepnumber(i*100000);
    mymcsystem.rand_update(0);
    mymcsystem.run();
    double end=clock();
    cout<<" [time]="<<(end-start)/CLOCKS_PER_SEC<<"seconds. "<<endl;//16~18s
    mymcsystem.fout_conformation("result.pdb",false);
    mymcsystem.memo_free();
  }
  return 0;
}

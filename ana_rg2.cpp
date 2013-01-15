

#include "ssbf.h"
#include "aglcalc.h"
#include <cstdlib>
#include <algorithm>
#include <fstream>
using namespace std;

int main(void) {
  string tempstr;
  vector<string> tempvec;
  vector<double> data_x; 
  double data_x_sum=0.0;
  double data_x2_sum=0.0;
  int temp_sz;
  double tempnum;
  string filename;
  system("rm -fr rh2_file.lst");
  system("if [ `ls rh2*.dat | wc -l` -gt 0 ]; then ls rh2*.dat > rh2_file.lst; fi");
  ifstream rh2filelist("rh2_file.lst");
  if(rh2filelist==NULL) {
    cout<<" Error: Can not open file: [ rh2_file.lst ]"<<endl;
    exit(IOERROR);
  }
  cout<<endl<<" Loading filename from file [ rh2_file.lst ] ..."<<endl;
  ofstream rh2datalist("rh2_data.lst");
  if(rh2datalist==NULL) {
    cout<<" Error: Can not open file: [ rh2_data.lst ]"<<endl;
    exit(IOERROR);
  }
  cout<<endl<<" writing filename into file [ rh2_data.lst ] ..."<<endl;
  while( getline(rh2filelist, tempstr) ) {
    tempvec=Split(BFilter(tempstr)); // a file
    if(tempvec.size()!=0) {
      filename=tempvec[0];
      ifstream rh2file(filename.c_str());
      if(rh2file==NULL) {
        cout<<" Error: Can not open file: [ "<<filename<<" ]"<<endl;
        exit(IOERROR);
      }
      cout<<" loading data from:[ "<<filename<<" ] ..."<<endl;
      data_x_sum=0.0;
      data_x2_sum=0.0;
      data_x.clear();
      while( getline(rh2file, tempstr) ) {
        tempvec=Split(BFilter(tempstr));
        if(tempvec.size()!=0) {
          tempnum=atof(tempvec[0].c_str());
          data_x2_sum+=tempnum;
          tempnum=sqrt(tempnum);
          data_x_sum+=tempnum;
          data_x.push_back(tempnum);
        }
      }//all data loaded;
      rh2file.close();
      temp_sz=data_x.size();
      cout<<" "<<temp_sz<<" data loaded... "<<endl;
      rh2datalist<<atoi(filename.substr(9,5).c_str())<<" "
                 <<data_x_sum/temp_sz<<" "
                 <<*min_element(data_x.begin(),data_x.end())<<" "
                 <<*max_element(data_x.begin(),data_x.end())<<" "
                 <<data_x2_sum/temp_sz-(data_x_sum/temp_sz)*(data_x_sum/temp_sz)<<endl;
      
    }
  } 
  rh2datalist.close();
  rh2filelist.close(); 
  system("echo \"fn=\'rh2_data\'\" > tempfn.gpl"); 
  system("echo \"yname=\'log(R_{H}) [0.1nm]\'\" >> tempfn.gpl");
  system("gnuplot < rg2list.gpl"); 
  /////////////////////////////////////
  system("rm -fr rg2_file.lst");
  system("if [ `ls rg2*.dat | wc -l` -gt 0 ]; then ls rg2*.dat > rg2_file.lst; fi");
  ifstream rg2filelist("rg2_file.lst");
  if(rg2filelist==NULL) {
    cout<<" Error: Can not open file: [ rg2_file.lst ]"<<endl;
    exit(IOERROR);
  }
  cout<<endl<<" Loading filename from file [ rg2_file.lst ] ..."<<endl;
  ofstream rg2datalist("rg2_data.lst");
  if(rg2datalist==NULL) {
    cout<<" Error: Can not open file: [ rg2_data.lst ]"<<endl;
    exit(IOERROR);
  }
  cout<<endl<<" writing filename into file [ rg2_data.lst ] ..."<<endl;
  while( getline(rg2filelist, tempstr) ) {
    tempvec=Split(BFilter(tempstr)); // a file
    if(tempvec.size()!=0) {
      filename=tempvec[0];
      ifstream rg2file(filename.c_str());
      if(rg2file==NULL) {
        cout<<" Error: Can not open file: [ "<<filename<<" ]"<<endl;
        exit(IOERROR);
      }
      cout<<" loading data from:[ "<<filename<<" ] ..."<<endl;
      data_x_sum=0.0;
      data_x2_sum=0.0;
      data_x.clear();
      while( getline(rg2file, tempstr) ) {
        tempvec=Split(BFilter(tempstr));
        if(tempvec.size()!=0) {
          tempnum=atof(tempvec[0].c_str())+atof(tempvec[1].c_str())+atof(tempvec[2].c_str());
          data_x2_sum+=tempnum;
          tempnum=sqrt(tempnum);
          data_x_sum+=tempnum;
          data_x.push_back(tempnum);
        }
      }//all data loaded;
      rg2file.close();
      temp_sz=data_x.size();
      cout<<" "<<temp_sz<<" data loaded... "<<endl;
      rg2datalist<<atoi(filename.substr(9,5).c_str())<<" "
                 <<data_x_sum/temp_sz<<" "
                 <<*min_element(data_x.begin(),data_x.end())<<" "
                 <<*max_element(data_x.begin(),data_x.end())<<" "
                 <<data_x2_sum/temp_sz-(data_x_sum/temp_sz)*(data_x_sum/temp_sz)<<endl;
      //mean, low, high, covariance;
    }
  } 
  rg2datalist.close();
  rg2filelist.close();
  system("echo \"fn=\'rg2_data\'\" > tempfn.gpl");   
  system("echo \"yname=\'log(R_{G}) [0.1nm]\'\" >> tempfn.gpl");  
  system("gnuplot < rg2list.gpl"); 
  system("./creat_html.bash");
  return 0;
}

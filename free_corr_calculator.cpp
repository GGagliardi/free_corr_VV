#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <string>
#include <sys/stat.h>
#include <boost/math/quadrature/gauss_kronrod.hpp>

using namespace std;

using Gauss_Kronrod= boost::math::quadrature::gauss_kronrod<double, 15>;

const double tol=1e-12;


template <typename T>
string to_string_with_precision(T a_value, const int n)
{
    ostringstream out;
    out.precision(n);
    out << fixed << a_value;
    return out.str();
}


void crash(string Message) {
  cout << Message <<endl;
  exit(-1);
  return ;
}

void Print_To_File(const vector<string>& row_id, const vector<vector<double>> &data, string Path, string MODE, string Header) {

  ofstream Print;
  if(MODE=="APP") Print.open(Path, ofstream::app);
  else Print.open(Path, ofstream::out);

  Print.precision(10);

  if(data.size() == 0) crash("In Print_To_File an empty vectors has been provided");

  unsigned int rows= data[0].size();

  if((signed)row_id.size() != 0 && row_id.size() != rows) crash(" In Print_To_File data size and data_id size do not match"); 

  Print<<Header<<endl;
  for(unsigned int j=0; j < rows; j++) {

    Print<<j;
    if(row_id.size() != 0) Print<<setw(20)<<row_id[j];
    
    for(unsigned int i=0; i<data.size(); i++) {
      if(data[i].size() != rows) crash("In Print_To_File, number of measurements is non-constant over columns. Exiting...");
      Print<<setw(20)<<data[i][j];
    }
    Print<<endl;
  }
  Print.close();
  return;
}



void Compute_free_corr(double am, int Tmax) {

 
  double Nc=3; //three colors

  auto C_cont = [&Nc, &am](int t) -> double {

		 
		  auto f = [&am, &t, &Nc](double x) {  return (Nc*2.0/pow(M_PI,2))*exp(-2.0*t*sqrt( pow(x,2) + pow(am,2)))*pow(x,2)*( 1.0/3 + pow(am,2)/( 6.0*( pow(am,2) + pow(x,2))));};

		  return Gauss_Kronrod::integrate( f, 0, numeric_limits<double>::infinity(), 5,tol);
		  
		};

  auto ptm2 = [](double p1,double p2, double p3) { return pow(sin(p1),2)+ pow(sin(p2),2)+ pow(sin(p3),2);};
  auto phm2 = [](double p1,double p2, double p3) { return pow(2*sin(p1/2),2) + pow(2*sin(p2/2),2) + pow(2*sin(p3/2),2);};
  auto fa = [&phm2](double p1, double p2, double p3) { return 1.0 + 0.5*phm2(p1,p2,p3);};
  auto fb = [&phm2](double p1, double p2, double p3) { return phm2(p1,p2,p3) + 0.5*( pow(2*sin(p1/2),2)*(pow(2*sin(p2/2),2)+pow(2*sin(p3/2),2)) + pow(2*sin(p2/2),2)*pow(2*sin(p3/2),2));};
  auto D2 = [&am, &fa, &fb](double p1,double p2, double p3) { return ( fb(p1,p2,p3) + pow(am,2))*( 4*fa(p1,p2,p3) + fb(p1,p2,p3) + pow(am,2));};
  auto W = [&am, &phm2, &fa, &fb](double p1, double p2, double p3) { return 0.25*pow( phm2(p1,p2,p3) - (fb(p1,p2,p3) + pow(am,2))/fa(p1,p2,p3)  ,2);};
  auto shaEp2 = [&am, &fa, &fb](double p1,double p2, double p3) { return pow( (fb(p1,p2,p3) + pow(am,2))/(2*fa(p1,p2,p3)),2) + (fb(p1,p2,p3)+pow(am,2))/fa(p1,p2,p3);};
  auto shaEp =[&shaEp2](double p1,double p2, double p3) { return sqrt(shaEp2(p1,p2,p3));};
  auto Ep =[&shaEp](double p1,double p2, double p3) {return asinh(shaEp(p1,p2,p3));};

  auto corr = [&am, &Nc, &ptm2, &D2, &W, &shaEp2, &Ep](double p1,double p2,double p3, double t, int r) -> double {

		return (32.0*Nc/(pow(2.0*M_PI,3)))*exp(-2*Ep(p1,p2,p3)*t)*( shaEp2(p1,p2,p3) +(1.0/3.0)*ptm2(p1,p2,p3) +pow(am,2)+ r*W(p1,p2,p3))/D2(p1,p2,p3);
	      };

  vector<double> corr_pert_res_tm(2*Tmax,0.0);
  vector<double> corr_pert_res_OS(2*Tmax, 0.0);

  for(int i=0;i<2;i++) { //loop over r

    int r= 2*i-1; //set Wilson parameter 

    int time;
     
    for(int t=1;t<=Tmax;t++) { //loop over time

      cout<<"r: "<<r<<" t: "<<t<<endl;
      time =t;
    
      auto corrp1= [&corr, &time, &r](double p1) -> double {  //boost only performs 1d integrals. 
     
		     auto corrp2 = [&corr, &p1, &time, &r](double p2) {

				     auto corrp3 = [&corr, &p1, &p2, &time, &r](double p3) {
						     return corr(p1, p2, p3, time,r);
						   };
				     return Gauss_Kronrod::integrate( corrp3, 0, M_PI, 5,tol);
				   };
		     return Gauss_Kronrod::integrate( corrp2, 0, M_PI, 5,tol);
		   };
  
      if(r==1)  corr_pert_res_OS[t] = Gauss_Kronrod::integrate(corrp1, 0, M_PI, 5, tol);
      else if(r==-1) corr_pert_res_tm[t] = Gauss_Kronrod::integrate(corrp1, 0, M_PI, 5, tol);
      else crash("Wilson parameter r: "+to_string(r)+" not recognized");

    }

  }


  //compute correlator in the continuum limit
  vector<double> corr_pert_cont(2*Tmax,0.0);
  for(int t=1;t<=Tmax;t++) corr_pert_cont[t] = C_cont(t);


  //take difference between corr_pert_cont and corr_pert_res_tm(OS)
  vector<double> a2corr_pert_res_tm(2*Tmax,0.0);
  vector<double> a2corr_pert_res_OS(2*Tmax,0.0);
  
  for(int t=0;t<2*Tmax;t++) {
    a2corr_pert_res_tm[t] = corr_pert_cont[t] - corr_pert_res_tm[t];
    a2corr_pert_res_OS[t] = corr_pert_cont[t] - corr_pert_res_OS[t];
  }
 
  //Print the result

  //create directory
  string dir_name= "vkvk_corr";
  string path="vkvk_corr/"+to_string(Tmax)+"_m"+to_string_with_precision(am,3);
  mkdir(dir_name.c_str(),0777);
  mkdir(path.c_str(),0777);


  Print_To_File({}, {a2corr_pert_res_tm, corr_pert_res_tm, corr_pert_cont}, path+"/OPPOR", "" , "#t  a2Corr Corr continuum");
  Print_To_File({}, {a2corr_pert_res_OS, corr_pert_res_OS, corr_pert_cont}, path+"/SAMER", "" , "#t  a2Corr Corr continuum");
  

  

  


  return;
}


int main() {




  //Example
  double am= 0.480;
  double Tmax= 24;
  Compute_free_corr(am, Tmax); //computes the free theory correlator up to Tmax (t=0 is never computed and C(0) is always set to zero)



  return 1;
}

#include <string>
#include <chrono>

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#include "TROOT.h"
#include "TCanvas.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussLegendreIntegrator.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TF2.h"
#include "TRandom3.h"

using namespace std; 
using namespace std::chrono; 

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

int symbol_newtonaA(int a, int b);
long double QnkA(int n, int k, long double P);
long double FjEA(long double Emin, long double Emax_temp);
long double array_avrA(double array[], int n);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/* 	This macro evaluates number of signals caused
	by cascades.

	The meaning of each starting parameters of macro:
	- I_criterium - list of particles IDs,

        - E_max - upper energy limit for integration [log(TeV)] (default = 6),
        - R_maximum - maximum radius for integration [cm] (deafault = 100 000),

	- Aprime - surface area of detectors [cm^2] (default = 25),
        - eff - efficiency of the detector [%] (default = 95),
	- n - number of detectors (default = 4),
	- Tdet - detector's registration time [s] (default = 10^(-7)),
	- f_fake - frequency of fake signals from the device [1/s] (default = 0.1),

	- Tp - time of measurement [s] (default = 3600*24*7 = 1 week),
	- fi - maximum vertical angle of arriving particles [degrees] (default = 90).

        - e_max - number of staps in integral over energy E (default = 50), 
        - eps - number of staps in integral over distance r (default = 1000),
        - eps_phi - number of staps in integral over angle theta (default = 100),
        - n_rand - number of draws of number of particles in the shower (default = 50),

        - our_l - lenght of particular arrays,
*/

void CascadeSignalEstimationApprox(int I_criterium[], string params_dir = "/home/jerzy/CREDO/Analiza/ResultCalculations/Mu", string plots_dir = get_current_dir_name(), double E_max = 6, double R_maximum = pow(10, 5), int e_max = 100, double Aprime = 25, double eff = 0.95, int n = 4, double Tdet = pow(10, -7), double f_fake = 0.1, int Tp = 3600*24*7, double fi = 90, int eps = 1500, int eps_phi = 100, int n_rand = 10, int our_l = 200)

{
  //-----------------------------------------------------------------

  auto start = high_resolution_clock::now(); 				//Time of computations

  //-----------------------------------------------------------------

  /*	In this part, previously generated files
	are read and several parameters are extracted
	from them.
  */

  double meanN_part_params[3], sigmaN_part_params[3]; 
  double meanRprc_params[3]; 
  double meanRrho_params[3];
  double Ntheta_params[4];
  double RhoN_param0[our_l], RhoN_param1[our_l];  

  //-------------------------

  string mN_param_file_name = params_dir;
  mN_param_file_name.append("/meanN(E)_fit_params.txt");

  ifstream mN_param_fit_file(mN_param_file_name);
  string mN_line_fit;

  int ipar_mN = 0;

  if(mN_param_fit_file.is_open()) {

    while(getline(mN_param_fit_file, mN_line_fit)) {  			//Parameters of meanN_part function
      meanN_part_params[ipar_mN] = stod(mN_line_fit);
      //cout << meanN_part_params[ipar_mN] <<endl;
      ipar_mN++;
    }

  }

  else cout<< "ERROR: COULD NOT FIND THE meanN(E) FIT PARAMETERS FILE" <<endl;

  mN_param_fit_file.close();

  //-------------------------

  string sN_param_file_name = params_dir;
  sN_param_file_name.append("/sigmaN(E)_fit_params.txt");

  ifstream sN_param_fit_file(sN_param_file_name);
  string sN_line_fit;

  int ipar_sN = 0;

  if(sN_param_fit_file.is_open()) {

    while(getline(sN_param_fit_file, sN_line_fit)) {  			//Parameters of sigmaN_part function
      sigmaN_part_params[ipar_sN] = stod(sN_line_fit);
      //cout << sigmaN_part_params[ipar_sN] <<endl;
      ipar_sN++;
    }

  }

  else cout<< "ERROR: COULD NOT FIND THE sigmaN(E) FIT PARAMETERS FILE" <<endl;

  sN_param_fit_file.close();

  //-------------------------

  string Rrho_param_file_name = params_dir;
  Rrho_param_file_name.append("/R_rho(E)_fit_params.txt");

  ifstream Rrho_param_fit_file(Rrho_param_file_name);
  string Rrho_line_fit;

  int ipar_Rrho = 0;

  if(Rrho_param_fit_file.is_open()) {

    while(getline(Rrho_param_fit_file, Rrho_line_fit)) {  		//Parameters of meanRrho function
      meanRrho_params[ipar_Rrho] = stod(Rrho_line_fit);
      //cout << meanRrho_params[ipar_Rrho] <<endl;
      ipar_Rrho++;
    }
  }

  else cout<< "ERROR: COULD NOT FIND THE R_rho(E) FIT PARAMETERS FILE" <<endl;

  Rrho_param_fit_file.close();

  //-------------------------

  string Rprc_param_file_name = params_dir;
  Rprc_param_file_name.append("/R_prc(E)_fit_params.txt");

  ifstream Rprc_param_fit_file(Rprc_param_file_name);
  string Rprc_line_fit;

  int ipar_Rprc = 0;

  if(Rprc_param_fit_file.is_open()) {

    while(getline(Rprc_param_fit_file, Rprc_line_fit)) {  		//Parameters of meanRprc function
      meanRprc_params[ipar_Rprc] = stod(Rprc_line_fit);
      //cout << meanRprc_params[ipar_Rprc] <<endl;
      ipar_Rprc++;
    }

  }

  else cout<< "ERROR: COULD NOT FIND THE Rprc(E) FIT PARAMETERS FILE" <<endl;

  Rprc_param_fit_file.close();

  //-------------------------

  string Ntheta_param_file_name = params_dir;
  Ntheta_param_file_name.append("/Npart(theta)_fit_params.txt");

  ifstream Ntheta_param_fit_file(Ntheta_param_file_name);
  string Ntheta_line_fit;

  int ipar_Ntheta = 0;

  if(Ntheta_param_fit_file.is_open()) {

    while(getline(Ntheta_param_fit_file, Ntheta_line_fit)) {		//Parameters of Ntheta function
      Ntheta_params[ipar_Ntheta] = stod(Ntheta_line_fit);
      //cout << Ntheta_params[ipar_Ntheta] <<endl;
      ipar_Ntheta++;
    }

  }

  else cout<< "ERROR: COULD NOT FIND THE Npart(theta) FIT PARAMETERS FILE" <<endl;

  Ntheta_param_fit_file.close();

  //-------------------------

  string RhoN_param_file_name = params_dir;
  RhoN_param_file_name.append("/Rho(r,N)_fit_params.txt");

  ifstream RhoN_param_fit_file(RhoN_param_file_name);
  string RhoN_line_fit;
  
  int iparN = 0;

  for(int i = 0; i < our_l; i++) {
    RhoN_param0[i] = 0.0;
    RhoN_param1[i] = 0.0;
  }

  if(RhoN_param_fit_file.is_open()) {

    while(getline(RhoN_param_fit_file, RhoN_line_fit)) {
      									//Parameter [0] of Rho(N/N0) function
      RhoN_param0[iparN] = stod(RhoN_line_fit);
      //cout << RhoN_param0[iparN] << " : ";
      getline(RhoN_param_fit_file, RhoN_line_fit);			//Parameter [1] of Rho(N/N0) function
      RhoN_param1[iparN] = stod(RhoN_line_fit);
      //cout << RhoN_param1[iparN] <<endl;
      iparN++;
    }

  }

  else cout<< "ERROR: COULD NOT FIND THE Rho(N/N0) FIT PARAMETERS FILE" <<endl;

  RhoN_param_fit_file.close();

  double RhoN_p0 = array_avrA(RhoN_param0, our_l);
  double RhoN_p1 = array_avrA(RhoN_param1, our_l);

  //-----------------------------------------------------------------

  /*	This part of code calculates flux of
	cosmic rays background.
  */

  double phi = (fi/180.0)*M_PI;					//Angle into radians
  double omega = 2*M_PI*(1 - cos(phi));				//Spherical angle

  TF1 f_bg("f_tlo","pow(cos(x), 3)*sin(x)", 0, phi);		//Relationship between effective 

  ROOT::Math::WrappedTF1 wf1(f_bg);
  ROOT::Math::GaussLegendreIntegrator ig;
  ig.SetFunction(wf1);
  ig.SetNumberPoints(1000);
  ig.SetRelTolerance(0.0001);
  double integral1 = ig.Integral(0, phi);			//Integral of previous function

  double I[3] = {70, 32, 2};						//Primary fluxes for different particle types
  string I_names[3] = {"muons", "electrons & positons", "photons"};	//Names of different particle types

  int I_len = sizeof(I)/sizeof(I[0]);				//Lengths of previous arrays			
  int Ic_len = sizeof(I_criterium)/sizeof(I_criterium[0]);

  double Ibg = 0;						//Background flux

  cout<<endl;
  cout<< "INCLUDED PARTICLES TYPES" <<endl;
  cout<<"________________________________________________" <<endl;
  cout<<endl;

  for(int j = 0; j < I_len; j++) {
  bool if_I = false;

    for(int jc = 0; jc <= Ic_len; jc++) {			//Checking if this type is in the list
      //cout<< I_criterium[jc] <<endl;
      if(j == I_criterium[jc]) {
        if_I = true;
        break;
      }
    }

    if(if_I == true) {						//Adding each type contribution to the flux
      Ibg = Ibg + I[j]*2*M_PI*integral1;
      cout<< "Including: " << I_names[j] <<endl;
    }
  }

  cout<<endl;
  cout<< "BACKGROUND LEVEL" <<endl;
  cout<<"________________________________________________" <<endl;
  cout<<endl;
  cout<< "Background flux					I_bg = "<< Ibg << " [m^(-2)s^(-1)]"<<endl;
  cout<< "Spherical angle of arriving particles	 	Omega = " << omega << " [rad]" <<endl;

  //-----------------------------------------------------------------

  /*	Here informations about detectors
	are printed and probability of detection
	is calculated.
  */

  long double Pbg;

  Pbg = 1 - exp(- Tdet * eff * ((Aprime/10000) * Ibg + f_fake));

  cout<<endl;
  cout<<"PARAMETERS OF THE SYSTEM" <<endl;
  cout<<"________________________________________________" <<endl;
  cout<<endl;

  cout<<"Number of detectors	 			n = " << n <<endl;
  cout<<"Detector's surface	 			Adet = " << Aprime << " [cm^2]" <<endl;
  cout<<"Detector's registration time 			Tdet = " << Tdet << " [s]" <<endl;
  cout<<"Detector's efficiency 				e_det = " << 100*eff << " %" <<endl;
  cout<<"Fake signals frequency 				f_fake = " << f_fake << " [s^(-1)]" <<endl;
  cout<<"\nProbability of detection in Tdet	 	Pbg = " << Pbg <<endl;
  cout<<"\nTime of measurement	 			Tp = " << Tp << " [s]" <<endl;
  cout<<"						Tp = " << Tp/60 << " [min]" <<endl;
  cout<<"						Tp = " << (Tp/60)/60 << " [h]"  <<endl;
  cout<<"						Tp = " << ((Tp/60)/60)/24 << " [days]" <<endl;
  cout<<"						Tp = " << (((Tp/60)/60)/24)/7 << " [weeks]" <<endl;

  //-----------------------------------------------------------------

  /*	In this part some arrays are created
  */

  long double Nsum[n+1];					//An array of final estimations for all eneries
  for(int ns = 0; ns <= n; ns++) {
    Nsum[ns] = 0;
  }


  long double** NiK = 0;					//An arry of coincidance signals for each energy sub-range
  NiK = new long double*[e_max];

  for(int i = 0; i < e_max; i++) {
    NiK[i] = new long double[n+1];

    for(int j = 0; j < n+1; j++) {
      NiK[i][j] = 0;
    }

  }

  //------------------------------------------------------------- 

  /*	Here energy ranges array is set.
  */

  long double E_lim[e_max];
  long double E_lim_temp = 1;
  long double kn = pow(10.0, E_max/double(e_max));

  for(int en = 0; en < e_max; en++) {

    E_lim[en] = E_lim_temp;							//Equally distributed on log scale
    E_lim_temp = E_lim_temp * kn;

  } 

  //------------------------------------------------------------

  /*	Here the functions for maximal radius (prc)
	are declared.
  */

  TF1 * f_meanR_prc = new TF1("f_meanR_prc" , " [0]*pow(x, [1]) + [2] ");	//meanR_rho(E) function
  f_meanR_prc->SetParameters(meanRprc_params[0], meanRprc_params[1], meanRprc_params[2]);

  //-----------------------------------------------------------------

  /*	Here the functions for maximal radius (rho)
	are declared.
  */

  TF1 * f_meanR_rho = new TF1("f_meanR_rho" , " [0]*pow(x, [1]) + [2] ");	//meanR_prc(E) function
  f_meanR_rho->SetParameters(meanRrho_params[0], meanRrho_params[1], meanRrho_params[2]);

  //-----------------------------------------------------------------

  /*	Here the functions for number of particles
	and factor for density are declared.
  */

  TF1 * f_meanNE = new TF1("f_meanNE" , " [0]*pow(x, [1]) + [2] ");		//meanN_part(E) function
  f_meanNE->SetParameters(meanN_part_params[0], meanN_part_params[1], meanN_part_params[2]);

  TF1 * f_sigmaNE = new TF1("f_sigmaNE" , " [0]*pow(x, [1]) + [2] "); 		//sigmaN_part(E) function
  f_sigmaNE->SetParameters(sigmaN_part_params[0], sigmaN_part_params[1], sigmaN_part_params[2]);


  TF1 F_fluct_N("F_fluct_N"," [0]*x + [2] ");  					//Density factor of fluctuations (F_N)
  F_fluct_N.SetParameters(RhoN_p0, RhoN_p1);

  //-----------------------------------------------------------------

  /* 	Here the F(theta) function is declared
  */

  TF1 F_theta("F_N_theta"," [0]*pow(cos([1]*(x + [2])), [3]) ");  		//Density factor for angle dependency (F_theta)
  F_theta.SetParameters(Ntheta_params[0], Ntheta_params[1], Ntheta_params[2], Ntheta_params[3]);

  //-----------------------------------------------------------------

  long double Gamma_const = 0.906402;

  TF2 F_N_Approx("F_N_Approx"," [2] *([0]*x) * (pow((y/100.0), -0.75)) * (pow((1 + [1]*(y/100.0)), -2.5)) ");			//Particles density function
  TF1 Fg_Approx("Fg_Approx"," [3] * (2*[2]) * x * ([0]) * (pow((x/100.0), -0.75)) * (pow((1 + [1]*(x/100.0)), -2.5)) ");	//Function for integration

  F_N_Approx.SetParameters( (1.25/(2.0*M_PI*Gamma_const)) * (pow((1.0/320.0), 1.25)), (1.0/320.0), (1.0/10000.0));
  Fg_Approx.SetParameters( (1.25/(2.0*M_PI*Gamma_const)) * (pow((1.0/320.0), 1.25)), (1.0/320.0), M_PI, (1.0/10000.0));

  //-----------------------------------------------------------------	

  for(int e = 0; e < e_max; e++) {						//Start ENERGY loop

  long double E_e = E_lim[e];

  //-----------------------------------------------------------------

  long double Rmax = R_maximum;
  if(R_maximum == 0) Rmax = f_meanR_prc->Eval(E_e);
  if(R_maximum == 1) Rmax = f_meanR_rho->Eval(E_e);

  //-----------------------------------------------------------------

  cout<<" \n" <<endl;
  cout<<"----------ENERGY: "<< E_lim[e] << " TeV----------" <<endl;
  cout<<" " <<endl;

  //-----------------------------------------------------------------							

  /*	This part calculates the integral over
	particles density function to see if number
	of particles in the shower is of proper order.
  */

  ROOT::Math::WrappedTF1 wf2(Fg_Approx);
  ROOT::Math::GaussLegendreIntegrator ig2;
  ig2.SetFunction(wf2);
  ig2.SetNumberPoints(10000);
  ig2.SetRelTolerance(0.0001);
  double Npart_max = ig2.Integral(0, Rmax);

  TRandom3 gen;
  long double Npart = f_meanNE->Eval(E_e);					//Expected number of particles in the function

  //-----------------------------------------------------------------

  cout<<"PARAMETERS OF THE CASCADE" <<endl;
  cout<<" " <<endl;
  cout<<"Expected number of particles for considered cascade		N = " << Npart <<endl;
  cout<<"Number of particles for found density function			N_part = " << Npart * Npart_max <<endl;
  cout<<"Radius in which 95% particles included				R_prc = " << f_meanR_prc->Eval(E_e) <<endl;
  cout<<"Radius in which particles density is greater than background	R_rho = " << f_meanR_rho->Eval(E_e) <<endl;


  //-----------------------------------------------------------------

  /*	Here the main "integral" is computed.
	A 3D matrix of probability of coincidence
	is calculated here. 	
  */  

  int Rcascn[eps];
  long double Aring[eps];

  long double*** PcascNR = 0;						//Probability of detections array
  PcascNR = new long double**[eps];

  for(int i = 0; i < eps; i++) {
    PcascNR[i] = new long double*[eps_phi];
    for(int rand = 0; rand < eps_phi; rand++) {
      PcascNR[i][rand] = new long double[n_rand];
    }
  }

  long double**** Qcasc = 0;						//Probability of coincidance array
  Qcasc = new long double***[n+1];

  for(int i = 0; i < n+1; i++) {
    Qcasc[i] = new long double**[eps];
    for(int j = 0; j < eps; j++) {
      Qcasc[i][j] = new long double*[eps_phi];
      for(int rand = 0; rand < eps_phi; rand++) {
        Qcasc[i][j][rand] = new long double[n_rand];
      }
    }
  }

  long double omega_c[eps_phi];
  long double Adet_phi[eps_phi];
  long double angle, sum_omega;


  long double Qtemp = 1.0;
  long double Q_n_temp = 1.0;


  for(int rka = 0; rka < eps; rka++) {					//Start DISTANCE loop

    Rcascn[rka] = (rka+1)*Rmax/eps;
    if(rka == 0) Aring[rka] = M_PI*pow(Rcascn[rka], 2);
    if(rka != 0) Aring[rka] = M_PI*(pow(Rcascn[rka], 2) - pow(Rcascn[rka-1], 2));

    sum_omega = 0;
    for(int phi_c = 0; phi_c < eps_phi; phi_c++) {			//Start ANGLE loop

      angle = double((phi_c+1)/double(eps_phi))*phi;
      long double F_theta_czyn = F_theta(angle);

      Adet_phi[phi_c] = Aprime*cos(angle);

      if(phi_c == 0) omega_c[phi_c] = 2*M_PI*(1 - cos(angle));
      if(phi_c != 0) omega_c[phi_c] = 2*M_PI*(1 - cos(angle)) - sum_omega;

      sum_omega = sum_omega + omega_c[phi_c];


      for(int f = 0; f < n_rand; f++) {					//Start FLUCTUATIONS loop 

      long double N_part = gen.Gaus(f_meanNE->Eval(E_e), f_sigmaNE->Eval(E_e));
      if(N_part < 0.0) N_part = f_meanNE->Eval(E_e);

      long double rho_temp = 0;
      if(F_theta(angle) >= 0.0) rho_temp = F_N_Approx(N_part, Rcascn[rka]) * F_theta_czyn;

      PcascNR[rka][phi_c][f] = Adet_phi[phi_c] * eff * rho_temp;
      PcascNR[rka][phi_c][f] = 1 - exp( - PcascNR[rka][phi_c][f]);	//Poisson distribution


      for(int k_temp = n; k_temp > 0; k_temp--) {			//Start K loop
  
        //--------------------------------------------------		//No coincidance with background - old method

        Qcasc[k_temp][rka][phi_c][f] = QnkA(n, k_temp, PcascNR[rka][phi_c][f]) * QnkA(n - k_temp, 0, Pbg);

        //--------------------------------------------------
    
      }  								//End K loop

      }									//End FLUCTUATIONS loop

    }									//End Angle loop

  }									//End DISTANCE loop

  //-----------------------------------------------------------------

  /*	This part calculates frequencies of
	the showers for given energy range
	and spherical angle.
  */

  cout<<" \n" <<endl;
  cout<<"CASCADES FREQUENCIES " <<endl;
  cout<<endl;

  long double E0, Emax, omega, FE;

  omega = 2*M_PI*(1 - cos(phi));

  E0 = 0.7*E_lim[e];							//Minimal energy of cascades [TeV]
  Emax = 0.7*E_lim[e+1];						//Maximal energy of cascades [TeV]

  if(e == e_max - 1) Emax = pow(10, E_max);

  FE = FjEA(E0*pow(10, 3), Emax*pow(10, 3));


  cout<< "Energy range from 		Emin =  " << E0 << " [TeV] to " << " Emax = " << Emax << " [TeV]" <<endl;
  cout<< "Spherical angle 		Omega(" << phi*180/M_PI << " degrees) = " << omega <<endl;
  cout<< "Cascades frequencies 		F(E) =  " << omega*FE << " [s^{-1}m^{-2}]" <<endl;
  cout<< "				     = " << 60*omega*FE << " [min^{-1}m^{-2}]" <<endl;
  cout<< "				     = " << 3600*omega*FE << " [h^{-1}m^{-2}]" <<endl;
  cout<< "				     = " << 3600*24*omega*FE << " [d^{-1}m^{-2}]" <<endl;
  cout<< "				     = " << 3600*24*7*omega*FE << " [week^{-1}m^{-2}]" <<endl;
  cout<< "				     = " << 3600*24*365*omega*FE << " [yr^{-1}m^{-2}]" <<endl;
  cout<<endl;

  //-----------------------------------------------------------------

  /*	Computing expected values and
	putting them into an array.
  */

  long double Ncasc_phi = 0;
  long double Ncasc_rka = 0;
  long double Nkcasc[n], Nk[n];

  Nkcasc[0] = 0;

  for(int k = 1; k <= n; k++) {						//Start K loop

    Nkcasc[k] = 0;

    for(int rka = 0; rka < eps; rka++) {
 
      Ncasc_rka = 0;
     
      for(int phi_c = 0; phi_c < eps_phi; phi_c++) {			//Start ANGLE loop

        for(int f = 0; f < n_rand; f++) {				//Start FLUCTUATIONS lopp
      
          Ncasc_phi = Qcasc[k][rka][phi_c][f] * (Aring[rka]/10000.0) * (1.0/n_rand) * Tp * FE * omega_c[phi_c];
          Ncasc_rka = Ncasc_rka + Ncasc_phi;

        }								//End FLUCTUATIONS lopp

      }									//End ANGLE loop

      Nkcasc[k] = Nkcasc[k] + Ncasc_rka;
    }

    cout<< "Avereage number of events for k = " << k <<" 	<N_casc(k)> = " << Nkcasc[k] <<endl;

    NiK[e][k] = Nkcasc[k];
    Nsum[k] = Nsum[k] + NiK[e][k];

  }									//End K loop

  //-----------------------------------------------------------------

  }									//End ENERGY loop


  //-----------------------------------------------------------------
  
  /*	Printing expected number of coincidences
	for all together and for certain energy ranges.
  */

  cout<<"\n" <<endl;
  cout<<"EXPECTED NUMBER OF COINCIDANCE SIGNALS" <<endl;
  cout<<"________________________________________________" <<endl;

  cout<<"k	|	N(k)		|	f(k)" <<endl;
  cout<<"________________________________________________" <<endl;

  for(int nk = 1; nk < n+1; nk++) {

    cout<< nk << "	|	" << Nsum[nk];
    cout<< "	|	" << Nsum[nk]/Tp << " [s^(-1)]";
    cout<< "	|	" << 60*Nsum[nk]/Tp << " [min^(-1)]";
    cout<< "	|	" << 60*60*Nsum[nk]/Tp << " [h^(-1)]";
    cout<< "	|	" << 60*60*24*Nsum[nk]/Tp << " [day^(-1)]" <<endl;

  }

  cout<<"________________________________________________" <<endl;
  cout<<"\nEnergy range |				";
  for(int k = 1; k < n+1; k++) cout<<" |	N(E, k = " << k << ") ";
  cout<<endl;
  cout<<"________________________________________________" <<endl;

  for(int ne = 0; ne < e_max; ne++) {
    if(ne < e_max - 1) cout<< "Energy range ( " << E_lim[ne] << " - " << E_lim[ne + 1] << " ) TeV :		";
    if(ne == e_max - 1) cout<< "Energy range ( " << E_lim[ne] << " - " << " 10e" << E_max << " ) TeV :		";

    for(int nk = 1; nk < n+1; nk++) cout<< NiK[ne][nk] << " |		";
    cout<<endl;
    
  }

  //-----------------------------------------------------------------

  /*	Creating histograms of coincidance signals
	and energy dependence.
  */

  TCanvas A("D");
  A.SetLogx(); 

  TMultiGraph * nE_g_All = new TMultiGraph();

  for(int nk = 1; nk < n+1; nk++) {

    TGraph * nE_g = new TGraph();

    for(int ne = 0; ne < e_max - 1; ne++) { 

      nE_g->SetPoint(ne + 1, E_lim[ne], NiK[ne][nk]/Nsum[nk]);
    
    }

    nE_g_All->Add(nE_g);

    string title_nE_g_All = " #LT N_{casc}(k, E) #GT for k: 1 - ";
    title_nE_g_All.append(to_string(n));

    nE_g_All->SetTitle(title_nE_g_All.c_str());
    nE_g_All->GetXaxis()->SetTitle("E [TeV]");
    nE_g_All->GetYaxis()->SetTitle("#LT N_{casc} #GT [%]");

    nE_g_All->GetXaxis()->SetLimits(0.7, pow(10, E_max));

    nE_g_All->Draw("AP");
  
    cout<<endl;
    A.SaveAs("N_cascApprox_All.png");					//Saving into .png file

    string title_nE_g = " #LT N_{casc}(k, E) #GT for k = ";
    title_nE_g.append(to_string(nk));

    nE_g->SetTitle(title_nE_g.c_str());
    nE_g->GetXaxis()->SetTitle("E [TeV]");
    nE_g->GetYaxis()->SetTitle("#LT N_{casc} #GT [%]");

    nE_g->SetMarkerStyle(8);
    nE_g->SetMarkerSize(1);
    nE_g->SetMarkerColor(nk);

    string name_nE_g = plots_dir;
    name_nE_g.append("/N_cascApprox(");
    name_nE_g.append(to_string(nk));
    name_nE_g.append(",E).png");
    const char* nE_g_name = name_nE_g.c_str();

    nE_g->Draw("AP");
    A.SaveAs(nE_g_name);						//Saving into .png file

    nE_g_All->Add(nE_g);

  }

  string title_nE_g_All = " #LT N_{casc}(k, E) #GT for k: 1 - ";
  title_nE_g_All.append(to_string(n));

  nE_g_All->SetTitle(title_nE_g_All.c_str());
  nE_g_All->GetXaxis()->SetTitle("E [TeV]");
  nE_g_All->GetYaxis()->SetTitle("#LT N_{casc} #GT [%]");

  nE_g_All->GetXaxis()->SetLimits(0.7, pow(10, E_max));

  nE_g_All->Draw("AP");
  
  cout<<endl;
  A.SaveAs("N_cascApprox_All.png");					//Saving into .png file

  nE_g_All->Delete();

  //-----------------------------------------------------------------

  auto stop = high_resolution_clock::now(); 				//Time of computations
  auto duration = duration_cast<microseconds>(stop - start); 

  cout<<"\n" <<endl;
  cout<< "Time of the computations: " << duration.count() << " * 10^{-6} s = " <<  duration.count()/pow(10, 6) << " s = "  << duration.count()/(pow(10, 6)*60) << " min" <<endl;

  //-----------------------------------------------------------------

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~				//End MACRO

int symbol_newtonaA(int a, int b)
{
  return (TMath::Factorial(a)/(TMath::Factorial(b)*TMath::Factorial(a - b)));
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

long double QnkA(int n, int k, long double P)
{
   return (symbol_newtonaA(n, k) * pow(P, k) * pow((1 - P), (n - k)));
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

long double FjEA(long double Emin, long double Emax)
{

  long double jE1 = 2.107 * pow(10, 4);
  long double jE2 = 3.739 * pow(10, 7); 
  long double jE3 = 7.494 * pow(10, 5);
  long double jE4 = 4.952 * pow(10, 9); 

  double gam1 = 2.648;
  double gam2 = 3.138;
  double gam3 = 2.903;
  double gam4 = 3.374;

  long double E0 = pow(10, 0);
  long double E1 = pow(10, 6.7);
  long double E2 = pow(10, 7.25);
  long double E3 = pow(10, 8.1);
  long double E4 = pow(10, 12);

  long double integralE1, integralE2, integralE3, integralE4;
  integralE1 = 0;
  integralE2 = 0;
  integralE3 = 0;
  integralE4 = 0;


  TF1 j_E1("j_E1", " [0]*pow(x, - [1]) ", E0, E1);		//Spektrum functions of cascades 1

  j_E1.SetParameter(0, jE1);
  j_E1.SetParameter(1, gam1);


  if((Emin >= E0 && Emin <= E1) || (Emax >= E0 && Emax <= E1)) {
    if(Emin >= E0) E0 = Emin;
    if(Emax <= E1) E1 = Emax;

    ROOT::Math::WrappedTF1 wfE1(j_E1);
    ROOT::Math::GaussLegendreIntegrator igE1;
    igE1.SetFunction(wfE1);
    igE1.SetNumberPoints(1000);
    igE1.SetRelTolerance(0.0001);

    integralE1 = igE1.Integral(E0, E1);
  }



  TF1 j_E2("j_E2", " [0]*pow(x, - [1]) ", E1, E2);		//Spektrum functions of cascades 2

  j_E2.SetParameter(0, jE2);
  j_E2.SetParameter(1, gam2);

  if((Emin >= E1 && Emin <= E2) || (Emax >= E1 && Emax <= E2)) {
    if(Emin >= E1) E1 = Emin;
    if(Emax <= E2) E2 = Emax;

    ROOT::Math::WrappedTF1 wfE2(j_E2);
    ROOT::Math::GaussLegendreIntegrator igE2;
    igE2.SetFunction(wfE2);
    igE2.SetNumberPoints(1000);
    igE2.SetRelTolerance(0.0001);

    integralE2 = igE2.Integral(E1, E2);
  }


  TF1 j_E3("j_E3", " [0]*pow(x, - [1]) ", E2, E3);		//Spektrum functions of cascades 3

  j_E3.SetParameter(0, jE3);
  j_E3.SetParameter(1, gam3);

  if((Emin >= E2 && Emin <= E3) || (Emax >= E2 && Emax <= E3)) {
    if(Emin >= E2) E2 = Emin;
    if(Emax <= E3) E3 = Emax;

    ROOT::Math::WrappedTF1 wfE3(j_E3);
    ROOT::Math::GaussLegendreIntegrator igE3;
    igE3.SetFunction(wfE3);
    igE3.SetNumberPoints(1000);
    igE3.SetRelTolerance(0.0001);

    integralE3 = igE3.Integral(E2, E3);
  }



  TF1 j_E4("j_E4", " [0]*pow(x, - [1]) ", E3, E4);		//Spektrum functions of cascades 4

  j_E4.SetParameter(0, jE4);
  j_E4.SetParameter(1, gam4);

  if((Emin >= E3 && Emin <= E4) || (Emax >= E3 && Emax <= E4)) {
    if(Emin >= E3) E3 = Emin;
    if(Emax <= E4) E4 = Emax;

    ROOT::Math::WrappedTF1 wfE4(j_E4);
    ROOT::Math::GaussLegendreIntegrator igE4;
    igE4.SetFunction(wfE4);
    igE4.SetNumberPoints(1000);
    igE4.SetRelTolerance(0.0001);

    integralE4 = igE4.Integral(E3, E4);

  }

  return (integralE1 + integralE2 + integralE3 + integralE4);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

long double array_avrA(double array[], int n)
{

  int N = 0;
  long double sum = 0.0;

  for(int i = 0; i < n; i++) {
    if(array[i] == 0.0) break;
    sum = sum + array[i];
    N++;
  }

  return (sum/double(N));

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


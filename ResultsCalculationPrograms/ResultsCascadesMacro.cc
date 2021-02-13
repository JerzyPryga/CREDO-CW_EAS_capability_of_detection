#include <string>
#include <iostream>
#include <vector>
#include <experimental/filesystem>
#include <dirent.h>
#include <codecvt>
#include <chrono>
#include <limits.h>
#include <unistd.h>
#include <bits/stdc++.h> 
#include <iostream> 
#include <sys/stat.h> 
#include <sys/types.h> 

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

#include "CascadeSignalEstimation.cc"						//inclding CREDO macro

using namespace std::chrono;
using namespace std;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*	Macro with all inputs included which
	executes previously written macros.
*/

void ResultsCascadesMacro(string current_dir) {					//Start MACRO

  //-------------------------------------------------------------

  auto start = high_resolution_clock::now(); 					//Starting counting time of computations

  //-------------------------------------------------------------

  /*	Here the starting parameters for the program
	are inputed. 

	The meaning of each parameter:

	- I_criterium - list of particles IDs,

	- params_dir - directory with parameters .txt files,
	- plots_dir - directory for storing plots and results,

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
	
	When user clicks ENTER the
	default value is set.
  */
  
  int Tp, n;
  int I_criterium[3] = {0, 1, 2};
  double Aprime, eff, Tdet, f_fake, fi;

  string params_dir = "/home/jerzy/CREDO/Analiza/ResultCalculations/Mu";
  string plots_dir = current_dir;

  double E_max, R_maximum;
  int e_max, eps, eps_phi, n_rand;

  int id_len = 3;

  string input;

  Tp = 7;
  n = 4;

  Aprime = 25.0;
  eff = 0.95;
  Tdet = 2*pow(10, -7);
  f_fake = 0.1;
  fi = 90.0;

  E_max = 6;
  R_maximum = 100000;

  e_max = 100;
  eps = 1000;
  eps_phi = 100;
  n_rand = 10;


  cout<<endl;
  cout<< "PARAMETERS OF THE SYSTEM" <<endl;

  cout<< "INPUT: number of detectors n (default = " << n << ") : ";				//Inputing n

  getline(cin, input);
  if(!input.empty()) {
    istringstream stream(input);
    stream >> n;
  }

  cout<< "INPUT: detectors area of the surface (each) Aprime [cm^2] (default = " << Aprime << ") : ";	//Inputing Aprime

  getline(cin, input);
  if(!input.empty()) {
    istringstream stream(input);
    stream >> Aprime;
  }

  cout<< "INPUT: detectors efficiency eff [%] (default = " << eff << ") : ";			//Inputing eff

  getline(cin, input);
  if(!input.empty()) {
    istringstream stream(input);
    stream >> eff;
  }

  cout<< "INPUT: time of coincidence Tdet [s] (default = " << Tdet << ") : ";			//Inputing Tdet

  getline(cin, input);
  if(!input.empty()) {
    istringstream stream(input);
    stream >> eff;
  }

  cout<< "INPUT: freqancy of fake signals [1/s] (default = " << f_fake << ") : ";		//Inputing f_fake

  getline(cin, input);
  if(!input.empty()) {
    istringstream stream(input);
    stream >> f_fake;
  }




  cout<<endl;
  cout<< "TIME OF MEASUREMENTS" <<endl;

  cout<< "INPUT: time of measurements Tp [days] (default = " << Tp << ") : ";//Inputing Tp

  getline(cin, input);
  if(!input.empty()) {
    istringstream stream(input);
    stream >> Tp;
  }

  Tp = Tp*60*60*24;


  cout<<endl;
  cout<< "CASCADE ASSUMPTIONS" <<endl;

  cout<< "INPUT: maximum vertical angle fi [degrees] (default = " << fi << ") : ";			//Inputing fi

  getline(cin, input);
  if(!input.empty()) {
    istringstream stream(input);
    stream >> fi;
  }


  cout<< "INPUT: maximum energy of considered showers [log(E/TeV)] (default = " << E_max << ") : ";		//Inputing E_max

  getline(cin, input);
  if(!input.empty()) {
    istringstream stream(input);
    stream >> E_max;
  }




  cout<<endl;
  cout<< "INTEGRATION OPTIONS" <<endl;

  cout<< "Integrate to R_prc (in which 95% of particles are included): " <<endl;			//Inputing If_Rprc
  cout<< "Press ENTER - YES" <<endl;
  cout<< "Type anything - NO" <<endl;

  bool If_Rprc;
  bool If_Rrho;

  getline(cin, input);
  if(!input.empty()) If_Rprc = false;
  else If_Rprc = true;

  input = "";

  if(If_Rprc == false) {

  cout<< "Integrate to R_rho (in which particles density is greater than background): " <<endl;	//Inputing If_Rrho
  cout<< "Press ENTER - YES" <<endl;
  cout<< "Type anything - NO" <<endl;

  getline(cin, input);
  if(!input.empty()) If_Rrho = false;
  else If_Rrho = true;

  input = "";

  }

  if((If_Rprc == false) && (If_Rrho == false)) {

  cout<< "INPUT: maximum distance to which it will be integrated R_maximum [cm] (default = " << R_maximum << ") : ";		//Inputing R_maximum

  getline(cin, input);
  if(!input.empty()) {
    istringstream stream(input);
    stream >> R_maximum;
  }

  }

  if(If_Rprc == true) R_maximum = 0;
  if(If_Rrho == true) R_maximum = 1;

  cout<< "INPUT: number of steps in integral over energy e_max (default = " << e_max << ") : ";		//Inputing e_max

  getline(cin, input);
  if(!input.empty()) {
    istringstream stream(input);
    stream >> e_max;
  }

  cout<< "INPUT: number of steps in integral over distance eps (default = " << eps << ") : ";		//Inputing eps

  getline(cin, input);
  if(!input.empty()) {
    istringstream stream(input);
    stream >> eps;
  }

  cout<< "INPUT: number of steps in integral over spherical angle eps_phi (default = " << eps_phi << ") : ";		//Inputing eps_phi

  getline(cin, input);
  if(!input.empty()) {
    istringstream stream(input);
    stream >> eps_phi;
  }

  cout<< "INPUT: number of draws for fluctuations in each step n_rand (default = " << n_rand << ") : ";		//Inputing n_rand

  getline(cin, input);
  if(!input.empty()) {
    istringstream stream(input);
    stream >> n_rand;
  }




  cout<<endl;
  cout<< "CRITERIA OF PARTICLE CHOOSING" <<endl;

  cout<< "\nID | Particle" <<endl;
  cout<< "0 - muons: mi+, mi-" <<endl;
  cout<< "1 - electrons & positons: e+, e-" <<endl;
  cout<< "2 - photons: gamma" <<endl;

  cout<< "\nINPUT: IDs of particles of cosmic background (default " << I_criterium[0] << " - " << I_criterium[id_len - 1] << ") : " <<endl;
											//Inputing I_criterium
  int id_iter = 0;

  getline(cin, input);
  while(!input.empty()) {
    istringstream stream(input);
    stream >> I_criterium[id_iter];

    getline(cin, input);

    id_iter++;
    if(id_iter >= id_len) break;
  }

  if(id_iter > 0) for(int id = id_iter; id < id_len; id++) I_criterium[id] = 0;

  //-------------------------------------------------------------

  /*	This part of code runs previously
	created macro.
  */

  CascadeSignalEstimation(I_criterium, params_dir, plots_dir, E_max, R_maximum, e_max, Aprime, eff, n, Tdet, f_fake, Tp, fi, eps, eps_phi, n_rand);	//Here the macro CascadeSignalEstimation is kick-started

  //-------------------------------------------------------------

  auto stop = high_resolution_clock::now(); 
  auto duration = duration_cast<microseconds>(stop - start); 			//Stoping counting time of computations

  cout<<"\n" <<endl;
  cout<< "Time of the computations: " << duration.count() << " * 10^{-6} s = " <<  duration.count()/pow(10, 6) << " s = "  << duration.count()/(pow(10, 6)*60) << " min" <<endl;

  //-------------------------------------------------------------

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~				//End MACRO


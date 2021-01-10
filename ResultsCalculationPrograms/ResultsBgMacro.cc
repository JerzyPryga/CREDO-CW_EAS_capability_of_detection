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

#include "BgSignalEstimation.cc"						//inclding CREDO macro

using namespace std::chrono;
using namespace std;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*	Macro with all inputs included which
	executes previously written macros.
*/

void ResultsBgMacro() {								//Start MACRO

  //-------------------------------------------------------------

  auto start = high_resolution_clock::now(); 					//Starting counting time of computations

  //-------------------------------------------------------------

  /*	Here the starting parameters for the program
	are inputed. 

	The meaning of each parameter:

	- I_criterium - list of particles IDs,

	- Aprime - surface area of detectors [cm^2] (default = 25),
        - eff - efficiency of the detector [%] (default = 95),
	- n - number of detectors (default = 4),
	- Tdet - detector's registration time [s] (default = 10^(-7)),
	- f_fake - frequency of fake signals from the device [1/s] (default = 0.1),

	- Tp - time of measurement [s] (default = 3600*24*7 = 1 week),
	- fi - maximum vertical angle of arriving particles [degrees] (default = 90).
	
	When user clicks ENTER the
	default value is set.
  */
  
  int Tp, n;
  int I_criterium[3] = {0, 1, 2};
  double Aprime, eff, Tdet, f_fake, fi;

  int id_len = 3;

  string input;

  Tp = 7;
  n = 4;

  Aprime = 25.0;
  eff = 0.95;
  Tdet = 2*pow(10, -7);
  f_fake = 0.1;
  fi = 90.0;


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
  cout<< "CCOSMIC RAYS ASSUMPTIONS" <<endl;

  cout<< "INPUT: maximum vertical angle fi [degrees] (default = 90) : ";			//Inputing fi

  getline(cin, input);
  if(!input.empty()) {
    istringstream stream(input);
    stream >> fi;
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

  BgSignalEstimation(I_criterium, Aprime, eff, n, Tdet, f_fake, Tp, fi);	//Here the macro BgSignalEstimation is kick-started


  //-------------------------------------------------------------

  auto stop = high_resolution_clock::now(); 
  auto duration = duration_cast<microseconds>(stop - start); 			//Stoping counting time of computations

  cout<<"\n" <<endl;
  cout<< "Time of the computations: " << duration.count() << " * 10^{-6} s = " <<  duration.count()/pow(10, 6) << " s = "  << duration.count()/(pow(10, 6)*60) << " min" <<endl;

  //-------------------------------------------------------------

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~				//End MACRO


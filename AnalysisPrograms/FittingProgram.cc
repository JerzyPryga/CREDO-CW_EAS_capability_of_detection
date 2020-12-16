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
#include "TF1.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TFitResult.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TLeafI.h"
#include "TFile.h"

#include "MacroNRe.cc"								//inclding CREDO macros
#include "MacroRhoNE.cc"							//inclding CREDO macros
#include "MacroRhoRTheta.cc"							//inclding CREDO macros

using namespace std::chrono;
using namespace std;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/* 	Used functions */

string get_current_dir();							//Function returning path to current directory

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*	The main C++ program is here and
	it just reads parameters from user
	and executes previously written macros.
*/

int main() {									//Start PROGRAM 

/* 	Here the starting parameters for the program
	are inputed.

	The meaning of each starting parameters of macro:

	- ID_criterium - list of particles IDs,

	- sim_dir_vertical - path to directory with files containing vertical CORSIKA simulations,
	- sim_dir_theta - path to directory with files containing angular CORSIKA simulations,
	- plots_dir - path to directory where all plots are saved (default = current directory),
        - norm_file - name of the file for parametrization (default = last one found),

	- r_number - number of points on axis of distance from cascade center r (default = 10)

	- p_part_min - minimum momenta of wanted particles [GeV] (default 1),
	- p_part_max - maximum momenta of wanted particles [GeV] (default pow(10, 10) ~ inf),

	- bg - background level [parts/m^2 s^2] (default = 70),
	- prc - percent of particles in radius R_prc (default = 95),
	- Tdet - detector's registration time Tdet [s] (default = pow(10, -7)),

	- bins_n - number of bins in histograms of other quantities (default = 40),
	- bins_N_part - number of bins in histograms of number of particles N_part (default = 40),
	- bins_r - number of bins in histogram of particles density Rho(r) (default = 40),
	- bins_R_prc - number of bins in histogram of R_prc (default = 40),
	- bins_R_rho - number of bins in histogram of R_rho (default = 40).

	- max_theta - maximum angle of primary particle [degrees] (default = 70).

	- n_datas_vertical - number of vertical simulations (default = 18).
	- n_datas_theta - number of angular simulations (default = 8).

	When user clicks ENTER the
	default value is set.
*/

  //-------------------------------------------------------------

  /*	Here the current directory is found. */

  string current_dir = get_current_dir();
  string input;

  //-------------------------------------------------------------

  cout<<endl;
  cout<< "---------- Number of particles | Maximum Radius ----------" <<endl;
  cout<< "Do You want to perform analysis of:" <<endl;
  cout<< " - Number of particles: N_part," <<endl;
  cout<< " - Radius for given percent of particles: R_prc," <<endl;
  cout<< " - Radius for given minimal particles density: R_rho." <<endl;
  cout<< " Press ENTER - YES" <<endl;
  cout<< " Type anything - NO" <<endl;

  bool NRe;

  getline(cin, input);
  if(!input.empty()) NRe = false;
  else NRe = true;

  input = "";

  if(NRe == true) MacroNRe(current_dir);		//Start NRe analysis

  //-------------------------------------------------------------

  cout<<endl;
  cout<< "---------- Density | Particles number | Energy ----------" <<endl;
  cout<< "Do You want to perform analysis of:" <<endl;
  cout<< " - Particles density as a function of distance r: Rho(r)," <<endl;
  cout<< " - Particles density as a function of energy E: Rho(E)," <<endl;
  cout<< " - Particles density as a function of number of particles N_part r: Rho( N_part)," <<endl;
  cout<< " and combination of above." <<endl;
  cout<< " Press ENTER - YES" <<endl;
  cout<< " Type anything - NO" <<endl;

  bool RhoNE;

  getline(cin, input);
  if(!input.empty()) RhoNE = false;
  else RhoNE = true;

  input = "";

  if(RhoNE == true) MacroRhoNE(current_dir);		//Start RhoNE analysis

  //-------------------------------------------------------------

  cout<<endl;
  cout<< "---------- Angle ----------" <<endl;
  cout<< "Do You want to perform analysis of:" <<endl;
  cout<< " - Particles density as a function of distance r: Rho(r)," <<endl;
  cout<< " - Particles density as a function of primary particle angle Theta: Rho(Theta)," <<endl;
  cout<< " and combination of above." <<endl;
  cout<< " Press ENTER - YES" <<endl;
  cout<< " Type anything - NO" <<endl;

  bool RhoRTheta;

  getline(cin, input);
  if(!input.empty()) RhoRTheta = false;
  else RhoRTheta = true;

  input = "";

  if(RhoRTheta == true) MacroRhoRTheta(current_dir);	//Start RhoRTheta analysis

  //-------------------------------------------------------------

  cout<<endl;
  cout<< "END OF PROGRAM" <<endl;

  return 0;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~				//End PROGRAM

/*	Used functions definitions
*/

std::string get_current_dir() {							//Function returning path to current directory

   char buff[FILENAME_MAX]; 							
   GetCurrentDir( buff, FILENAME_MAX );
   std::string current_working_dir(buff);
   return current_working_dir;

}

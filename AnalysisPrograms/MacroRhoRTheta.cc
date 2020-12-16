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

#include "Rho_theta_fit_macro.cc"						//inclding CREDO macros
#include "Rho_theta_macro.cc"							//inclding CREDO macros

using namespace std::chrono;
using namespace std;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*	The main C++ program is here and
	it just reads parameters from user
	and executes previously written macros.
*/

void MacroRhoRTheta(string current_dir) {					//Start PROGRAM

  //-------------------------------------------------------------

  auto start = high_resolution_clock::now(); 					//Starting counting time of computations

  //-------------------------------------------------------------

/* 	Here the starting parameters for the program
	are inputed.

	The meaning of each starting parameters of macro:

	- ID_criterium - list of particles IDs,

	- sim_dir - path to directory with files containing CORSIKA simulations,
	- plots_dir - path to directory where all plots are saved (default = current directory),

	- bins_rho_N_part - number of bins in histograms of number of particles N_part (default = 40),
	- bins_rho_r - number of bins in histogram of particles density Rho(r) (default = 40),

	- r_number - number of points on axis of distance from cascade center r (default = 10)

	- p_part_min - minimum momenta of wanted particles [GeV] (default 1),
	- p_part_max - maximum momenta of wanted particles [GeV] (default pow(10, 10) ~ inf),

	- n_datas - number of simulations (default = 8).

	When user clicks ENTER the
	default value is set.
*/
  
  string sim_dir = "/home/jerzy/CREDO/Analiza/Angle_distribution/Simulations";
  string plots_dir = current_dir;

  int bins_rho_r, bins_rho_N_part, r_number, ID_min, ID_max;
  int ID_criterium[9] = {5, 6, 0, 0, 0 ,0 ,0 ,0, 0};
  double p_part_min, p_part_max, max_theta;

  string input;
  int id_len = sizeof(ID_criterium)/sizeof(ID_criterium[0]);

  bins_rho_r = 40;
  bins_rho_N_part = 40;

  r_number = 10;

  p_part_min = 1.0;
  p_part_max = pow(10, 10);

  max_theta = 70.0;


  cout<<endl;
  cout<< "DIRECTORIES FOR SIMULATIONS" <<endl;

  cout<< "INPUT: path to directory with files containing CORSIKA simulations : ";		//Inputing sim_dir

  getline(cin, input );
  if(!input.empty()) {
    istringstream stream(input);
    stream >> sim_dir;
  }




  cout<<endl;
  cout<< "DIRECTORIES FOR PLOTS" <<endl;

  cout<< "INPUT: path to directory where all plots are saved (default: current directory) : ";	//Inputing plots_dir
	
  getline(cin, input );
  if(!input.empty()) {
    istringstream stream(input);
    stream >> plots_dir;
  }




  cout<<endl;
  cout<< "DRAWING PARAMETERS" <<endl;

  cout<< "INPUT: number of bins in histogram of number of particles N_part (default = 40) : ";	//Inputing bins_rho_N_part

  getline(cin, input );
  if(!input.empty()) {
    istringstream stream(input);
    stream >> bins_rho_N_part;
  }

  cout<< "INPUT: number of bins in histogram of particles density Rho(r) (default = 40) : ";	//Inputing bins_rho_r

  getline(cin, input );
  if(!input.empty()) {
    istringstream stream(input);
    stream >> bins_rho_r;
  }

  bins_rho_r = bins_rho_r + 1;									//Adding one due to different iterations in
  bins_rho_N_part = bins_rho_N_part + 1;							//arrays and histograms


  cout<<endl;
  cout<< "FITTING PARAMETER" <<endl;

  cout<< "INPUT: number of distance points on X-axis for fitting f(r) (default = 10) : ";	//Inputing bins_R_rho

  getline(cin, input );
  if(!input.empty()) {
    istringstream stream(input);
    stream >> r_number;
  }




  cout<<endl;
  cout<< "CRITERIA OF PARTICLE CHOOSING" <<endl;

  cout<< "\nCorsica ID | Particle" <<endl;
  cout<< "1 - Photon gamma" <<endl;
  cout<< "2 - Positon e+" <<endl;
  cout<< "3 - Electron e-" <<endl;

  cout<< "\n5 - Muon mi+" <<endl;
  cout<< "6 - Muon mi-" <<endl;
  cout<< "7 - Pi meson pi0 " <<endl;
  cout<< "8 - Pi meson pi+" <<endl;
  cout<< "9 - Pi meson pi-" <<endl;

  cout<< "\nWhole list of particles IDs can be found on: " <<endl;
  cout<< "https://web.ikp.kit.edu/corsika/usersguide/usersguide.pdf" <<endl;

  cout<< "\nINPUT: IDs of particles in the shower (default " << ID_criterium[0] << " & " << ID_criterium[1] << ") : " <<endl;		//Inputing ID_min

  int id_iter = 0;

  getline(cin, input);
  while(!input.empty()) {
    istringstream stream(input);
    stream >> ID_criterium[id_iter];

    getline(cin, input);

    id_iter++;
    if(id_iter >= id_len) break;
  }

  if(id_iter > 0) for(int id = id_iter; id < id_len; id++) ID_criterium[id] = 0;

  cout<< "INPUT: minimum momentum p_part of particles of the shower [GeV] (default " << p_part_min << " GeV) : ";		//Inputing p_part_min

  getline(cin, input );
  if(!input.empty()) {
    istringstream stream(input);
    stream >> p_part_min;
  }

  cout<< "INPUT: maximum momentum p_part of particles of the shower [GeV] (default " << p_part_max << " GeV ~ inf) : ";		//Inputing p_part_max

  getline(cin, input );
  if(!input.empty()) {
    istringstream stream(input);
    stream >> p_part_max;
  }

  cout<< "INPUT: maximum angle of primary particle of the shower [degrees] (default " << max_theta << " degrees) : ";		//Inputing max_theta

  getline(cin, input );
  if(!input.empty()) {
    istringstream stream(input);
    stream >> max_theta;
  }


  DIR *dir_plots;
  struct dirent *ent_plots;

  if((dir_plots = opendir(plots_dir.c_str())) == NULL) {			//Checking the existence of directory for plots
    plots_dir = current_dir;
    cout<< "\nERROR: PLOTS DIRECTORY SET TO CURRENT DIRECTORY" <<endl;
  }

  //-------------------------------------------------------------

  /*	This part of code counts how many
	files (not directories) are in
	designated directory.
  */
  
  int n_datas = 0;

  DIR *dir;
  struct dirent *ent;

  if((dir = opendir(sim_dir.c_str())) != NULL) {				//Opening the directory
    while((ent = readdir(dir)) != NULL) {					//Checking if there is something
      if(opendir(ent->d_name) == NULL) {					//Checking if it is a dirtectory

        n_datas++;								//Counting the files

      }
    }
    cout<<endl;
    cout<< "NUMBER OF FILES FOUND: " << n_datas <<endl;

    closedir(dir);

    Rho_theta_macro(ID_criterium, sim_dir, plots_dir, bins_rho_N_part, bins_rho_r, r_number, p_part_min, p_part_max, n_datas, max_theta);
										//Here the macro Rho_theta_macro is kick-started
    Rho_theta_fit_macro(plots_dir, r_number);					//Here the macro Rho_theta_fit_macro is kick-started

  } 

  else {
    cout<< "ERROR: COULD NOT OPEN THE SIMULATIONS DIRECTORY" <<endl;
  }

  //-------------------------------------------------------------

  auto stop = high_resolution_clock::now(); 
  auto duration = duration_cast<microseconds>(stop - start); 			//Stoping counting time of computations

  cout<<"\n" <<endl;
  cout<< "Time of the computations: " << duration.count() << " * 10^{-6} s = " <<  duration.count()/pow(10, 6) << " s = "  << duration.count()/(pow(10, 6)*60) << " min" <<endl;

  //-------------------------------------------------------------

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~				//End PROGRAM


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


#include "NR_E_fit_macro.cc"							//inclding CREDO macros
#include "NR_E_macro.cc"							//inclding CREDO macros

#include "Rho_rNE_fit_macro.cc"							//inclding CREDO macros
#include "Rho_rNE_macro.cc"							//inclding CREDO macros

#include "Rho_theta_fit_macro.cc"						//inclding CREDO macros
#include "Rho_theta_macro.cc"							//inclding CREDO macros

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

  //-------------------------------------------------------------
  
  string sim_dir_vertical = "/home/jerzy/CREDO/Analiza/More_simulations";
  string sim_dir_theta = "/home/jerzy/CREDO/Analiza/Angle_distribution/Simulations";
  string plots_dir = current_dir;
  string norm_file = "";

  int bins_n, bins_N_part, bins_r, bins_R_prc, bins_R_rho, r_number;
  int ID_criterium[9] = {5, 6, 0, 0, 0 ,0 ,0 ,0, 0};
  double p_part_min, p_part_max, max_theta;
  double bg, prc, Tdet;
  int n_datas_vertical, n_datas_theta;

  string input;
  int id_len = sizeof(ID_criterium)/sizeof(ID_criterium[0]);

  p_part_min = 1.0;
  p_part_max = pow(10, 10);

  bg = 70.0;
  prc = 95;
  Tdet = pow(10, -7);

  r_number = 10;

  bins_n = 40;
  bins_N_part = 40;
  bins_r = 40;
  bins_R_prc = 40;
  bins_R_rho = 40;

  max_theta = 70.0;

  bool NRe = true;
  bool RhoNE = true;
  bool RhoRTheta = true;

  //-------------------------------------------------------------

  cout<<endl;
  cout<< "---------- Number of particles | Maximum Radius ----------" <<endl;
  cout<< "Do You want to perform analysis of:" <<endl;
  cout<< " - Number of particles: N_part," <<endl;
  cout<< " - Radius for given percent of particles: R_prc," <<endl;
  cout<< " - Radius for given minimal particles density: R_rho." <<endl;
  cout<< " Press ENTER - YES" <<endl;
  cout<< " Type anything - NO" <<endl;

  getline(cin, input );
  if(!input.empty()) {
    NRe = false;
  }


  if(NRe == true) {					//Start NRe analysis

  //-------------------------------------------------------------

  cout<<endl;
  cout<< "DIRECTORIES FOR SIMULATIONS" <<endl;

  cout<< "INPUT: path to directory with files containing CORSIKA simulations : ";		//Inputing sim_dir

  getline(cin, input);
  if(!input.empty()) {
    istringstream stream(input);
    stream >> sim_dir_vertical;
  }


  cout<<endl;
  cout<< "DIRECTORIES FOR PLOTS" <<endl;

  cout<< "INPUT: path to directory where all plots are saved (default: current directory) : ";	//Inputing plots_dir
	
  getline(cin, input);
  if(!input.empty()) {
    istringstream stream(input);
    stream >> plots_dir;
  }




  cout<<endl;
  cout<< "ASSUMPTIONS PARAMETERS" <<endl;

  cout<< "INPUT: background level [parts/m^2 s^2] (default = 70) : ";				//Inputing bg

  getline(cin, input);
  if(!input.empty()) {
    istringstream stream(input);
    stream >> bg;
  }

  cout<< "INPUT: percent of particles in radius R_prc (default = 95) : ";			//Inputing prc

  getline(cin, input);
  if(!input.empty()) {
    istringstream stream(input);
    stream >> prc;
  }

  cout<< "INPUT: detector's registration time Tdet [s] (default = pow(10, -7)) : ";		//Inputing Tdet

  getline(cin, input);
  if(!input.empty()) {
    istringstream stream(input);
    stream >> Tdet;
  }




  cout<<endl;
  cout<< "DRAWING PARAMETERS" <<endl;

  cout<< "INPUT: number of bins in histogram of number of particles N (default = 40) : ";	//Inputing bins_N_part

  getline(cin, input);
  if(!input.empty()) {
    istringstream stream(input);
    stream >> bins_N_part;
  }

  cout<< "INPUT: number of bins in histogram of R_prc (default = 40) : ";			//Inputing bins_R_prc

  getline(cin, input);
  if(!input.empty()) {
    istringstream stream(input);
    stream >> bins_R_prc;
  }

  cout<< "INPUT: number of bins in histogram of R_rho (default = 40) : ";			//Inputing bins_R_rho

  getline(cin, input);
  if(!input.empty()) {
    istringstream stream(input);
    stream >> bins_R_rho;
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

  cout<< "\nINPUT: minimum IDs of particles of the shower (default " << ID_criterium[0] << " & " << ID_criterium[1] << ") : " <<endl;		//Inputing ID array

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

  if((dir = opendir(sim_dir_vertical.c_str())) != NULL) {			//Opening the directory
    while((ent = readdir(dir)) != NULL) {					//Checking if there is something
      if(opendir(ent->d_name) == NULL) {					//Checking if it is a dirtectory

        n_datas++;								//Counting the files

      }
    }
    cout<<endl;
    cout<< "NUMBER OF FILES FOUND: " << n_datas <<endl;

    closedir(dir);

    NR_E_macro(ID_criterium, sim_dir_vertical, plots_dir, prc, bins_N_part, bins_R_prc, bins_R_rho, bg, Tdet, p_part_min, p_part_max, n_datas);
										//Here the macro NR_E is kick-started
    NR_E_fit_macro(plots_dir);							//Here the macro NR_E_fit is kick-started

  } 

  else {
    cout<< "ERROR: COULD NOT OPEN THE SIMULATIONS DIRECTORY" <<endl;
  }

  //-------------------------------------------------------------

  }							//Stop NRe analysis

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

  getline(cin, input );
  if(!input.empty()) {
    RhoNE = false;
  }


  if(RhoNE == true) {					//Start RhoNE analysis

  //-------------------------------------------------------------
  
  cout<<endl;
  cout<< "DIRECTORIES FOR SIMULATIONS" <<endl;

  cout<< "INPUT: path to directory with files containing CORSIKA simulations : ";		//Inputing sim_dir_vertical

  getline(cin, input );
  if(!input.empty()) {
    istringstream stream(input);
    stream >> sim_dir_vertical;
  }




  cout<<endl;
  cout<< "FILE FOR PARAMETRISATION" <<endl;

  cout<< "INPUT: name of .root file which contain simulation of cascades " <<endl;
  cout<< "for which fitted functions will be parametrised " <<endl;
  cout<< "(default: the last file read from the dictionary) :";					//Inputing norm_file
	
  getline(cin, input );
  if(!input.empty()) {
    istringstream stream(input);
    stream >> norm_file;
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

  cout<< "INPUT: number of bins in histogram of partices density Rho(r) (default = 40) : ";	//Inputing bins_N_part

  getline(cin, input );
  if(!input.empty()) {
    istringstream stream(input);
    stream >> bins_N_part;
  }

  cout<< "INPUT: number of bins in histogram of relative partices density Rho(N)/Rho(N0) (default = 40) : ";	//Inputing bins_N_part

  getline(cin, input );
  if(!input.empty()) {
    istringstream stream(input);
    stream >> bins_N_part;
  }

  bins_r = bins_r + 1;										//Adding one due to different iterations in
  bins_N_part = bins_N_part + 1;								//arrays and histograms


  cout<<endl;
  cout<< "FITTING PARAMETER" <<endl;

  cout<< "INPUT: number of distance points on X-axis for fitting f(r) (default = 10) : ";	//Inputing r_number

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
   
  n_datas_vertical = 0;

  DIR *dir;
  struct dirent *ent;

  if((dir = opendir(sim_dir_vertical.c_str())) != NULL) {			//Opening the directory
    while((ent = readdir(dir)) != NULL) {					//Checking if there is something
      if(opendir(ent->d_name) == NULL) {					//Checking if it is a dirtectory

        n_datas_vertical++;							//Counting the files

      }
    }
    cout<<endl;
    cout<< "NUMBER OF FILES FOUND: " << n_datas_vertical <<endl;

    closedir(dir);

    Rho_rNE_macro(ID_criterium, sim_dir_vertical, plots_dir, norm_file, bins_N_part, bins_r, r_number, p_part_min, p_part_max, n_datas_vertical);
										//Here the macro Rho_rNE is kick-started
    Rho_rNE_fit_macro(plots_dir, r_number);					//Here the macro Rho_rNE_fit is kick-started

  } 

  else {
    cout<< "ERROR: COULD NOT OPEN THE SIMULATIONS DIRECTORY" <<endl;
  }
  
  //-------------------------------------------------------------

  }							//Stop RhoNE analysis

  //-------------------------------------------------------------

  cout<<endl;
  cout<< "---------- Angle ----------" <<endl;
  cout<< "Do You want to perform analysis of:" <<endl;
  cout<< " - Particles density as a function of distance r: Rho(r)," <<endl;
  cout<< " - Particles density as a function of primary particle angle Theta: Rho(Theta)," <<endl;
  cout<< " and combination of above." <<endl;
  cout<< " Press ENTER - YES" <<endl;
  cout<< " Type anything - NO" <<endl;

  getline(cin, input );
  if(!input.empty()) {
    RhoRTheta = false;
  }


  if(RhoRTheta == true) {				//Start RhoRTheta analysis

  //-------------------------------------------------------------

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

  cout<< "INPUT: number of bins in histogram of number of particles N_part (default = 40) : ";	//Inputing bins_N_part

  getline(cin, input );
  if(!input.empty()) {
    istringstream stream(input);
    stream >> bins_N_part;
  }

  cout<< "INPUT: number of bins in histogram of particles density Rho(r) (default = 40) : ";	//Inputing bins_r

  getline(cin, input );
  if(!input.empty()) {
    istringstream stream(input);
    stream >> bins_r;
  }

  bins_r = bins_r + 1;										//Adding one due to different iterations in
  bins_N_part = bins_N_part + 1;								//arrays and histograms


  cout<<endl;
  cout<< "FITTING PARAMETER" <<endl;

  cout<< "INPUT: number of distance points on X-axis for fitting f(r) (default = 10) : ";	//Inputing r_number

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

  cout<< "\nINPUT: IDs of particles in the shower (default " << ID_criterium[0] << " & " << ID_criterium[1] << ") : " <<endl;	//Inputing IDs

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

  cout<< "INPUT: maximum angle of primary particle of the shower [degrees] (default " << max_theta << " degrees : ";		//Inputing max_theta

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
  
  n_datas_theta = 0;

  DIR *dir;
  struct dirent *ent;

  if((dir = opendir(sim_dir_theta.c_str())) != NULL) {				//Opening the directory
    while((ent = readdir(dir)) != NULL) {					//Checking if there is something
      if(opendir(ent->d_name) == NULL) {					//Checking if it is a dirtectory

        n_datas_theta++;							//Counting the files

      }
    }
    cout<<endl;
    cout<< "NUMBER OF FILES FOUND: " << n_datas_theta <<endl;

    closedir(dir);

    Rho_theta_macro(ID_criterium, sim_dir_theta, plots_dir, bins_N_part, bins_r, r_number, p_part_min, p_part_max, n_datas_theta, max_theta);
										//Here the macro Rho_theta_macro is kick-started
    Rho_theta_fit_macro(plots_dir, r_number);					//Here the macro Rho_theta_fit_macro is kick-started

  } 

  else {
    cout<< "ERROR: COULD NOT OPEN THE SIMULATIONS DIRECTORY" <<endl;
  }

  //-------------------------------------------------------------

  }							//Stop RhoTheta analysis

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

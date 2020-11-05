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

using namespace std::chrono;
using namespace std;

//-------------------------------------------------------------

/* 	Used functions */

string get_current_dir();							//Function returning path to current directory

//-------------------------------------------------------------

/* 	The majority of the code is
	a ROOT macro called "NR_E". All computations
	are performed by this macro.
*/

void NR_E()
{

  //-------------------------------------------------------------

  auto start = high_resolution_clock::now(); 					//Starting counting time of computations

  //-------------------------------------------------------------  

  /*	Here the current directory is found. */

  string current_dir = get_current_dir();

  //-------------------------------------------------------------

  /*	Here the starting parameters for the program
	are inputed. 

	The meaning of each parameter:

	- sim_dir - path to directory with files containing CORSIKA simulations,
	- plots_dir - path to directory where all plots are saved,

	- bg - background level [parts/m^2 s^2] (default = 70),
	- prc - percent of particles in radius R_prc (default = 95),
	- Tdet - detector's registration time Tdet [s] (default = pow(10, -7)),

	- bins_N_part - number of bins in histogram of number of particles N (default = 40),
	- bins_R_prc - number of bins in histogram of R_prc (default = 40),
	- bins_R_rho - number of bins in histogram of R_rho (default = 40).

	- ID_min - minimum ID of wanted particles (defoult 5),
	- ID_max - maximum ID of wanted particles (defoult 6),
	- p_part_min - minimum momenta of wanted particles [GeV] (default 1),
	- p_part_max - maximum momenta of wanted particles [GeV] (default pow(10, 10) ~ inf).
	
	When user clicks ENTER the
	default value is set.
  */

  string sim_dir = "/home/jerzy/CREDO/Analiza/More_simulations";
  string plots_dir = current_dir;

  int prc, bins_N_part, bins_R_prc, bins_R_rho;
  int ID_criterium[9] = {5, 6, 0, 0, 0 ,0 ,0 ,0, 0};
  double bg, Tdet, p_part_min, p_part_max;

  string input;

  int id_len = sizeof(ID_criterium)/sizeof(ID_criterium[0]);

  bg = 70.0;
  prc = 95;
  Tdet = pow(10, -7);

  bins_N_part = 40;
  bins_R_prc = 40;
  bins_R_rho = 40;

  p_part_min = 1.0;
  p_part_max = pow(10, 10);




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
  cout<< "ASSUMPTIONS PARAMETERS" <<endl;

  cout<< "INPUT: background level [parts/m^2 s^2] (default = 70) : ";				//Inputing bg

  getline(cin, input );
  if(!input.empty()) {
    istringstream stream(input);
    stream >> bg;
  }

  cout<< "INPUT: percent of particles in radius R_prc (default = 95) : ";			//Inputing prc

  getline(cin, input );
  if(!input.empty()) {
    istringstream stream(input);
    stream >> prc;
  }

  cout<< "INPUT: detector's registration time Tdet [s] (default = pow(10, -7)) : ";		//Inputing Tdet

  getline(cin, input );
  if(!input.empty()) {
    istringstream stream(input);
    stream >> Tdet;
  }




  cout<<endl;
  cout<< "DRAWING PARAMETERS" <<endl;

  cout<< "INPUT: number of bins in histogram of number of particles N (default = 40) : ";	//Inputing bins_N_part

  getline(cin, input );
  if(!input.empty()) {
    istringstream stream(input);
    stream >> bins_N_part;
  }

  cout<< "INPUT: number of bins in histogram of R_prc (default = 40) : ";			//Inputing bins_R_prc

  getline(cin, input );
  if(!input.empty()) {
    istringstream stream(input);
    stream >> bins_R_prc;
  }

  cout<< "INPUT: number of bins in histogram of R_rho (default = 40) : ";			//Inputing bins_R_rho

  getline(cin, input );
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

  cout<< "\nINPUT: minimum IDs of particles of the shower (default " << ID_criterium[0] << " & " << ID_criterium[1] << ") : " <<endl;		//Inputing ID_min

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

  if((dir = opendir(sim_dir.c_str())) != NULL) {				//Opening the directory
    while((ent = readdir(dir)) != NULL) {					//Checking if there is something
      if(opendir(ent->d_name) == NULL) {					//Checking if it is a dirtectory

        n_datas++;								//Counting the files

      }
    }
    cout<<endl;
    cout<< "NUMBER OF FILES FOUND: " << n_datas <<endl;

    closedir(dir);
  } 

  else {
    cout<< "ERROR: COULD NOT OPEN THE SIMULATIONS DIRECTORY" <<endl;
  }

  //-------------------------------------------------------------

  /*	Here are only definitions of some
	arrays and variables used later.
  */

  double meanR_prc_lst[n_datas], sigmaR_prc_lst[n_datas], cR_prc_lst[n_datas];		//An arrays of parameters of fitted functions: meanR_prc, sigmaR_prc and cR_prc
  double meanR_rho_lst[n_datas], sigmaR_rho_lst[n_datas], cR_rho_lst[n_datas];		//An arrays of parameters of fitted functions: meanR_rho, sigmaR_rho and cR_rho
  double meanN_part_lst[n_datas], sigmaN_part_lst[n_datas], cN_part_lst[n_datas];	//An arrays of parameters of fitted functions: meanN_part, sigmaN_part and cN_part	
  double E_lst[n_datas]; 								//An array of energies

  double E_min;									//Minimum E
  double E_max;									//Maximum E
  double R_prc_maximum = 0;							//Maximum R_prc - global
  double R_rho_maximum = 0;							//Maximum R_rho - global
  int N_part_maximum = 0;							//Maximum N_part - global
  

  //-------------------------------------------------------------

  /*	Here are defined some directories
	in which plots are stored.
  */

  string N_part_E_dir_name = "/N_part(E)";					//Naming directories for plots for energy E relationship
  string N_part_E_dir = plots_dir;
  N_part_E_dir.append(N_part_E_dir_name);
  
  string R_prc_E_dir_name = "/R_prc(E)";
  string R_prc_E_dir = plots_dir;
  R_prc_E_dir.append(R_prc_E_dir_name);

  string R_rho_E_dir_name = "/R_rho(E)";
  string R_rho_E_dir = plots_dir;
  R_rho_E_dir.append(R_rho_E_dir_name);


  string RN_prc_dir_name = "/n(N_part, R_prc)";					//Naming directories for histograms
  string RN_prc_dir = plots_dir;
  RN_prc_dir.append(RN_prc_dir_name);

  string N_part_dir_name = "/n(N_part)";
  string N_part_dir = plots_dir;
  N_part_dir.append(N_part_dir_name);
  
  string R_prc_dir_name = "/n(R_prc)";
  string R_prc_dir = plots_dir;
  R_prc_dir.append(R_prc_dir_name);

  string R_rho_dir_name = "/n(R_rho)";
  string R_rho_dir = plots_dir;
  R_rho_dir.append(R_rho_dir_name);


  int check_N_part_E = mkdir(N_part_E_dir.c_str(), 0777);			//Creating directories for plots for energy E relationship
  int check_R_prc_E = mkdir(R_prc_E_dir.c_str(), 0777);
  int check_R_rho_E = mkdir(R_rho_E_dir.c_str(), 0777);

  int check_N_part = mkdir(N_part_dir.c_str(), 0777);				//Creating directories for histograms
  int check_R_prc = mkdir(R_prc_dir.c_str(), 0777);
  int check_R_rho = mkdir(R_rho_dir.c_str(), 0777);
  int check_RN_prc = mkdir(RN_prc_dir.c_str(), 0777);

  //-------------------------------------------------------------

  /*	This part of code calculate level of background.
	Further in the script rho_R is evaluated which is
	a distance from the center of the cascade to the point
	in which partcile density from the shower drops to
	the level of background density.
  */ 

  double fi = M_PI/2;

  TF1 f_bg("f_bg","pow(cos(x), 3)*sin(x)", 0, fi);				//Function of background particles distribution
										//and geometrical property of the detector over zenith angle - theta
  double integral_bg = f_bg.Integral(0,fi);					//Integral over all angles up to 90 degrees

  double Im0 = bg;
  double Isr = Im0*2*M_PI*integral_bg;
 
  cout<<endl;
  cout<< "BACKGROUND" <<endl;
  cout<<endl;
  cout<< "Background flux = "<< Isr << " m^(-2)s^(-1)" <<endl;

  double rho_bg = (Isr*Tdet/10000);
  cout<< "Background density (for Tdet = " << Tdet << ") Rho_bg = " << rho_bg << " cm^(-2)" <<endl;				//Printing background level


  //-------------------------------------------------------------


  int nf = 0;									//Files iterator

  if((dir = opendir(sim_dir.c_str())) != NULL) {				
  while((ent = readdir(dir)) != NULL) {						//A loop over all files in directory - Start FILES LOOP
  if(opendir(ent->d_name) == NULL) {

  //-------------------------------------------------------------

  /*	This part of code reads each file,
	finds energy of primary particle for
	simulated cascades and number of simulations
	in the file.
  */


  cout<< "\nFile name: " <<endl;
  cout<< ent->d_name <<endl;

  string file_sim = sim_dir;
  file_sim.append("/");						
  file_sim.append(ent->d_name);							//Adding name of the file into directory path

  cout<< "\nOpening file nr " << nf+1 << " : " <<endl;
  cout<< file_sim <<endl;

  TFile *file = new TFile(file_sim.c_str(), "READ");				//Opening each file
  TTree *sim = (TTree*) file->Get("sim");

  
  TLeafI* l_id = (TLeafI*)sim->GetLeaf("particle..ParticleID");			//Reading needed data
  TLeaf* l_x = (TLeaf*)sim->GetLeaf("particle..x");
  TLeaf* l_y = (TLeaf*)sim->GetLeaf("particle..y");
  TLeaf* l_pz = (TLeaf*)sim->GetLeaf("particle..Pz");
  TLeaf* l_px = (TLeaf*)sim->GetLeaf("particle..Px");
  TLeaf* l_py = (TLeaf*)sim->GetLeaf("particle..Py");

  TLeaf* E_sim = (TLeaf*)sim->GetLeaf("shower.Energy");

  //-------------------------------------------------------------

  /*	This part of code reads energy of
	the cascades in each file
	(assuming that all showers within the file
	have the same energy). It also finds
	its maximum and minimal values E_min and E_max.
  */

  sim->GetEntry(0);
  //E_lst[nf] = E_sim->GetValue(0);						//Filling an array with energy in GeV
  E_lst[nf] = E_sim->GetValue(0)/1000.0;					//Filling an array with energy in TeV

  if(nf == 0) {
    E_min = E_lst[nf];
    E_max = E_lst[nf];
  }
  if(E_lst[nf] <= E_min) E_min = E_lst[nf];					//Finding minimal energy E_min

  if(E_lst[nf] >= E_max) E_max = E_lst[nf];					//Finding maximal energy E_max


  //cout<< "\nEnergy of primary particle: " << E_lst[nf] << " GeV" <<endl;	//Printing energy in GeV 
  cout<< "\nEnergy of primary particle: " << E_lst[nf] << " TeV" <<endl;	//Printing energy in TeV

  //-------------------------------------------------------------

  /*	Here, simulations are counted.
  */

  int n_sim = 0;
  while(sim->GetEntry(n_sim) != 0) {
    n_sim++;
  }										//Counting number of simulations in the file		

  cout<< "Number of simulations in the shower: " << n_sim <<endl;		//Printing number of simulations in the file
  cout<< "\n" <<endl;

  //-------------------------------------------------------------

  /*
	There are some definitions of
	variables and lists used later.
  */

  long double R_prc_temp_lst[n_sim];						//An array of R_prc for each simulation in the file
  double R_prc_max = 0;								//Maximum R_prc - temporary

  long double R_rho_temp_lst[n_sim];						//An array of R_rho for each simulation in the file
  double R_rho_max = 0;								//Maximum R_rho - temporary


  int N_part_temp_lst[n_sim];							//An array of N for each simulation in the file
  int N_part_max = 0;								//Maximum N_part - temporary

  //-------------------------------------------------------------

  for(int n = 0; n < n_sim; n++) {						//A loop over all simulations in the file - Start SIMULATIONS LOOP

  /*	In this part for each particle in the simulation
	that fulfil certain criteria
	distance from the canter of the cacade is
	calculataed and placed in an array.

  */
  
  sim->GetEntry(n);								//For each simulation
  int len0 = l_pz->GetLen();							//Setting number of all particles of the shower

  double r_lst_temp[len0];							//An array of all particles distances from the canter of the cacade - r
  int ID_lst_temp[len0];							//An array of all particles identificators - ID

  int N_mi = 0;									//Number of particles that filfil given criteria
  double x, y;									//x and y coordinates of the paricle
  double p_x, p_y, p_z;								//Values of momenta of the particle in each direction
  double p_mi;									//Total momenta of a particle that filfil given criteria


  for(int i = 0; i < len0; i++)  {						//A loop over all particles of the simulations  - Start PARICLES LOOP
    r_lst_temp[i] = 0;
    ID_lst_temp[i] = 0;
    
    int ID = l_id->GetValue(i);
    //p_mi = l_pz->GetValue(i);							//Approximating the momenta of the particle by its z component

    p_x = l_px->GetValue(i);
    p_y = l_py->GetValue(i);    
    p_z = l_pz->GetValue(i);       
    p_mi = sqrt(pow(p_x, 2) + pow(p_y, 2) + pow(p_z, 2));			//Calculating the momenta of the particle

    int is_ID = 0;

    for(int id = 0; id < id_len; id++) {
      if(ID == ID_criterium[id]) is_ID = 1;
    }
    
    if((p_mi >= p_part_min) && (p_mi <= p_part_max) && (is_ID != 0)) {		//Checking the criteria 
									
      x = l_x->GetValue(i);
      y = l_y->GetValue(i);

      r_lst_temp[N_mi] = sqrt((x*x) + (y*y));					//Filling an array of r
      ID_lst_temp[N_mi] = ID;							//Filling an array of ID
      N_mi++;									//Counting number of particles that filfil given criteria
    }
  }										//End PARTICLES LOOP

  N_part_temp_lst[n] = N_mi;							//Filling an array of particles in each simulation with found number of particles
  if(N_part_temp_lst[n] > N_part_max) N_part_max = N_part_temp_lst[n];		//Setting it as maximum if it is greater than previus one - temporary
  if(N_part_max > N_part_maximum) N_part_maximum = N_part_max;			//Setting it as maximum if it is greater than previus one - global
 
  //-------------------------------------------------------------

  /*	This part of code just copies previous arrays
	but it skips the zeros.
  */
	

  double r_lst[N_mi];								//New array of r without the zeros
  int ID_lst[N_mi];								//New array of ID's without the zeros

  int i_mi = 0;

  for(int i = 0; i < len0; i++)  {
    if(r_lst_temp[i] != 0) {
      r_lst[i_mi] = r_lst_temp[i];						//Filling a new array or r
      ID_lst[i_mi] = ID_lst_temp[i];						//Filling a new array or ID
      i_mi++;
    }
  }

  //-------------------------------------------------------------

  sort(r_lst, r_lst + N_mi, greater<double>());					//Sorting previously defined and filled arrays r and ID (in descending order)

  //-------------------------------------------------------------

  /*	Here the R_prc is calculated.
  */

  double maxR[100]; 
 
  for(int j = 0; j < 100; j++) {
    int sum = 0;
    double n_prc = (j+1)*N_mi/100.0;
    maxR[j] = 0;

    for(int i = N_mi-1; i >= 0; i--) {

      if (sum >= n_prc) break;
      sum = sum + 1;
      maxR[j] = r_lst[i];
    }
  }

  R_prc_temp_lst[n] = maxR[prc];						//Filling an array of R_prc
  if(R_prc_temp_lst[n] > R_prc_max) R_prc_max = R_prc_temp_lst[n];		//Setting it as maximum if it is greater than previus one - temporary

  //-------------------------------------------------------------

  /*	In this part of the script the temporary
	histograms, that are used later are defined.
  */

  string title_rho_avr = " #rho_{avr}(r) distribution for E_{casc} = ";		//Setting the title
  title_rho_avr.append(to_string(E_lst[nf]));					//For each energy
  //title_rho_avr.append(" GeV ; r [cm] ; #rho [particles/cm^{2}]");		//Setting title in GeV
  title_rho_avr.append(" TeV ; r [cm] ; #rho [particles/cm^{2}]");		//Setting title in TeV

  const char* rho_avr_title = title_rho_avr.c_str();

  TH1 * h_nr_temp = new TH1F("h_nr_temp", " ; ; ", bins_R_rho, 0.0, maxR[99]);
  TH1 * h_rho_temp = new TH1F("h_rho_temp", rho_avr_title , bins_R_rho, 0.0, maxR[99]);


  //-------------------------------------------------------------

  /*	Here the script find R_rho and
	it fills the array of R_rho with found value.
  */

  for(int ir = 0; ir < N_mi; ir++) {
    h_nr_temp->Fill(r_lst[ir]);							//Filling the histogram with r of each particle
  }

  for(int b = 1; b <= bins_R_rho; b++) {

    double rho_temp = h_nr_temp->GetBinContent(b)/((2*M_PI)*h_nr_temp->GetBinWidth(b)*h_nr_temp->GetBinCenter(b));	
    h_rho_temp->SetBinContent(b, rho_temp);					//Filling the histogram with calculated particles density
    			
    if(rho_temp <= rho_bg) {							//Checking if the particles density is greater than the BG level
										//If it finally drops such low r is saved as R_rho in an array
      R_rho_temp_lst[n] = h_rho_temp->GetBinCenter(b);				//and replaced on each iteration than the loop stops
      break;						

    }
  }

  if(R_rho_temp_lst[n] > R_rho_max) R_rho_max = R_rho_temp_lst[n]; 		//Setting it as maximum if it is greater than previus one - temporary


  h_nr_temp->Delete();								//Delating the histograms
  h_rho_temp->Delete();

  }										//End SIMULATIONS LOOP

  //-------------------------------------------------------------

  /*	Here the histograms of n(N_part) are created,
	Gaussian distributions are fitted and than
	they are saved.
  */

  string title_N_part = " n (cascades) in N_{part} distribution for E_{casc} = ";
  title_N_part.append(to_string(E_lst[nf]));
  //title_N_mi.append(" GeV ; N_{part} ; n (cascades)");
  title_N_part.append(" TeV ; N_{part} ; n (cascades)");

  const char* N_part_title = title_N_part.c_str();

  TH1 * h_N_part = new TH1F("h_N_part", N_part_title , bins_N_part, 0, 1.1*N_part_max);

  for(int n = 0; n < n_sim; n++) {
    h_N_part->Fill(N_part_temp_lst[n]);						//Filling tha histogram with N_part
  }

  TCanvas C("C");  

  cout<<endl;
  TFitResultPtr result = h_N_part->Fit("gaus", "S");				//Fitting the Gaussian distribution

  cN_part_lst[nf] = result->Parameter(0);					//Saving its parameter C
  meanN_part_lst[nf] = result->Parameter(1);					//Saving its parameter Mean
  sigmaN_part_lst[nf] = result->Parameter(2);					//Saving its parameter Sigma

  string path_N_part = N_part_dir;
  path_N_part.append("/n(N_part_");
  path_N_part.append(to_string(E_lst[nf]));
  path_N_part.append(").png");

  const char* N_part_path = path_N_part.c_str();

  cout<<endl;
  h_N_part->Draw();
  C.SaveAs(N_part_path);							//Saving it in the created directory
  h_N_part->Delete();								//Delating the histogram

  //-------------------------------------------------------------


  /*	Here the histograms of n(R_prc) are created,
	Gaussian distributions are fitted and than
	they are saved.
  */

  string title_R_prc = " n (cascades) in R_{prc}_{avr} distribution for E_{casc} = ";
  title_R_prc.append(to_string(E_lst[nf]));
  //title_R_prc.append(" GeV ; R_{prc} [cm]; n (cascades)");
  title_R_prc.append(" TeV ; R_{prc} [cm]; n (cascades)");

  const char* R_prc_title = title_R_prc.c_str();

  TH1 * h_R_prc = new TH1F("h_r", R_prc_title, bins_R_prc, 0, 1.1*R_prc_max);

  for(int n = 0; n < n_sim; n++) {
    h_R_prc->Fill(R_prc_temp_lst[n]);						//Filling tha histogram with R_prc
  }


  TCanvas B("B");  

  cout<<endl;
  TFitResultPtr resultR_prc = h_R_prc->Fit("gaus", "S");  			//Fitting the Gaussian distribution

  cR_prc_lst[nf] = resultR_prc->Parameter(0);					//Saving its parameter C
  meanR_prc_lst[nf] = resultR_prc->Parameter(1);				//Saving its parameter Mean
  sigmaR_prc_lst[nf] = resultR_prc->Parameter(2);				//Saving its parameter Sigma

  if(meanR_prc_lst[nf] > R_prc_maximum) R_prc_maximum = meanR_prc_lst[nf]; 	//Setting it as maximum if it is greater than previus one - global

  string path_R_prc = R_prc_dir;
  path_R_prc.append("/n(R_prc_");
  path_R_prc.append(to_string(E_lst[nf]));
  path_R_prc.append(").png");

  const char* R_prc_path = path_R_prc.c_str();

  cout<<endl;
  h_R_prc->Draw();
  B.SaveAs(R_prc_path);								//Saving it in the created directory
  h_R_prc->Delete();								//Delating the histogram

  //-------------------------------------------------------------

  /*	Here the histograms of n(R_rho) are created,
	Gaussian distributions are fitted and than
	they are saved.
  */

  string title_R_rho = " n (cascades) in R_{rho}_{avr} distribution for E_{casc} = ";
  title_R_rho.append(to_string(E_lst[nf]));
  //title_R_rho.append(" GeV ; R_{rho} [cm]; n (cascades)");
  title_R_rho.append(" TeV ; R_{rho} [cm]; n (cascades)");

  const char* R_rho_title = title_R_rho.c_str();

  TH1 * h_R_rho = new TH1F("h_R_rho", R_rho_title, bins_R_rho, 0, 1.1*R_rho_max);

  for(int n = 0; n < n_sim; n++) {
    h_R_rho->Fill(R_rho_temp_lst[n]);						//Filling tha histogram with R_prc
  }

  TCanvas H("H");

  cout<<endl;
  TFitResultPtr resultR_rho = h_R_rho->Fit("gaus", "S");			//Fitting the Gaussian distribution  

  cR_rho_lst[nf] = resultR_rho->Parameter(0);					//Saving its parameter C
  meanR_rho_lst[nf] = resultR_rho->Parameter(1);				//Saving its parameter Mean
  if(meanR_rho_lst[nf] <= 0) meanR_rho_lst[nf] = h_R_rho->GetMean();	
  sigmaR_rho_lst[nf] = resultR_rho->Parameter(2);				//Saving its parameter Sigma

  if(meanR_rho_lst[nf] > R_rho_maximum) R_rho_maximum = meanR_rho_lst[nf]; 	//Setting it as maximum if it is greater than previus one - global

  string path_R_rho = R_rho_dir;
  path_R_rho.append("/n(R_rho_");
  path_R_rho.append(to_string(E_lst[nf]));
  path_R_rho.append(").png");

  const char* R_rho_path = path_R_rho.c_str();

  cout<<endl;
  h_R_rho->Draw();
  H.SaveAs(R_rho_path);								//Saving it in the created directory
  h_R_rho->Delete();								//Delating the histogram

  //-------------------------------------------------------------

  /*	Here the 2D histograms of n(N_part, R_prc) are created,
	Gaussian distributions are fitted and than
	they are saved.
  */

  string title_RN_prc = " R_{prc} and N_{part} for E_{casc} = ";
  title_RN_prc.append(to_string(E_lst[nf]));
  //title_RN_prc.append(" GeV ; N_{part} ; R_{prc} [cm]; n (cascades)");
  title_RN_prc.append(" TeV ; N_{part} ; R_{prc} [cm]; n (cascades)");

  const char* RN_prc_title = title_RN_prc.c_str();

  TH2 * h_R_N = new TH2F("h_R_N", RN_prc_title , bins_N_part, 0, 1.1*N_part_max, bins_N_part, 0, 1.1*R_prc_max);

  for(int n = 0; n < n_sim; n++) {	
    h_R_N->Fill(N_part_temp_lst[n], R_prc_temp_lst[n]);				//Filling tha histogram with N_part and R_prc
  }

  TCanvas I("I");

  string path_RN_prc = RN_prc_dir;
  path_RN_prc.append("/n(N_part, R_prc)_");
  path_RN_prc.append(to_string(E_lst[nf]));
  path_RN_prc.append(".png");

  const char* RN_prc_path = path_RN_prc.c_str();

  h_R_N->Draw("colz");
  I.SaveAs(RN_prc_path);							//Saving it in the created directory
  h_R_N->Delete();								//Delating the histogram


  nf++;

  }
  }
  }										//End FILES LOOP

  //-------------------------------------------------------------

  if((dir = opendir(sim_dir.c_str())) != NULL) {				//Continue only if the path to directory with simulations is correct			

  //-------------------------------------------------------------

  /*	In this part the functions of meanN_part(E),
	sigmaN_part(E) and CN_part(E) are defined.
  */


  TF1 * f_meanN_part_E = new TF1("f_meanN_part_E" , " [0]*pow(x, [1]) + [2]" , 0, E_lst[n_datas]);
  f_meanN_part_E->SetParameters(15.5473, 0.91507, 0.0322778);

  TF1 * f_sigmaN_part_E = new TF1("f_sigmaN_part_E" , " [0]*pow(x, [1]) + [2]" , 0, E_lst[n_datas]);
  f_sigmaN_part_E->SetParameters(4.71603, 0.859976, 5.07782);

  TF1 * f_cN_part_E = new TF1("f_cN_part_E" , " [0]*pow(x, [1]) + [2]" , 0, E_lst[n_datas]);
  f_cN_part_E->SetParameters(0.00736672, 0.884137, 105.498);

  //-------------------------------------------------------------

  /*	Here the graphs of meanN_part(E), sigmaN_part(E),
	cN_part(E) and sigmaN_part/meanN_part(E) are created.
	Functions defined earlier are fitted
	and than they are saved.
  */

  TCanvas D("D");
  D.SetLogy(); 
  D.SetLogx();   

  TGraph * g_meanN_part = new TGraph(n_datas, E_lst, meanN_part_lst);
  TGraph * g_sigmaN_part = new TGraph(n_datas, E_lst, sigmaN_part_lst);
  TGraph * g_cN_part = new TGraph(n_datas, E_lst, cN_part_lst);


  double N_frac[n_datas];
  for(int nf = 0; nf < n_datas; nf++) {
    N_frac[nf] = sigmaN_part_lst[nf]/meanN_part_lst[nf];			//Calculating ratio sigmaN_part/meanN_part
  }

  TGraph * g_N_part_frac = new TGraph(n_datas, E_lst, N_frac);

  //-------------------------------------------------------------

  /*	In this part functions to fitted and
	their parameters are saved into ".txt" file.
  */

  ofstream param_file;
  param_file.open("N(E)_params.txt");

  cout<< "\n Fitting meanN_part(E)" <<endl;					//Fitting functions to the graphs
  TFitResultPtr result_meanN_part = g_meanN_part->Fit(f_meanN_part_E, "S");
  param_file << result_meanN_part->Parameter(0) << "\n";
  param_file << result_meanN_part->Parameter(1) << "\n";
  param_file << result_meanN_part->Parameter(2) << "\n";

  cout<< "\n Fitting sigma_N_part(E)" <<endl;
  TFitResultPtr result_sigmaN_part = g_sigmaN_part->Fit(f_sigmaN_part_E, "S");
  param_file << result_sigmaN_part->Parameter(0) << "\n";
  param_file << result_sigmaN_part->Parameter(1) << "\n";
  param_file << result_sigmaN_part->Parameter(2) << "\n";
  
  cout<< "\n Fitting c_N_part(E)" <<endl;
  TFitResultPtr result_cN_part = g_cN_part->Fit(f_cN_part_E, "S");
  param_file << result_cN_part->Parameter(0) << "\n";
  param_file << result_cN_part->Parameter(1) << "\n";
  param_file << result_cN_part->Parameter(2) << "\n";

  cout<< "\n Fitting sigma_N_part/meanN_part (E)" <<endl;
  g_N_part_frac->Fit("expo", "S");


  //-------------------------------------------------------------

  /*	Here also maximum and minimum energies
	E_min adn E_max are saved into a file,
	as well as distances R_prc_max, R_rho_max
	and maximum number of parts N_part_max.
  */

  param_file << E_min << "\n";
  param_file << E_max << "\n";

  param_file << R_prc_maximum << "\n";
  param_file << R_rho_maximum << "\n";

  param_file << N_part_maximum << "\n";

  param_file.close();

  //-------------------------------------------------------------

  g_meanN_part->SetTitle("#LT N_{part} #GT(E)");				//Setting titles of plot, axises etc.
  g_sigmaN_part->SetTitle("#sigma N_{part}(E)");
  g_cN_part->SetTitle("C(N_{part})(E)");
  g_N_part_frac->SetTitle("#sigma N_{part}/#LT N_{part} #GT (E)");

  //g_meanN_part->GetXaxis()->SetTitle("E_{casc} [GeV]");
  //g_sigmaN_part->GetXaxis()->SetTitle("E_{casc} [GeV]");
  //g_cN_part->GetXaxis()->SetTitle("E_{casc} [GeV]");
  //g_N_part_frac->GetXaxis()->SetTitle("E_{casc [GeV]
  g_meanN_part->GetXaxis()->SetTitle("E_{casc} [TeV]");
  g_sigmaN_part->GetXaxis()->SetTitle("E_{casc} [TeV]");
  g_cN_part->GetXaxis()->SetTitle("E_{casc} [TeV]");
  g_N_part_frac->GetXaxis()->SetTitle("E_{casc} [TeV]");

  g_meanN_part->GetYaxis()->SetTitle("#LT N_{part} #GT");
  g_sigmaN_part->GetYaxis()->SetTitle("#sigma N_{part}");
  g_cN_part->GetYaxis()->SetTitle("C(N_{part})");
  g_N_part_frac->GetYaxis()->SetTitle("#sigma N_{part}/#LT N_{part} #GT");

  g_meanN_part->SetMarkerStyle(8);
  g_sigmaN_part->SetMarkerStyle(8);
  g_cN_part->SetMarkerStyle(8);
  g_N_part_frac->SetMarkerStyle(8);

  g_meanN_part->SetMarkerSize(1);
  g_sigmaN_part->SetMarkerSize(1);
  g_cN_part->SetMarkerSize(1);
  g_N_part_frac->SetMarkerSize(1);


  string path_N_part_E = N_part_E_dir;						//Saving plots
  cout<<endl;

  string path_meanN_part_E = path_N_part_E;
  path_meanN_part_E.append("/meanN_part(E).png");
  const char* meanN_part_E_path = path_meanN_part_E.c_str();

  g_meanN_part->Draw("AP");							
  D.SaveAs(meanN_part_E_path);
  g_meanN_part->Delete();


  string path_sigmaN_part_E = path_N_part_E;
  path_sigmaN_part_E.append("/sigmaN_part(E).png");
  const char* sigmaN_part_E_path = path_sigmaN_part_E.c_str();

  g_sigmaN_part->Draw("AP");
  D.SaveAs(sigmaN_part_E_path);
  g_sigmaN_part->Delete();


  string path_cN_part_E = path_N_part_E;
  path_cN_part_E.append("/cN_part(E).png");
  const char* cN_part_E_path = path_cN_part_E.c_str();

  g_cN_part->Draw("AP");
  D.SaveAs(cN_part_E_path);
  g_cN_part->Delete();


  string path_N_sigma_frac_E = path_N_part_E;
  path_N_sigma_frac_E.append("/N_sigma_frac(E).png");
  const char* N_sigma_frac_E_path = path_N_sigma_frac_E.c_str();

  g_N_part_frac->Draw("AP");
  D.SaveAs(N_sigma_frac_E_path);
  g_N_part_frac->Delete();						


  //-------------------------------------------------------------


  /*	In this part the functions of meanR_prc(E),
	sigmaR_prc(E) are defined.
  */

  TF1 * f_meanR_prc_E = new TF1("f_meanR_prc_E" , " [0]*pow(x, [1]) + [2]" , 0, E_lst[n_datas]);
  f_meanR_prc_E->SetParameters(-3.32045e+07, 0.000957474, 4.77929e+06);

  TF1 * f_sigmaR_prc_E = new TF1("f_sigmaR_prc_E" , " [0]*pow(x, [1]) + [2]" , 0, E_lst[n_datas]);
  f_sigmaR_prc_E->SetParameters(-3.32045e+07, 0.000137065, 3.32486e+07);

  //-------------------------------------------------------------

  /*	Here the graphs of meanR_prc(E), sigmaR_prc(E),
	cR_prc(E) are created.
	Functions defined earlier are fitted
	and than they are saved.
  */

  TCanvas F("F");
  F.SetLogy(); 
  F.SetLogx();   


  TGraph * g_meanR_prc = new TGraph(n_datas, E_lst, meanR_prc_lst);
  TGraph * g_sigmaR_prc = new TGraph(n_datas, E_lst, sigmaR_prc_lst);
  TGraph * g_cR_prc = new TGraph(n_datas, E_lst, cR_prc_lst);

  cout<< "\n Fitting meanR_prc(E)" <<endl;					//Fitting functions to the graphs
  g_meanR_prc->Fit("expo");
  cout<< "\n Fitting sigma_R_prc(E)" <<endl;
  g_sigmaR_prc->Fit(f_sigmaR_prc_E);
  cout<< "\n Fitting cR_prc(E)" <<endl;
  g_cR_prc->Fit("pol1");


  g_meanR_prc->SetTitle("#LT R_{prc} #GT(E)");					//Setting titles of plot, axises etc.
  g_sigmaR_prc->SetTitle("#sigmaR_{prc}(E)");
  g_cR_prc->SetTitle("C(R_{prc})(E)");

  //g_meanR_prc->GetXaxis()->SetTitle("E_{casc} [GeV]");
  //g_sigmaR_prc->GetXaxis()->SetTitle("E_{casc} [GeV]");
  //g_cR_prc->GetXaxis()->SetTitle("E_{casc} [GeV]");
  g_meanR_prc->GetXaxis()->SetTitle("E_{casc} [TeV]");
  g_sigmaR_prc->GetXaxis()->SetTitle("E_{casc} [TeV]");
  g_cR_prc->GetXaxis()->SetTitle("E_{casc} [TeV]");

  g_meanR_prc->GetYaxis()->SetTitle("#LT R_{prc} #GT [cm]");
  g_sigmaR_prc->GetYaxis()->SetTitle("#sigma R_{prc} [cm]");
  g_cR_prc->GetYaxis()->SetTitle("C(R_{prc}) [cm]");

  g_meanR_prc->SetMarkerStyle(8);
  g_sigmaR_prc->SetMarkerStyle(8);
  g_cR_prc->SetMarkerStyle(8);

  g_meanR_prc->SetMarkerSize(1);
  g_sigmaR_prc->SetMarkerSize(1);
  g_cR_prc->SetMarkerSize(1);


  string path_R_prc_E = R_prc_E_dir;						//Saving plots
  cout<<endl;

  string path_meanR_prc_E = path_R_prc_E;
  path_meanR_prc_E.append("/meanR_prc(E).png");
  const char* meanR_prc_E_path = path_meanR_prc_E.c_str();

  g_meanR_prc->Draw("AP");							
  F.SaveAs(meanR_prc_E_path);
  g_meanR_prc->Delete();


  string path_sigmaR_prc_E = path_R_prc_E;
  path_sigmaR_prc_E.append("/sigmaR_prc(E).png");
  const char* sigmaR_prc_E_path = path_sigmaR_prc_E.c_str();

  g_sigmaR_prc->Draw("AP");
  F.SaveAs(sigmaR_prc_E_path);
  g_sigmaR_prc->Delete();


  string path_cR_prc_E = path_R_prc_E;
  path_cR_prc_E.append("/cR_prc(E).png");
  const char* cR_prc_E_path = path_cR_prc_E.c_str();

  g_cR_prc->Draw("AP");
  F.SaveAs(cR_prc_E_path);
  g_cR_prc->Delete();


  //-------------------------------------------------------------

  /*	Here the graphs of meanR_rho(E), sigmaR_rho(E),
	cR_rho(E) are created.
	Functions defined earlier are fitted
	and than they are saved.
  */

  TCanvas G("G"); 
  G.SetLogy(); 
  G.SetLogx();  

  TGraph * g_meanR_rho = new TGraph(n_datas, E_lst, meanR_rho_lst);
  TGraph * g_sigmaR_rho = new TGraph(n_datas, E_lst, sigmaR_rho_lst);
  TGraph * g_cR_rho = new TGraph(n_datas, E_lst, cR_rho_lst);

  cout<< "\n Fitting meanR_rho(E)" <<endl;					//Fitting functions to the graphs
  g_meanR_rho->Fit("expo");
  cout<< "\n Fitting sigma_R_rho(E)" <<endl;
  g_sigmaR_rho->Fit("expo");
  cout<< "\n Fitting cR_rho(E)" <<endl;
  g_cR_rho->Fit("pol1");



  g_meanR_rho->SetTitle("#LT R_{rho} #GT(E)");					//Setting titles of plot, axises etc.
  g_sigmaR_rho->SetTitle("#sigma R_{rho}(E)");
  g_cR_rho->SetTitle("C(R_{rho})(E)");

  //g_meanR_rho->GetXaxis()->SetTitle("E_{casc} [GeV]");
  //g_sigmaR_rho->GetXaxis()->SetTitle("E_{casc} [GeV]");
  //g_cR_rho->GetXaxis()->SetTitle("E_{casc} [GeV]");
  g_meanR_rho->GetXaxis()->SetTitle("E_{casc} [TeV]");
  g_sigmaR_rho->GetXaxis()->SetTitle("E_{casc} [TeV]");
  g_cR_rho->GetXaxis()->SetTitle("E_{casc} [TeV]");

  g_meanR_rho->GetYaxis()->SetTitle("#LT R_{rho} #GT [cm]");
  g_sigmaR_rho->GetYaxis()->SetTitle("#sigma R_{rho} [cm]");
  g_cR_rho->GetYaxis()->SetTitle("C(R_{rho}) [cm]");

  g_meanR_rho->SetMarkerStyle(8);
  g_sigmaR_rho->SetMarkerStyle(8);
  g_cR_rho->SetMarkerStyle(8);

  g_meanR_rho->SetMarkerSize(1);
  g_sigmaR_rho->SetMarkerSize(1);
  g_cR_rho->SetMarkerSize(1);


  string path_R_rho_E = R_rho_E_dir;						//Saving plots
  cout<<endl;

  string path_meanR_rho_E = path_R_rho_E;
  path_meanR_rho_E.append("/meanR_rho(E).png");
  const char* meanR_rho_E_path = path_meanR_rho_E.c_str();

  g_meanR_rho->Draw("AP");							
  G.SaveAs(meanR_rho_E_path);
  g_meanR_rho->Delete();


  string path_sigmaR_rho_E = path_R_rho_E;
  path_sigmaR_rho_E.append("/sigmaR_rho(E).png");
  const char* sigmaR_rho_E_path = path_sigmaR_rho_E.c_str();

  g_sigmaR_rho->Draw("AP");
  G.SaveAs(sigmaR_rho_E_path);
  g_sigmaR_rho->Delete();


  string path_cR_rho_E = path_R_rho_E;
  path_cR_rho_E.append("/cR_rho(E).png");
  const char* cR_rho_E_path = path_cR_rho_E.c_str();

  g_cR_rho->Draw("AP");
  G.SaveAs(cR_rho_E_path);
  g_cR_rho->Delete();

  //-------------------------------------------------------------

  }  

  //-------------------------------------------------------------


  auto stop = high_resolution_clock::now(); 
  auto duration = duration_cast<microseconds>(stop - start); 			//Stoping counting time of computations

  cout<<"\n" <<endl;
  cout<< "Time of the computations: " << duration.count() << " * 10^{-6} s = " <<  duration.count()/pow(10, 6) << " s = "  << duration.count()/(pow(10, 6)*60) << " min" <<endl;

  //-------------------------------------------------------------

}

//-------------------------------------------------------------

/*	The main C++ program is here and
	it just executes previously written macro.
*/

int main() {
  NR_E();
  return 0;
}

//-------------------------------------------------------------

/* Used functions definitions */

string get_current_dir() {							//Function returning path to current directory

   char buff[FILENAME_MAX]; 							//Create string buffer to hold path

   GetCurrentDir( buff, FILENAME_MAX );
   string current_working_dir(buff);

   return current_working_dir;
}



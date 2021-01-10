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
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TFitResult.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TLeafI.h"
#include "TFile.h"

using namespace std::chrono;
using namespace std;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/* 	This macro is responsible for making histograms,
	graphs and saving them.

	The meaning of each starting parameters of macro:

	- ID_criterium - list of particles IDs,

	- sim_dir - path to directory with files containing CORSIKA simulations,
	- plots_dir - path to directory where all plots are saved (default = current directory),

	- norm_file - name of file for parametrization (default = last file),

	- bins_rho_N_part - number of bins in histograms of number of particles N_part (default = 40),
	- bins_rho_r - number of bins in histogram of particles density Rho(r) (default = 100),

	- r_number - number of points on axis of distance from cascade center r (default = 10)

	- p_part_min - minimum momenta of wanted particles [GeV] (default 1),
	- p_part_max - maximum momenta of wanted particles [GeV] (default pow(10, 10) ~ inf),

	- n_datas - number of simulations (default = 18).
*/


void Rho_rNE_macro(int ID_criterium[], string sim_dir = "/home/jerzy/CREDO/Analiza/More_simulations", string plots_dir = get_current_dir_name(), string norm_file = "", int bins_rho_N_part = 40, int bins_rho_r = 40, int r_number = 10, double p_part_min = 1.0, double p_part_max = pow(10, 10), int n_datas = 18)
{

  //-------------------------------------------------------------

  auto start = high_resolution_clock::now(); 					//Starting counting time of computations

  //------------------------------------------------------------- 

  /*	In this part, previously generated files
	are read and several parameters are extracted
	from them.
  */

  double meanN_part_params[3], sigmaN_part_params[3]; 
  double E_min, E_max;
  double r_prc_max, r_rho_max;
  double N_part_min, N_part_max;  

  //-------------------------

  ifstream mN_param_fit_file("meanN(E)_fit_params.txt");
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

  ifstream sN_param_fit_file("sigmaN(E)_fit_params.txt");
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

  ifstream param_file("N(E)_params.txt");
  string line;

  if(param_file.is_open()) {

    getline(param_file, line);
    E_min = stod(line) ;				//Minimum energy E on the graph [TeV]

    getline(param_file, line);
    E_max = stod(line) ;				//Maximum energy E on the graph [TeV]


    getline(param_file, line);
    r_prc_max = stod(line) ;				//Maximum distance r_prc (from the center of the cascade) on the graph

    getline(param_file, line);
    r_rho_max = stod(line) ;				//Maximum distance r_rho (from the center of the cascade) on the graph (NOT USED)


    N_part_min = 1;					//Minimum number of partciles N_part on the graph

    getline(param_file, line);
    N_part_max = stod(line) ;				//Maximum number of partciles N_part on the graph
    
  }

  else cout<< "ERROR: COULD NOT FIND THE N(E) PARAMETERS FILE" <<endl;

  param_file.close();

  //------------------------------------------------------------- 

  /*	In this part, functions of
	mean_N_part(E), sigma_N_part(E), c_N_part(E)
	are defined.
  */

  TF1 * f_meanN_part_E = new TF1("f_meanN_part_E" , " [0]*pow(x, [1]) + [2] ");
  f_meanN_part_E->SetParameters(meanN_part_params[0], meanN_part_params[1], meanN_part_params[2]);		//Function found earlier - meanN_part(E)

  TF1 * f_sigmaN_part_E = new TF1("f_sigmaN_part_E" , " [0]*pow(x, [1]) + [2] ");
  f_sigmaN_part_E->SetParameters(sigmaN_part_params[0], sigmaN_part_params[1], sigmaN_part_params[2]);		//Function found earlier - sigmaN_part(E)

  //------------------------------------------------------------- 

  /*	Here distance points for fitting
	Rho(r, E) are set.
  */

  long double R_fit_lst[r_number];
  long double r_fit_temp = r_prc_max/(pow(2, r_number - 1));

  for(int rp = 0; rp < r_number; rp++) {

    //R_fit_lst[rp] = r_fit_temp;						//Equally distributed for log scale
    //r_fit_temp = r_fit_temp * 2;

    R_fit_lst[rp] = ((double(rp+1)/double(r_number))) * pow(10, (double(rp+1)/double(r_number)) - 1) * r_prc_max;
										//Almost equally distributed for log scale (no repetitions in bins for small r_number)
    if(rp == r_number - 1) R_fit_lst[rp] = R_fit_lst[rp] - 100.0;		//Modification preventing from getting out of histogram range

  } 

  //------------------------------------------------------------- 

  /*	This part creates log bining for
	Rho(r) histograms.
  */

  double Rr_bin_lst[bins_rho_r];
  Rr_bin_lst[0] = 0.0;
  double Rr_bin_temp = 100.0;
  double bin_w_scale = pow(10.0, log10(r_prc_max/Rr_bin_temp)/double(bins_rho_r - 1));
  
  cout<< Rr_bin_lst[0] <<endl;
  for(int rb = 1; rb < bins_rho_r - 1; rb++) {

    Rr_bin_lst[rb] = Rr_bin_temp;						//Equally distributed for log scale
    Rr_bin_temp = Rr_bin_temp * bin_w_scale;
    //cout<< Rr_bin_lst[rb] <<endl;

  } 

  Rr_bin_lst[bins_rho_r - 1] = r_prc_max - 1.0;					//Modification preventing from getting out of histogram range

  cout<< Rr_bin_lst[bins_rho_r - 1] <<endl;

  const double *Rr_bins = Rr_bin_lst;

  //-------------------------------------------------------------

  /*	Here are only definitions of some
	arrays and variables used later.
  */

  double Rho_NE_p0_lst[n_datas], Rho_NE_p1_lst[n_datas], Rho_NE_p2_lst[n_datas];			//An arrays of parameters of Rho(r, E) function
  double Rho_rN_p0_lst[n_datas], Rho_rN_p1_lst[n_datas], Rho_rN_p2_lst[n_datas];			//An arrays of parameters of Rho(r, N) function

  double N_part_lst[bins_rho_N_part], N_part_error[bins_rho_N_part];					//An arrays of N_part/N_avr ratio
  double n_N_part_counter[bins_rho_N_part], sum_n_N_part[bins_rho_N_part];				//An arrays of n(N_part) helping parameters

  double E_lst[n_datas];										//An array of energies of the simulations E
  double E_error[n_datas];										//An array of errors for E
  int n_sim_lst[n_datas];										//An array of number of simulations for each file

  double Rho_rE_min = 1.0;										//Minimum Rho(r, E) for the histogram
  double Rho_rE_max = 0.0;										//Maximum Rho(r, E) for the histogram


  long double** Rho_E_lst = 0;				//An array of particles density distribution Rho(r, E)
  Rho_E_lst = new long double*[n_datas];

  for(int i = 0; i < n_datas; i++) {
    Rho_E_lst[i] = new long double[r_number];

    for(int j = 0; j < r_number; j++) {
      Rho_E_lst[i][j] = 0;
    }
  }

  long double** Rho_E_error = 0;			//An array of errors of particles density distribution Rho(r, E)
  Rho_E_error = new long double*[n_datas];

  for(int i = 0; i < n_datas; i++) {
    Rho_E_error[i] = new long double[r_number];

    for(int j = 0; j < r_number; j++) {
      Rho_E_error[i][j] = 0;
    }
  }




  long double** Rho_N_lst = 0;				//An array of particles density distribution Rho(r, N)
  Rho_N_lst = new long double*[bins_rho_N_part];

  for(int i = 0; i < bins_rho_N_part; i++) {
    Rho_N_lst[i] = new long double[r_number];

    for(int j = 0; j < r_number; j++) {
      Rho_N_lst[i][j] = 0;
    }
  }

  long double** Rho_N_error = 0;			//An array of errors of particles density distribution Rho(r, N)
  Rho_N_error = new long double*[bins_rho_N_part];

  for(int i = 0; i < bins_rho_N_part; i++) {
    Rho_N_error[i] = new long double[r_number];

    for(int j = 0; j < r_number; j++) {
      Rho_N_error[i][j] = 0;
    }
  }

  //-------------------------------------------------------------  

  /*	Here are defined some directories
	in which plots are stored.
  */

  string Rho_rN_dir_name = "/Rho_(r,N)_2D";					//Naming directories for histograms
  string Rho_rN_dir = plots_dir;
  Rho_rN_dir.append(Rho_rN_dir_name);

  string Rho_rN_norm_dir_name = "/Rho_(r,N)_2D_norm";				//Naming directories for histograms
  string Rho_rN_norm_dir = plots_dir;
  Rho_rN_norm_dir.append(Rho_rN_norm_dir_name);


  string Rho_avr_r_dir_name = "/Rho_avr(r)";					//Naming directories for graphs
  string Rho_avr_r_dir = plots_dir;
  Rho_avr_r_dir.append(Rho_avr_r_dir_name);

  string Rho_rE_norm_dir_name = "/Rho_(r,E)_1D_norm";				//Naming directories for graphs
  string Rho_rE_norm_dir = plots_dir;
  Rho_rE_norm_dir.append(Rho_rE_norm_dir_name);

  string Rho_rN_1D_norm_dir_name = "/Rho_(r,N)_1D_norm";			//Naming directories for graphs
  string Rho_rN_1D_norm_dir = plots_dir;
  Rho_rN_1D_norm_dir.append(Rho_rN_1D_norm_dir_name);



  int check_Rho_avr_r = mkdir(Rho_avr_r_dir.c_str(), 0777);			//Creating directories for histograms
  int check_Rho_rN = mkdir(Rho_rN_dir.c_str(), 0777);				//Creating directories for histograms
  int check_Rho_norm_rN = mkdir(Rho_rN_norm_dir.c_str(), 0777);			//Creating directories for histograms

  int check_Rho_norm_rE = mkdir(Rho_rE_norm_dir.c_str(), 0777);			//Creating directories for graphs
  int check_Rho_norm_rE_1D = mkdir(Rho_rN_1D_norm_dir.c_str(), 0777);		//Creating directories for graphs

  //-------------------------------------------------------------  

  /*	This part of code defines
	a 2D histogram of particles density
	Rho(r, E).
  */

  TH2 * h_Rho_rE_2D = new TH2F("h_rho_rE", " #rho(r, E) distribution; r [cm]; File number ; #rho [cm^{-2}] ", bins_rho_r - 1, Rr_bins, n_datas, 0, n_datas);

  //-------------------------------------------------------------

  int nf = 0;									//Files iterator
  int n_norm = n_datas - 1;							//Number of file for normalization
  //int id_len = sizeof(ID_criterium)/sizeof(ID_criterium[0]);			//Length of IDs array
  int id_len = 12;								//Length of IDs array

  //-------------------------------------------------------------

  DIR *dir;
  struct dirent *ent;

  if((dir = opendir(sim_dir.c_str())) != NULL) {				
  while((ent = readdir(dir)) != NULL) {						//A loop over all files in directory - Start FILES LOOP
  if(opendir(ent->d_name) == NULL) {

  /*	This part of code reads each file,
	finds energy of primary particle for
	simulated cascades and number of simulations
	in the file.

	There are also some definitions of
	variables and lists used later.
  */

  cout<< "\nFile name: " <<endl;
  cout<< ent->d_name <<endl;

  string file_sim = sim_dir;
  file_sim.append("/");						
  file_sim.append(ent->d_name);							//Adding name of the file into directory path

  if(file_sim.find(".root") > file_sim.size()) continue;			//Checking if its a root file

  cout<< "\nOpening file nr " << nf+1 << " : " <<endl;
  cout<< file_sim <<endl;

  TFile *file = new TFile(file_sim.c_str(), "READ");				//Opening each file
  TTree *sim = (TTree*) file->Get("sim");
  
  TLeafI* l_id = (TLeafI*)sim->GetLeaf("particle..ParticleID");
  TLeaf* l_x = (TLeaf*)sim->GetLeaf("particle..x");
  TLeaf* l_y = (TLeaf*)sim->GetLeaf("particle..y");
  TLeaf* l_pz = (TLeaf*)sim->GetLeaf("particle..Pz");
  TLeaf* l_px = (TLeaf*)sim->GetLeaf("particle..Px");
  TLeaf* l_py = (TLeaf*)sim->GetLeaf("particle..Py");

  TLeaf* E_sim = (TLeaf*)sim->GetLeaf("shower.Energy");

  if(ent->d_name == norm_file) n_norm = nf;					//Checking if this is a file chosen for normalisation 

  //-------------------------------------------------------------  

  /*	Here some arrays are filled with zeros.
  */

  for(int n = 0; n < bins_rho_N_part; n++) {

    N_part_lst[n] = N_part_error[n] = 0;
    n_N_part_counter[n] = sum_n_N_part[n] = 0;


    for(int rf = 0; rf < r_number; rf++) {
      Rho_N_lst[n][rf] = Rho_N_error[n][rf] = 0;
    }
  }

  //-------------------------------------------------------------  

  /*	This part of code reads energy of
	the cascades in each file
	(assuming that all showers within the file
	have the same energy). It also atribute
	an error to the energy.
  */

  sim->GetEntry(0);
  //E_lst[nf] = E_sim->GetValue(0);						//Filling an array with energy in GeV
  E_lst[nf] = E_sim->GetValue(0)/1000.0;					//Filling an array with energy in TeV
  E_error[nf] = 0.05*E_lst[nf];							//Filling an array with energy error in TeV

  //cout<< "\nEnergy of primary particle: " << E_lst[nf] << " GeV" <<endl;	//Printing energy in GeV 
  cout<< "\nEnergy of primary particle: " << E_lst[nf] << " TeV" <<endl;	//Printing energy in TeV

  //-------------------------------------------------------------

  /*	Here, simulations are counted.
  */

  int n_sim = 0;
  while(sim->GetEntry(n_sim) != 0) {
    n_sim++;
  }										//Counting number of simulations in the file		

  n_sim_lst[nf] = n_sim;

  cout<< "Number of simulations in the file: " << n_sim <<endl;			//Printing number of simulations in the file
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

  double Rho_rN_min = 1.0;							//Minimum Rho(r, N) for the histogram
  double Rho_rN_max = 0.0;							//Maximum Rho(r, N) for the histogram


  int N_part_temp_lst[n_sim];							//An array of N_part for each simulation in the file
  int N_part_max = 0;								//Maximum N_part - temporary


  //-------------------------------------------------------------

  /*	In this part of code, temporary histogram
	needed to calculate particles density Rho(r)
	over distance r is created.
  */

  TH1 * h_nr_temp = new TH1F(" ", " ; ; ", bins_rho_r - 1, Rr_bins);

  //-------------------------------------------------------------

  /*	Here a temporary 2D histogram needed
	to create other histograms is declared.
  */

  double N_rN_max = f_meanN_part_E->Eval(E_lst[nf]) + 4*f_sigmaN_part_E->Eval(E_lst[nf]);			//Setting ranges using previousle found functions
  double N_rN_min = f_meanN_part_E->Eval(E_lst[nf]) - 4*f_sigmaN_part_E->Eval(E_lst[nf]);
  if(N_rN_min <= 0) N_rN_min = N_part_min;

  TH2 * h_nr_temp_2D = new TH2F(" ", " ; ; ; ", bins_rho_r - 1, Rr_bins, bins_rho_N_part - 1, N_rN_min, N_rN_max);

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

  int N_mi = 0;									//Number of particles that filfil given criteria
  double x, y;									//x and y coordinates of the paricle
  double p_x, p_y, p_z;								//Values of momenta of the particle in each direction
  double p_mi;									//Total momenta of a particle that filfil given criteria


  for(int i = 0; i < len0; i++)  {						//A loop over all particles of the simulations  - Start PARICLES LOOP
    r_lst_temp[i] = 0;
    
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
      N_mi++;									//Counting number of particles that filfil given criteria
    }
  }										//End PARTICLES LOOP

  N_part_temp_lst[n] = N_mi;							//Filling an array of particles in each simulation with found number of particles
 
  //-------------------------------------------------------------

  /*	This part of code just copies previous arrays
	but it skips the zeros.
  */

  double r_lst[N_mi];								//New array of r without the zeros
  int ID_lst[N_mi];								//New array of ID's without the zeros

  int i_mi = 0;

  for(int i = 0; i < len0; i++)  {
    if(r_lst_temp[i] != 0) {
      r_lst[i_mi] = r_lst_temp[i];						//Filling a new array of r
      i_mi++;
    }
  }

  //-------------------------------------------------------------

  for(int ir = 0; ir < N_mi; ir++) {
    h_nr_temp->Fill(r_lst[ir]);							//Histogram of particles in ring r for all simulations
  }

  //-------------------------------------------------------------

  /*	In this part, a particles density is
	computed and added into appropriate bin.
  */

  for(int ir = 0; ir < N_mi; ir++) {						//Histogram of particles in ring r and number of particles in the shower N_part
    h_nr_temp_2D->Fill(r_lst[ir], N_part_temp_lst[n]);				//in the shower N_part
  }

  int n_i = (h_nr_temp_2D->GetYaxis()->FindBin(N_part_temp_lst[n])) % bins_rho_r;

  sum_n_N_part[n_i] = sum_n_N_part[n_i] + N_part_temp_lst[n];			//Counting sum of all particles of all simulations in this bin
  n_N_part_counter[n_i]++;							//Counting number of simulations in this bin

  //-------------------------------------------------------------

  }										//End SIMULATIONS LOOP

  //-------------------------------------------------------------

  /*	In this part of code, Rho(r) histogram
	and graph over distance r is created.
  */

  //TH1F * h_Rho_avr = new TH1F("h_Rho_avr", " ; ; " , bins_rho_r - 1, 0, r_prc_max);
  TH1F * h_Rho_avr = new TH1F("h_Rho_avr", " ; ; " , bins_rho_r - 1, Rr_bins);

  TGraphErrors * g_Rho_avr = new TGraphErrors("Graph");

  string title_Rho_avr = " #rho_{avr}(r) distribution for E_{casc} = ";
  title_Rho_avr.append(to_string(E_lst[nf]));
  title_Rho_avr.append(" TeV");

  const char* Rho_avr_title = title_Rho_avr.c_str();

  g_Rho_avr->SetTitle(Rho_avr_title);

  g_Rho_avr->GetXaxis()->SetTitle(" r [cm] ");
  g_Rho_avr->GetYaxis()->SetTitle("#rho_{avr}(r) [particles/cm^{2}]");

  //-------------------------------------------------------------

  /*	Here the histograms of Rho_avr(r) and
	Rho_avr(r, E) as well as 2D and 1D graph are
	filled and saved.
  */

  for(int b = 0; b < bins_rho_r; b++) {

    double rho_avr = (1.0/double(n_sim)) * (h_nr_temp->GetBinContent(b)/((2*M_PI)*h_nr_temp->GetBinWidth(b)*h_nr_temp->GetBinCenter(b)));
										//Calculating average particles density Rho_avr for given r bin
    double rho_avr_error_y = (1.0/double(n_sim)) * (sqrt(h_nr_temp->GetBinContent(b))/((2*M_PI)*h_nr_temp->GetBinWidth(b)*h_nr_temp->GetBinCenter(b)));
    double rho_avr_error_x = h_nr_temp->GetBinWidth(b)/2.0;

    h_Rho_avr->SetBinContent(b, rho_avr);					//Filling histogram Rho_avr(r) 1D
    h_Rho_rE_2D->SetBinContent(b, nf, rho_avr);					//Filling histogram Rho_avr(r, E) 2D

    g_Rho_avr->SetPoint(b, h_nr_temp->GetBinCenter(b), rho_avr);		//Filling graph Rho_avr(r) 1D 
    g_Rho_avr->SetPointError(b, rho_avr_error_x, rho_avr_error_y);		//Filling graph errors
									 
    if(rho_avr > 0.0) {
      if(rho_avr < Rho_rE_min) Rho_rE_min = rho_avr;				//Finding minimum Rho(r, E)
      if(rho_avr > Rho_rE_max) Rho_rE_max = rho_avr;				//Finding maximum Rho(r, E)
    }

  }

  g_Rho_avr->SetMinimum(h_Rho_avr->GetBinContent(bins_rho_r - 1));
  g_Rho_avr->GetXaxis()->SetLimits(1.0, r_prc_max);

  TCanvas A("A");
  A.SetLogy();

  string path_Rho_avr = Rho_avr_r_dir;
  if(nf != n_norm) path_Rho_avr.append("/Rho_avr(r)_");
  if((nf == n_norm) || (nf == n_datas - 1)) path_Rho_avr.append("/Rho_avr_norm(r)_");
  path_Rho_avr.append(to_string(E_lst[nf]));
  path_Rho_avr.append(".root");

  const char* Rho_avr_path = path_Rho_avr.c_str();
  g_Rho_avr->SaveAs(Rho_avr_path);						//Saving it in the created directory
  g_Rho_avr->Delete();
  
  //-------------------------------------------------------------

  /*	This part of code fills arrays of
	particles density ratio Rho_(E)/Rho_param
	for different distances r.
  */

  for(int rf = 0; rf < r_number; rf++) {					

    int r_temp = h_Rho_avr->GetXaxis()->FindBin(R_fit_lst[rf]);			//Finding proper bin in Rho_avr(r, E) histogram

    Rho_E_lst[nf][rf] = h_Rho_avr->GetBinContent(r_temp); 			//Filling 2D array of Rho_avr(r, E)
    Rho_E_error[nf][rf] = h_Rho_avr->GetBinErrorLow(r_temp); 			//Filling 2D error array of Rho_avr(r, E)

  }										

  //-------------------------------------------------------------

  /*	2D graphs of Rho(r, N_part) (normalised and unnurmalised)
	are here created.
  */

  string title_Rho_rN_2D = " #rho_{avr}(r, N_{part}) distribution for E_{casc} = ";
  title_Rho_rN_2D.append(to_string(E_lst[nf]));
  title_Rho_rN_2D.append(" TeV ; r [cm]; N_{part} ; #rho(N_{part}) [particles/cm^{2}]");

  const char* Rho_rN_2D_title = title_Rho_rN_2D.c_str();

  string title_Rho_rN_2D_norm = " #rho_{avr}(r, N_{part}) distribution for E_{casc} = ";
  title_Rho_rN_2D_norm.append(to_string(E_lst[nf]));
  title_Rho_rN_2D_norm.append(" TeV ; r [cm]; N_{part}/N_{avr} ; #rho(N_{part})/#rho_{avr}");

  const char* Rho_rN_2D_title_norm = title_Rho_rN_2D_norm.c_str();
  
  long double N_rN_min_norm = (N_rN_min/f_meanN_part_E->Eval(E_lst[nf]));
  long double N_rN_max_norm = (N_rN_max/f_meanN_part_E->Eval(E_lst[nf]));

  TH2 * h_Rho_rN_2D = new TH2F("h_Rho_rN_2D", Rho_rN_2D_title, bins_rho_r - 1, Rr_bins, bins_rho_N_part - 1, N_rN_min, N_rN_max);
  TH2 * h_Rho_rN_2D_norm = new TH2F("h_Rho_rN_2D_norm", Rho_rN_2D_title_norm, bins_rho_r - 1, Rr_bins, bins_rho_N_part - 1, N_rN_min_norm, N_rN_max_norm);

  //-------------------------------------------------------------

  /*	This part of code fills histograms and arrays of
	particles density ratio Rho_(N_part)/Rho_avr
	for different distances r.
  */

  for(int bn = 0; bn < bins_rho_N_part; bn++) {

    N_part_lst[bn] = 0;
    if(n_N_part_counter[bn] > 0.0) N_part_lst[bn] = (sum_n_N_part[bn]/n_N_part_counter[bn]) / f_meanN_part_E->Eval(E_lst[nf]) ;	//Calculating ratio of N_part/N_avr
    N_part_error[bn] = ((N_rN_max_norm - N_rN_min_norm)/bins_rho_N_part)/2.0;							//Evaluating error on x axis (number of particles)



    for(int rn = 0; rn < bins_rho_r; rn++) {					//A loop for histograms - over all bins

    long double rho_temp = 0.0;
    if(n_N_part_counter[bn] > 0.0) rho_temp = (1.0/double(n_N_part_counter[bn]))*h_nr_temp_2D->GetBinContent(rn, bn)/((2*M_PI)*h_nr_temp->GetBinWidth(rn)*h_nr_temp->GetBinCenter(rn));
	    									//Calculating average particles density in the ring
    if(rho_temp > 0.0) {
      if(rho_temp < Rho_rN_min) Rho_rN_min = rho_temp;				//Finding minimum Rho(r, N)
      if(rho_temp > Rho_rN_max) Rho_rN_max = rho_temp;				//Finding maximum Rho(r, N)
    }

    long double rho_norm = 0.0;
    if(h_Rho_avr->GetBinContent(rn) > 0.0) rho_norm = rho_temp / h_Rho_avr->GetBinContent(rn);

    h_Rho_rN_2D->SetBinContent(rn, bn, rho_temp);
    h_Rho_rN_2D_norm->SetBinContent(rn, bn, rho_norm);

    }



    for(int rf = 0; rf < r_number; rf++) {					//A loop for arrays - over all distances in an R_fit_lst array

    int r_temp = h_Rho_avr->GetXaxis()->FindBin(R_fit_lst[rf]);

    long double rho_temp = 0.0;
    if(n_N_part_counter[bn] > 0.0) rho_temp = (1.0/double(n_N_part_counter[bn]))*h_nr_temp_2D->GetBinContent(r_temp, bn)/((2*M_PI)*h_nr_temp->GetBinWidth(r_temp)*h_nr_temp->GetBinCenter(r_temp));
	    									//Calculating average particles density in the ring
    long double rho_norm = 0.0;
    if(h_Rho_avr->GetBinContent(r_temp) > 0.0) rho_norm = rho_temp / h_Rho_avr->GetBinContent(r_temp);

    Rho_N_lst[bn][rf] = rho_norm; 						//Filling 2D array of Rho(N_part)/Rho_avr
    Rho_N_error[bn][rf] = h_Rho_rN_2D_norm->GetBinErrorLow(r_temp, bn);
    if(n_N_part_counter[bn] > 0.0) Rho_N_error[bn][rf] = h_Rho_rN_2D_norm->GetBinErrorLow(r_temp, bn)/sqrt(n_N_part_counter[bn]); 	//Filling 2D error array of Rho_avr(r, N) 

    }

  }

  //-------------------------------------------------------------

  /*	Here, the 2D histograms of Rho(r,N)
	(normalised and unnormalised) are ploted and saved.
  */

  //Rho(N_mi, r) - 2D
  TCanvas J("J");
  J.SetLogx();
  //J.SetLogy();
  J.SetLogz();

  string path_Rho_rN = Rho_rN_dir;
  path_Rho_rN.append("/Rho_(r,N)_");
  path_Rho_rN.append(to_string(E_lst[nf]));
  path_Rho_rN.append(".png");

  const char* Rho_rN_path = path_Rho_rN.c_str();

  h_Rho_rN_2D->SetMinimum(0.5*Rho_rN_min);
  h_Rho_rN_2D->SetMaximum(2.0*Rho_rN_max);

  cout<<endl;
  h_Rho_rN_2D->Draw("lego");
  J.SaveAs(Rho_rN_path);

  h_Rho_rN_2D->Delete();




  TCanvas Jn("Jn");
  Jn.SetLogx();
  //Jn.SetLogy();
  Jn.SetLogz();

  string path_Rho_rN_norm = Rho_rN_norm_dir;
  path_Rho_rN_norm.append("/Rho_(r,N)_norm_");
  path_Rho_rN_norm.append(to_string(E_lst[nf]));
  path_Rho_rN_norm.append(".png");

  const char* Rho_rN_norm_path = path_Rho_rN_norm.c_str();

  cout<<endl;
  h_Rho_rN_2D_norm->Draw("lego");
  Jn.SaveAs(Rho_rN_norm_path);

  h_Rho_rN_2D_norm->Delete();

  h_Rho_avr->Delete();								//Delating histograms
  h_nr_temp->Delete();

  //-------------------------------------------------------------

  /*	Here, the directories for each
	energy for 1D graphs of Rho_(N_part)/Rho_avr
	are created.
  */

  string Rho_rN_1D_norm_E = Rho_rN_1D_norm_dir;					//Naming directories for graphs
  Rho_rN_1D_norm_E.append("/E_");
  Rho_rN_1D_norm_E.append(to_string(E_lst[nf]));

  int check_Rho_norm_rE_1D_E = mkdir(Rho_rN_1D_norm_E.c_str(), 0777);		//Creating directories for graphs

  //-------------------------------------------------------------

  /*	This part creates 1D graphs of Rho_(N_part)/Rho_avr
	saves it into a .root file.
  */

  TCanvas L("L"); 
  L.SetLogy(); 
  //L.SetLogx();  

  for(int rf = 0; rf < r_number; rf++) {					//A loop over all distances r from the R_fit_lst array - Start DISTANCE LOOP

  TGraphErrors * g_Rho_rN_1D = new TGraphErrors("Graph");
  for(int bn = 0; bn < bins_rho_N_part; bn++) {					//Copying density fraction graph

    if(Rho_N_lst[bn][rf] != 0.0) {
      g_Rho_rN_1D->SetPoint(bn, N_part_lst[bn], Rho_N_lst[bn][rf]);
      g_Rho_rN_1D->SetPointError(bn, N_part_error[bn], Rho_N_error[bn][rf]);
    }

  }

  string title_Rho_rN_1D = "#rho(N)/#rho_{avr} distribution for r = ";
  title_Rho_rN_1D.append(to_string(R_fit_lst[rf]));
  title_Rho_rN_1D.append(" cm and E_{casc} = ");
  title_Rho_rN_1D.append(to_string(E_lst[nf]));
  title_Rho_rN_1D.append(" TeV");

  const char* Rho_rN_1D_title = title_Rho_rN_1D.c_str();

  g_Rho_rN_1D->SetTitle(Rho_rN_1D_title);

  g_Rho_rN_1D->SetMinimum(pow(10, -2));
  g_Rho_rN_1D->SetMaximum(pow(10, 2));

  //g_Rho_N->GetXaxis()->SetTitle("E_{casc} [GeV]");
  g_Rho_rN_1D->GetXaxis()->SetTitle("N_{part}/N_{part}_{avr}");
  g_Rho_rN_1D->GetXaxis()->SetLimits(0.0, 2.0);

  g_Rho_rN_1D->GetYaxis()->SetTitle("#rho(N)/#rho_{avr}");

  g_Rho_rN_1D->SetMarkerStyle(8);

  g_Rho_rN_1D->SetMarkerSize(1);



  string name_Rho_rN_1D = Rho_rN_1D_norm_E;
  name_Rho_rN_1D.append("/Rho_(r,N)_1D_norm_");
  name_Rho_rN_1D.append(to_string(R_fit_lst[rf]));
  name_Rho_rN_1D.append(".root");

  const char* Rho_rN_1D_name = name_Rho_rN_1D.c_str();
  g_Rho_rN_1D->SaveAs(Rho_rN_1D_name);						//Saving into .root file
  g_Rho_rN_1D->Delete();

  }										//End DISTANCES LOOP

  //-------------------------------------------------------------

  nf++;

  }
  }
  }										//End FILES LOOP

  //-------------------------------------------------------------

  if((dir = opendir(sim_dir.c_str())) != NULL) {				//Continue only if the path to directory with simulations is correct			

  //-------------------------------------------------------------

  /*	This part creates 1D graphs of Rho(r, E)/Rho_(r, E_norm)
	for given r. It creates a graph for each
	distance r from the R_fit_lst array.
  */

  int E_norm = n_datas - 1;
  E_norm = n_norm;

  for(int rf = 0; rf < r_number; rf++) {					//A loop over all distances r from the R_fit_lst array - Start DISTANCE LOOP

  double Rho_E_temp[n_datas], Rho_E_error_temp[n_datas];

  for(int e_temp = 0; e_temp < n_datas; e_temp++) {

    Rho_E_temp[e_temp] = Rho_E_lst[e_temp][rf]/Rho_E_lst[E_norm][rf];
    Rho_E_error_temp[e_temp] = (Rho_E_error[e_temp][rf]/Rho_E_error[E_norm][rf])/sqrt(n_sim_lst[e_temp]);

  }

  TCanvas N("N"); 
  N.SetLogy(); 
  N.SetLogx();  

  TGraphErrors * g_Rho_E = new TGraphErrors(n_datas, E_lst, Rho_E_temp, E_error, Rho_E_error_temp);

  string title_rho_fit = " #rho_{avr}(E) distribution for r = ";
  title_rho_fit.append(to_string(R_fit_lst[rf]));
  title_rho_fit.append(" cm");

  const char* rho_fit_title = title_rho_fit.c_str();

  g_Rho_E->SetTitle(rho_fit_title);

  //g_Rho_E->SetMinimum(0.5*Rho_E_temp[0]);
  //g_Rho_E->SetMaximum(2*Rho_E_temp[n_datas - 1]);


  //g_Rho_E->GetXaxis()->SetTitle("E_{casc} [GeV]");
  g_Rho_E->GetXaxis()->SetTitle("E_{casc} [TeV]");
  g_Rho_E->GetYaxis()->SetTitle("#rho_{avr}(E)/#rho_{avr}(E_{norm})");

  g_Rho_E->SetMarkerStyle(8);
  g_Rho_E->SetMarkerSize(1);


  string path_Rho_E_r_root = Rho_rE_norm_dir;
  path_Rho_E_r_root.append("/Rho_(E)_norm_");
  path_Rho_E_r_root.append(to_string(R_fit_lst[rf]));
  path_Rho_E_r_root.append(".root");

  const char* Rho_E_r_root_path = path_Rho_E_r_root.c_str();
  g_Rho_E->SaveAs(Rho_E_r_root_path);						//Saving into .root file
  g_Rho_E->Delete();

  }										//End DISTANCES LOOP

  //-------------------------------------------------------------

  /*	This part plots and saves
	2D histogram of Rho_avr(r, E).
  */

  TCanvas Oh("Oh");
  Oh.SetLogx(); 
  //Oh.SetLogy(); 
  Oh.SetLogz(); 

  h_Rho_rE_2D->Draw("lego");
  h_Rho_rE_2D->SetMinimum(0.5*Rho_rE_min);					//Setting minimum
  h_Rho_rE_2D->SetMaximum(2.0*Rho_rE_max);					//Setting maximum

  Oh.SaveAs("Rho_rE_2D.png");							//Saving histogram

  //-------------------------------------------------------------

  /*	Here the 2D graph of Rho_avr(r, E)/Rho_norm
	is filled, ploted and saved.
  */

  int rE_point = 0;								//2D graph iterator
  TGraph2D * g_Rho_rE_2D_norm = new TGraph2D();

  string Rho_rE_2D_tittle = "#rho_{avr}(r, E) distribution; r [cm]; E_{casc} [TeV]; #rho_{avr} [particles/cm^{2}]";

  for(int nf = 0; nf < n_datas; nf++) {
    for(int rf = 0; rf < r_number; rf++) {

      long double rho_frac = Rho_E_lst[nf][rf]/Rho_E_lst[E_norm][rf];		//Calculating average particles density fraction Rho_avr(r, E)/Rho_norm for given r

      g_Rho_rE_2D_norm->SetPoint(rE_point, R_fit_lst[rf], E_lst[nf], rho_frac);	//Filling 2D graph of Rho_avr(r, E)/Rho_norm

      rE_point++;								

    }
  }

  TCanvas O("O");
  O.SetLogx(); 
  O.SetLogy(); 
  O.SetLogz(); 


  g_Rho_rE_2D_norm->SetTitle("#rho_{avr}(r, E) distribution; r [cm]; E_{casc} [TeV]; #rho_{avr}(r, E)/#rho_{avr}(r, E_{norm})");

  cout<<endl;
  g_Rho_rE_2D_norm->Draw();
  O.SaveAs("Rho_rE_2D_norm.png");						//Saving graph
  g_Rho_rE_2D_norm->Delete();							//Delating graphs
  h_Rho_rE_2D->Delete();					
 
  //-------------------------------------------------------------

  }  

  //-------------------------------------------------------------

  auto stop = high_resolution_clock::now(); 
  auto duration = duration_cast<microseconds>(stop - start); 			//Stoping counting time of computations

  cout<<"\n" <<endl;
  cout<< "Time of the computations: " << duration.count() << " * 10^{-6} s = " <<  duration.count()/pow(10, 6) << " s = "  << duration.count()/(pow(10, 6)*60) << " min" <<endl;

  //-------------------------------------------------------------

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~				//End MACRO


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

	- bins_rho_N_part - number of bins in histograms of number of particles N_part (default = 40),
	- bins_rho_r - number of bins in histogram of particles density Rho(r) (default = 40),

	- r_number - number of points on axis of distance from cascade center r (default = 10)

	- p_part_min - minimum momenta of wanted particles [GeV] (default 1),
	- p_part_max - maximum momenta of wanted particles [GeV] (default pow(10, 10) ~ inf),

	- n_datas - number of simulations (default = 18).
*/

void Rho_theta_macro(int ID_criterium[], string sim_dir = "/home/jerzy/CREDO/Analiza/Angle_distribution/Simulations", string plots_dir = get_current_dir_name(), int bins_rho_N_part = 50, int bins_rho_r = 40, int r_number = 10, double p_part_min = 1.0, double p_part_max = pow(10, 10), int n_datas = 8, double max_theta = 70.0)

{

  //-------------------------------------------------------------

  auto start = high_resolution_clock::now(); 					//Starting counting time of computations

  //-------------------------------------------------------------  

  /*	In this part, previously generated file
	is read and several parameters are extracted
	from it.
  */

  double E_min, E_max;
  double r_prc_max, r_rho_max;
  double N_part_min, N_part_max;


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

  else cout<< "ERROR: COULD NOT FIND THE PARAMETERS FILE" <<endl;

  param_file.close();

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

  /*	Here are only definitions of some
	arrays and variables used later.
  */

  double Theta_lst[n_datas], Theta_error[n_datas];				//An array of arrival angle of the simulations Theta
  double E_lst[n_datas], E_error[n_datas];					//An array of energy of the simulations E
  double Npart_lst[n_datas];							//An array of number of parts in the simulations N_part
  double n_sim_lst[n_datas];							//An array of number of simulations in each file
  for(int t; t < n_datas; t++)  Npart_lst[t] = n_sim_lst[t] = 0;				

  long double** RhoTheta_lst = 0;						//An array of rho(theta) of the simulations
  RhoTheta_lst = new long double*[n_datas];

  for(int i = 0; i < n_datas; i++) {
    RhoTheta_lst[i] = new long double[r_number];

    for(int j = 0; j < r_number; j++) {
      RhoTheta_lst[i][j] = 0;
    }
  }

  long double** RhoTheta_error = 0;						//An array of error rho(theta) of the simulations
  RhoTheta_error = new long double*[n_datas];

  for(int i = 0; i < n_datas; i++) {
    RhoTheta_error[i] = new long double[r_number];

    for(int j = 0; j < r_number; j++) {
      RhoTheta_error[i][j] = 0;
    }
  }

  //-------------------------------------------------------------

  /*	Here are defined some directories
	in which plots are stored.
  */

  string Rho_theta_dir_name = "/Rho_(r,theta)_1D";				//Naming directories for histograms
  string Rho_theta_dir = plots_dir;
  Rho_theta_dir.append(Rho_theta_dir_name);

  string Npart_theta_dir_name = "/N_part_(theta)_1D";				//Naming directories for histograms
  string Npart_theta_dir = plots_dir;
  Npart_theta_dir.append(Npart_theta_dir_name);


  int check_Npart_theta = mkdir(Npart_theta_dir.c_str(), 0777);			//Creating directories for histograms
  int check_Rho_theta = mkdir(Rho_theta_dir.c_str(), 0777);			//Creating directories for histograms

  //-------------------------------------------------------------  

  /*	This part of code defines
	a 2D histograms of particles density
	Rho(r, Theta) and its error.
  */

  TH2 * h_Rho_rangle_2D = new TH2F("h_rho_rTheta", " #rho(r, #theta) distribution; r [cm]; #theta ; #rho [cm^{-2}] ", bins_rho_r, 0, r_prc_max, n_datas, 0, (2.0*M_PI/360.0)*max_theta);
  TH2 * h_Rho_rangle_error_2D = new TH2F("h_rho_rTheta_error", " #rho(r, #theta) distribution; r [cm]; #theta ; #rho [cm^{-2}] ", bins_rho_r, 0, r_prc_max, n_datas, 0, (2.0*M_PI/360.0)*max_theta);

  //-------------------------------------------------------------

  int nf = 0;									//Files iterator
  int id_len = sizeof(ID_criterium)/sizeof(ID_criterium[0]);

  //-------------------------------------------------------------

  DIR *dir;
  struct dirent *ent;

  if((dir = opendir(sim_dir.c_str())) != NULL) {				
  while((ent = readdir(dir)) != NULL) {						//A loop over all files in directory - Start FILES LOOP
  if(opendir(ent->d_name) == NULL) {

  //-------------------------------------------------------------

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
  TTree *run = (TTree*) file->Get("run");
  
  TLeafI* l_id = (TLeafI*)sim->GetLeaf("particle..ParticleID");
  TLeaf* l_x = (TLeaf*)sim->GetLeaf("particle..x");
  TLeaf* l_y = (TLeaf*)sim->GetLeaf("particle..y");
  TLeaf* l_pz = (TLeaf*)sim->GetLeaf("particle..Pz");
  TLeaf* l_px = (TLeaf*)sim->GetLeaf("particle..Px");
  TLeaf* l_py = (TLeaf*)sim->GetLeaf("particle..Py");

  TLeaf* E_sim = (TLeaf*)sim->GetLeaf("shower.Energy");
  TLeaf* theta_sim = (TLeaf*)run->GetLeaf("run.ZenithMax");

  //-------------------------------------------------------------   

  /*	This part of code reads angle of arrival of
	the cascades in each file
	(assuming that all showers within the file
	have the same angle) and also atribute
	an error to the energy. It also reads energy of
	the cascades in each file
	(assuming that all showers within the file
	have the same energy). It also atribute
	an error to the energy.
  */

  run->GetEntry(0);
  Theta_lst[nf] = (2.0*M_PI/360.0)*(theta_sim->GetValue(0));			//Filling an array with theta angle

  cout<< "\nTheta angle of showers in the file: " << theta_sim->GetValue(0) << " degrees" <<endl;//Printing angle
  cout<< "Theta angle of showers in the file: " << Theta_lst[nf] << " radians" <<endl;	//Printing angle

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

  /*	In this part of code a 2D and 1D histograms
	of number of particles for different energies
	and angles is created.
  */

  //TH2 * h_NEangle2D_temp = new TH1F(" ", " ; ; ", n_datas, 0, n_datas, n_datas, 0, n_datas);
  TH1 * h_nr_temp = new TH1F(" ", " ; ; ", bins_rho_r, 0, r_prc_max);
  //TH1 * h_Nangle_temp = new TH1F("h_Nangle_temp", " #LT N_{part} #GT (#theta) ; #LT N_{part} #GT; #theta [^o]", n_datas, 0, n_datas);

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

      Npart_lst[nf]++;

      r_lst_temp[N_mi] = sqrt((x*x) + (y*y));					//Filling an array of r
      N_mi++;									//Counting number of particles that filfil given criteria
    }
  }										//End PARTICLES LOOP
 
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

  /*	Here histograms of particles density
	and its error are filled.
  */

  for(int b = 1; b <= bins_rho_r; b++) {

    double rho_avr = (1.0/double(n_sim)) * (h_nr_temp->GetBinContent(b)/((2*M_PI)*h_nr_temp->GetBinWidth(b)*h_nr_temp->GetBinCenter(b)));
										//Calculating average particles density Rho_avr for given r bin
    double rho_avr_error_y = (1.0/double(n_sim)) * (sqrt(h_nr_temp->GetBinContent(b))/((2*M_PI)*h_nr_temp->GetBinWidth(b)*h_nr_temp->GetBinCenter(b)));
    double rho_avr_error_x = h_nr_temp->GetBinWidth(b)/2.0;

    int nt_temp = h_Rho_rangle_2D->GetYaxis()->FindBin(Theta_lst[nf]);		//Finding proper theta = 0 bin in Rho_avr(r,  Theta) histogram

    h_Rho_rangle_2D->SetBinContent(b, nt_temp, rho_avr);			//Filling histogram Rho_avr(r, Theta) 2D
    h_Rho_rangle_error_2D->SetBinContent(b, nt_temp, rho_avr_error_y);		//Filling histogram error Rho_avr(r, Theta) 2D
								 
  }										

  //-------------------------------------------------------------

  }										//End SIMULATIONS LOOP

  //-------------------------------------------------------------

  nf++;

  }
  }
  }										//End FILES LOOP

  //-------------------------------------------------------------

  /*	Here are defined some directories
	in which plots are stored.
  */

  string Rho_theta_dir_nameE = "/Rho_(r,theta)_";
  Rho_theta_dir_nameE.append(to_string(E_lst[0]));				//Naming directories for histograms
  string Rho_theta_dirE = Rho_theta_dir;
  Rho_theta_dirE.append(Rho_theta_dir_nameE);

  int check_Rho_thetaE = mkdir(Rho_theta_dirE.c_str(), 0777);			//Creating directories for histograms

  //-------------------------------------------------------------

  /*	This part of code creates a graph
	for each distance from the list. The graphs
	shows Rhop(theta) relationship.
  */

  for(int br = 0; br < r_number; br++) {

    double norm_max = 0.0;							//Maximum density for graph maximum

    TGraphErrors * g_Rho_theta = new TGraphErrors("Graph");

    cout<< R_fit_lst[br] <<endl;
    string title_Rho_theta = " #rho_{avr}(#theta) distribution for r = ";
    title_Rho_theta.append(to_string(R_fit_lst[br]));
    title_Rho_theta.append(" cm");

    const char* Rho_theta_title = title_Rho_theta.c_str();

    g_Rho_theta->SetTitle(Rho_theta_title);

    g_Rho_theta->GetXaxis()->SetTitle(" #theta [^o] ");
    g_Rho_theta->GetYaxis()->SetTitle(" #rho_{avr}(#theta) [particles/cm^{2}] ");

    int br_temp = h_Rho_rangle_2D->GetXaxis()->FindBin(R_fit_lst[br]);		//Finding proper r bin in Rho_avr(r,  Theta) histogram

    for(int nt = 0; nt < n_datas; nt++) {

      int nt_temp = h_Rho_rangle_2D->GetYaxis()->FindBin(Theta_lst[nt]);	//Finding proper theta bin in Rho_avr(r,  Theta) histogram
      int nt_param = h_Rho_rangle_2D->GetYaxis()->FindBin(0.0);			//Finding proper theta = 0 bin in Rho_avr(r,  Theta) histogram

      long double rho_norm = (h_Rho_rangle_2D->GetBinContent(br_temp, nt_temp))/(h_Rho_rangle_2D->GetBinContent(br_temp, nt_param));
      long double error_norm = (h_Rho_rangle_error_2D->GetBinContent(br_temp, nt_temp))/(h_Rho_rangle_2D->GetBinContent(br_temp, nt_temp));

      long double theta_err_temp = 1.0;
      for(int ner = 0; ner < n_datas; ner++) {					//Finding smallest error
        if((abs(Theta_lst[nt] - Theta_lst[ner]) < theta_err_temp) && (ner != nt)) theta_err_temp = abs(Theta_lst[nt] - Theta_lst[ner]);
      }

      Theta_error[nt] = theta_err_temp/2.0;					//Filling an array with theta angle error
     
      if(rho_norm > norm_max) norm_max = rho_norm;				//Finding maximum particles density

      g_Rho_theta->SetPoint(nt, Theta_lst[nt], rho_norm);			//Filling graph Rho_avr(r) 1D 
      g_Rho_theta->SetPointError(nt, Theta_error[nt], error_norm);		//Filling graph errors
  
    }

    g_Rho_theta->SetMinimum(0.0);
    g_Rho_theta->SetMaximum(1.1*norm_max);
    g_Rho_theta->GetXaxis()->SetLimits(-0.1, M_PI/2.0);

    TCanvas A("A");
    //A.SetLogy();

    string path_Rho_theta = Rho_theta_dirE;
    path_Rho_theta.append("/Rho_avr(theta)_");
    path_Rho_theta.append(to_string(R_fit_lst[br]));
    path_Rho_theta.append(".root");

    const char* Rho_theta_path = path_Rho_theta.c_str();
    g_Rho_theta->SaveAs(Rho_theta_path);					//Saving it in the created directory

  }

  //-------------------------------------------------------------

  /*	Here 1D graph of number
	of parts N_part(Theta) is created and
	filled.
  */

  TGraphErrors * g_Npart_theta = new TGraphErrors("Graph");

  g_Npart_theta->SetTitle("#LT N_{part} #GT (#theta)");

  g_Npart_theta->GetXaxis()->SetTitle("#theta [^o]");
  g_Npart_theta->GetYaxis()->SetTitle("#LT N_{part} #GT");

  double Nnorm_max = 0.0;							//Maximum Npart for graph maximum

  for(int nt = 0; nt < n_datas; nt++) {

    g_Npart_theta->SetPoint(nt, Theta_lst[nt], Npart_lst[nt]/Npart_lst[1]);			//Filling graph Rho_avr(r) 1D 
    g_Npart_theta->SetPointError(nt, Theta_error[nt], sqrt(Npart_lst[nt])/Npart_lst[1]);	//Filling graph errors

    if(Npart_lst[nt] > Nnorm_max) Nnorm_max = Npart_lst[nt];					//Finding maximum N_part
  
  }

  int nt_temp = h_Rho_rangle_2D->GetYaxis()->FindBin(0.0);			//Finding proper theta bin in Rho_avr(r,  Theta) histogram

  g_Npart_theta->SetMinimum(0.0);
  g_Npart_theta->SetMaximum(1.1*Nnorm_max/Npart_lst[nt_temp]);
  g_Npart_theta->GetXaxis()->SetLimits(-0.1, M_PI/2.0);

  TCanvas B("B");
  //A.SetLogy();

  string path_Npart_theta = Npart_theta_dir;
  path_Npart_theta.append("/Npart_(theta)_");
  path_Npart_theta.append(to_string(E_lst[0]));
  path_Npart_theta.append(".png");
  path_Npart_theta.append(".root");

  const char* Npart_theta_path = path_Npart_theta.c_str();
  g_Npart_theta->SaveAs(Npart_theta_path);					//Saving it in the created directory

  //-------------------------------------------------------------

  auto stop = high_resolution_clock::now(); 
  auto duration = duration_cast<microseconds>(stop - start); 			//Stoping counting time of computations

  cout<<"\n" <<endl;
  cout<< "Time of the computations: " << duration.count() << " * 10^{-6} s = " <<  duration.count()/pow(10, 6) << " s = "  << duration.count()/(pow(10, 6)*60) << " min" <<endl;

  //-------------------------------------------------------------

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~				//End MACRO


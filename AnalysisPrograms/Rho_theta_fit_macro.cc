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

/* 	This macro fits functions
	to graphs created later for later use.

	Its starting parameters are:
	- plots_dir - path to directory where the plots are stored,
	- r_number - number of distance points analysed.
*/

void Rho_theta_fit_macro(string plots_dir = get_current_dir_name(), int r_number = 10)	//Start MACRO
{

  //-------------------------------------------------------------

  auto start = high_resolution_clock::now(); 					//Starting counting time of computations

  //-------------------------------------------------------------

  /*	In this part the functions of rho(N)/rho_avr(E)
	for different distances r are defined.
  */

  TF1 * f_Rho_rtheta = new TF1("f_Rho_rtheta" , " [0]*pow(cos([1]*(x + [2])), [3]) ");

  //-------------------------------------------------------------

  string plots_dir_Rho_rtheta = plots_dir;
  plots_dir_Rho_rtheta.append("/Rho_(r,theta)_1D/");

  cout<< plots_dir_Rho_rtheta <<endl;

  int nf = 0;
  int nr = 0;

  //-------------------------------------------------------------

  /*	Here the graphs of Rho(theta)/Rho(0) are loaded
	from directories for each energy in which it were
	previously saved as .root files.
  */

  ofstream param_Rho_theta_file("Rho(r, theta)_fit_params.txt");

  DIR *dir;
  struct dirent *ent;

  if((dir = opendir(plots_dir_Rho_rtheta.c_str())) != NULL) {				
  while((ent = readdir(dir)) != NULL) {						//A loop over all directories in the directory - Start DIRECTORY LOOP

  string dir_name = plots_dir_Rho_rtheta;
  dir_name.append(ent->d_name);

  cout<< "\nDirectory or file name nr "<< nf+1 << " : " <<endl;
  cout<< dir_name <<endl;

  nr = 0;

  //-------------------------------------------------------------

  DIR * dir2;
  struct dirent * ent2;

  if((dir2 = opendir(dir_name.c_str())) != NULL) {				
  while((ent2 = readdir(dir2)) != NULL) {					//A loop over all files in the directory - Start FILES LOOP
  if(opendir(ent2->d_name) == NULL) {

  cout<< "\nFile name: " <<endl;
  cout<< ent2->d_name <<endl;

  string file_graph = dir_name;
  file_graph.append("/");
  file_graph.append(ent2->d_name);						//Adding name of the file into directory path

  if(file_graph.find(".root") > file_graph.size()) continue;			//Checking if its a root file

  double Rho_rtheta_p0_lst[r_number], Rho_rtheta_p1_lst[r_number], Rho_rtheta_p2_lst[r_number], Rho_rtheta_p3_lst[r_number];	//Arrays of parameters of all functions

  cout<< "\nOpening file nr " << nr+1 << " : " <<endl;
  cout<< file_graph <<endl;
  
  TFile *file = new TFile(file_graph.c_str(), "READ");				//Opening each file
  TGraphErrors * g_Rho_rtheta = (TGraphErrors *) file->Get("Graph");		//Getting a graph from the file

  //-------------------------------------------------------------

  /*	In this part the functions are fitted
	and their parameters are saved into
	a Rho(r, theta)_fit_params.txt file.
  */

    f_Rho_rtheta->SetParameters(1.0, 0.1, -0.1, 30.0);
    //if(nr > 0) f_Rho_rtheta->SetParameters(Rho_rtheta_p0_lst[nr-1], Rho_rtheta_p1_lst[nr-1], Rho_rtheta_p2_lst[nr-1], Rho_rtheta_p3_lst[nr-1]);

    TFitResultPtr result_Rho_rtheta = g_Rho_rtheta->Fit(f_Rho_rtheta, "SQ");  
    Rho_rtheta_p0_lst[nr] = result_Rho_rtheta->Parameter(0);
    Rho_rtheta_p1_lst[nr] = result_Rho_rtheta->Parameter(1);
    Rho_rtheta_p2_lst[nr] = result_Rho_rtheta->Parameter(2);
    Rho_rtheta_p3_lst[nr] = result_Rho_rtheta->Parameter(3);

    param_Rho_theta_file << Rho_rtheta_p0_lst[nr] << " ";
    param_Rho_theta_file << Rho_rtheta_p1_lst[nr] << " ";
    param_Rho_theta_file << Rho_rtheta_p2_lst[nr] << " ";
    param_Rho_theta_file << Rho_rtheta_p3_lst[nr] << "\n";


  //-------------------------------------------------------------

    TCanvas F1("F"); 
    //F1.SetLogy(); 

    string file_name = ent2->d_name;
    int size_name = file_name.size();
    file_name.resize(size_name - 5);

    cout<<endl;
    string path_Rho_rtheta_fit = dir_name;
    path_Rho_rtheta_fit.append("/");
    path_Rho_rtheta_fit.append(file_name);
    path_Rho_rtheta_fit.append("_fit.png");

    const char* Rho_rtheta_fit_path = path_Rho_rtheta_fit.c_str();

    g_Rho_rtheta->Draw("AP");
    F1.SaveAs(Rho_rtheta_fit_path);
    g_Rho_rtheta->Delete();

    nr++;

  //-------------------------------------------------------------

  }
  }
  }										//End FILES LOOP

  nf++;
  param_Rho_theta_file << "\n";

  //-------------------------------------------------------------

  }
  }										//End DIRECTORIES LOOP

  param_Rho_theta_file.close();

  //-------------------------------------------------------------

  /*	In this part the function of
	Npart(theta) is defined.
  */

  TF1 * f_Npart_theta = new TF1("f_Npart_theta" , " [0]*pow(cos([1]*(x + [2])), [3]) ");

  //-------------------------------------------------------------

  string plots_dir_Npart_theta = plots_dir;
  plots_dir_Npart_theta.append("/N_part_(theta)_1D/");

  nr = 0;

  double Npart_theta_p0_lst[r_number], Npart_theta_p1_lst[r_number], Npart_theta_p2_lst[r_number], Npart_theta_p3_lst[r_number]; //Arrays of parameters of all functions

  //-------------------------------------------------------------

  /*	Here the graphs of Npart(theta)/Npart(0) are
	loaded from the directory in
	which it were previously saved as .root files.
  */

  ofstream param_Npart_theta_file("Npart(theta)_fit_params.txt");

  if((dir = opendir(plots_dir_Npart_theta.c_str())) != NULL) {				
  while((ent = readdir(dir)) != NULL) {						//A loop over all files in directory - Start FILES LOOP
  if(opendir(ent->d_name) == NULL) {

								
  cout<< "\nFile name: " <<endl;
  cout<< ent->d_name <<endl;

  string file_graph = plots_dir_Npart_theta;
  file_graph.append(ent->d_name);						//Adding name of the file into directory path


  if(file_graph.find(".root") > file_graph.size()) continue;			//Checking if its a root file

  cout<< "\nOpening file nr " << nf+1 << " : " <<endl;
  cout<< file_graph <<endl;
  
  TFile * file = new TFile(file_graph.c_str(), "READ");				//Opening each file
  TGraphErrors* g_Npart_theta = (TGraphErrors *) file->Get("Graph");		//Getting a graph from the file

  //-------------------------------------------------------------

  /*	In this part the functions are fitted
	and their parameters are saved into
	a Npart(theta)_fit_params.txt file.
  */

    f_Npart_theta->SetParameters(1.0, 0.1, -0.1, 30.0);
    //if(nr > 0) f_Npart_theta->SetParameters(Npart_theta_p0_lst[nr-1], Npart_theta_p1_lst[nr-1], Npart_theta_p2_lst[nr-1], Npart_theta_p3_lst[nr-1]);

    TFitResultPtr result_Npart_theta = g_Npart_theta->Fit(f_Npart_theta, "SQ");  
    Npart_theta_p0_lst[nr] = result_Npart_theta->Parameter(0);
    Npart_theta_p1_lst[nr] = result_Npart_theta->Parameter(1);
    Npart_theta_p2_lst[nr] = result_Npart_theta->Parameter(2);
    Npart_theta_p3_lst[nr] = result_Npart_theta->Parameter(3);

    param_Npart_theta_file << Npart_theta_p0_lst[nr] << " ";
    param_Npart_theta_file << Npart_theta_p1_lst[nr] << " ";
    param_Npart_theta_file << Npart_theta_p2_lst[nr] << " ";
    param_Npart_theta_file << Npart_theta_p3_lst[nr] << "\n";

  //-------------------------------------------------------------

  TCanvas F2("F"); 
  //F2.SetLogy(); 
  //F2.SetLogx(); 

  string file_name = ent->d_name;
  int size_name = file_name.size();
  file_name.resize(size_name - 5);

  cout<<endl;
  string path_Npart_theta_fit = plots_dir_Npart_theta;
  path_Npart_theta_fit.append(file_name);
  path_Npart_theta_fit.append("_fit.png");

  const char* Npart_theta_fit_path = path_Npart_theta_fit.c_str();

  g_Npart_theta->Draw("AP");
  F2.SaveAs(Npart_theta_fit_path);
  g_Npart_theta->Delete();

  nr++;

  //-------------------------------------------------------------

  }
  }
  }										//End FILES LOOP

  param_Npart_theta_file.close();

  //-------------------------------------------------------------

  auto stop = high_resolution_clock::now(); 
  auto duration = duration_cast<microseconds>(stop - start); 			//Stoping counting time of computations

  cout<<"\n" <<endl;
  cout<< "Time of the computations: " << duration.count() << " * 10^{-6} s = " <<  duration.count()/pow(10, 6) << " s = "  << duration.count()/(pow(10, 6)*60) << " min" <<endl;

  //-------------------------------------------------------------

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~				//End MACRO


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

void Rho_rNE_fit_macro(string plots_dir = get_current_dir_name(), int r_number = 10)	//Start MACRO
{

  //-------------------------------------------------------------

  auto start = high_resolution_clock::now(); 					//Starting counting time of computations

  //-------------------------------------------------------------

  /*	In this part the functions of rho(N)/rho_avr(E)
	for different distances r are defined.
  */

  TF1 * f_Rho_rN = new TF1("f_Rho_rN" , " [0]*x + [1]");

  //-------------------------------------------------------------

  string plots_dir_Rho_rN = plots_dir;
  plots_dir_Rho_rN.append("/Rho_(r,N)_1D_norm/");

  cout<< plots_dir_Rho_rN <<endl;

  int nf = 0;
  int nr = 0;

  //-------------------------------------------------------------

  /*	Here the graphs of rho(N)/rho_avr(E) are loaded
	from directories for each energy in which it were
	previously saved as .root files.
  */

  ofstream param_rN_file("Rho(r,N)_fit_params.txt");

  DIR *dir;
  struct dirent *ent;

  if((dir = opendir(plots_dir_Rho_rN.c_str())) != NULL) {				
  while((ent = readdir(dir)) != NULL) {						//A loop over all directories in the directory - Start DIRECTORY LOOP

  string dir_name = plots_dir_Rho_rN;
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

  double Rho_rN_p0_lst[r_number], Rho_rN_p1_lst[r_number];			//Arrays of parameters of all functions

  cout<< "\nOpening file nr " << nr+1 << " : " <<endl;
  cout<< file_graph <<endl;
  
  TFile *file = new TFile(file_graph.c_str(), "READ");				//Opening each file
  TGraphErrors * g_Rho_rN = (TGraphErrors *) file->Get("Graph");		//Getting a graph from the file

  //-------------------------------------------------------------

  /*	In this part the functions are fitted
	and their parameters are saved into
	a Rho(r,N)_fit_params.txt file.
  */

    f_Rho_rN->SetParameters(1.0, 0.0);

    TFitResultPtr result_Rho_rN = g_Rho_rN->Fit(f_Rho_rN, "SQ");  
    Rho_rN_p0_lst[nr] = result_Rho_rN->Parameter(0);
    Rho_rN_p1_lst[nr] = result_Rho_rN->Parameter(1);

    param_rN_file << Rho_rN_p0_lst[nr] << " ";
    param_rN_file << Rho_rN_p1_lst[nr] << "\n";


  //-------------------------------------------------------------

    TCanvas F1("F"); 
    F1.SetLogy(); 

    string file_name = ent2->d_name;
    int size_name = file_name.size();
    file_name.resize(size_name - 5);

    cout<<endl;
    string path_Rho_rN_fit = dir_name;
    path_Rho_rN_fit.append("/");
    path_Rho_rN_fit.append(file_name);
    path_Rho_rN_fit.append("_fit.png");

    const char* Rho_rN_fit_path = path_Rho_rN_fit.c_str();

    g_Rho_rN->Draw("AP");
    F1.SaveAs(Rho_rN_fit_path);
    g_Rho_rN->Delete();

    nr++;

  //-------------------------------------------------------------

  }
  }
  }										//End FILES LOOP

  nf++;
  param_rN_file << "\n";

  //-------------------------------------------------------------

  }
  }										//End DIRECTORIES LOOP

  param_rN_file.close();

  //-------------------------------------------------------------

  /*	In this part the functions of meanR_rho(E),
	sigmaR_rho(E) and cR_rho(E) are defined.
  */

  TF1 * f_Rho_E = new TF1("f_Rho_E" , " [0]*pow(x, [1]) + [2]");

  //-------------------------------------------------------------

  string plots_dir_Rho_E = plots_dir;
  plots_dir_Rho_E.append("/Rho_(r,E)_1D_norm/");

  nr = 0;

  double Rho_E_p0_lst[r_number], Rho_E_p1_lst[r_number], Rho_E_p2_lst[r_number]; //Arrays of parameters of all functions

  //-------------------------------------------------------------

  /*	Here the graphs of rho(E)/rho_(E_param) are
	loaded from the directory in
	which it were previously saved as .root files.
  */

  ofstream param_rE_file("Rho(r,E)_fit_params.txt");

  if((dir = opendir(plots_dir_Rho_E.c_str())) != NULL) {				
  while((ent = readdir(dir)) != NULL) {						//A loop over all files in directory - Start FILES LOOP
  if(opendir(ent->d_name) == NULL) {

								
  cout<< "\nFile name: " <<endl;
  cout<< ent->d_name <<endl;

  string file_graph = plots_dir_Rho_E;
  file_graph.append(ent->d_name);						//Adding name of the file into directory path


  if(file_graph.find(".root") > file_graph.size()) continue;			//Checking if its a root file

  cout<< "\nOpening file nr " << nf+1 << " : " <<endl;
  cout<< file_graph <<endl;
  
  TFile * file = new TFile(file_graph.c_str(), "READ");				//Opening each file
  TGraphErrors* g_Rho_E = (TGraphErrors *) file->Get("Graph");			//Getting a graph from the file

  //-------------------------------------------------------------

  /*	In this part the functions are fitted
	and their parameters are saved into
	a Rho(r,E)_fit_params.txt file.
  */

    f_Rho_E->SetParameters(3e-04, 1.0, -1e-04);

    TFitResultPtr result_Rho_E = g_Rho_E->Fit(f_Rho_E, "SQ");  
    Rho_E_p0_lst[nr] = result_Rho_E->Parameter(0);
    Rho_E_p1_lst[nr] = result_Rho_E->Parameter(1);
    Rho_E_p2_lst[nr] = result_Rho_E->Parameter(2);

    param_rE_file << Rho_E_p0_lst[nr] << " ";
    param_rE_file << Rho_E_p1_lst[nr] << " ";
    param_rE_file << Rho_E_p2_lst[nr] << "\n";

  //-------------------------------------------------------------

  TCanvas F2("F"); 
  F2.SetLogy(); 
  F2.SetLogx(); 

  string file_name = ent->d_name;
  int size_name = file_name.size();
  file_name.resize(size_name - 5);

  cout<<endl;
  string path_Rho_E_fit = plots_dir_Rho_E;
  path_Rho_E_fit.append(file_name);
  path_Rho_E_fit.append("_fit.png");

  const char* Rho_E_fit_path = path_Rho_E_fit.c_str();

  g_Rho_E->Draw("AP");
  F2.SaveAs(Rho_E_fit_path);
  g_Rho_E->Delete();

  nr++;

  //-------------------------------------------------------------

  }
  }
  }										//End FILES LOOP

  param_rE_file.close();

  //-------------------------------------------------------------

  /*	In this part the function of Rho_avr(r)
	is defined.
  */

  TF1 * f_Rho_avr = new TF1("f_Rho_avr" , " [0]/(pow((x + [1]), [2])) ");

  //-------------------------------------------------------------

  string plots_dir_Rho_avr = plots_dir;
  plots_dir_Rho_avr.append("/Rho_avr(r)/");

  nf = 0;

  double Rho_avr_p0_lst[r_number], Rho_avr_p1_lst[r_number], Rho_avr_p2_lst[r_number]; //Arrays of parameters of all functions

  //-------------------------------------------------------------

  /*	Here the graphs of Rho_avr are
	loaded from the directory in
	which it were previously saved as .root files.
  */

  ofstream param_avr_file("Rho(r)_fit_params.txt");

  if((dir = opendir(plots_dir_Rho_avr.c_str())) != NULL) {				
  while((ent = readdir(dir)) != NULL) {						//A loop over all files in directory - Start FILES LOOP
  if(opendir(ent->d_name) == NULL) {

								
  cout<< "\nFile name: " <<endl;
  cout<< ent->d_name <<endl;

  string file_graph = plots_dir_Rho_avr;
  file_graph.append(ent->d_name);						//Adding name of the file into directory path


  if(file_graph.find(".root") > file_graph.size()) continue;			//Checking if its a root file

  cout<< "\nOpening file nr " << nf+1 << " : " <<endl;
  cout<< file_graph <<endl;
  
  TFile * file = new TFile(file_graph.c_str(), "READ");				//Opening each file
  TGraphErrors * g_Rho_avr = (TGraphErrors *) file->Get("Graph");		//Getting a graph from the file


  //-------------------------------------------------------------

  /*	In this part the functions are fitted
	and their parameters are saved into
	a Rho(r)_fit_params.txt file.
  */

    f_Rho_avr->SetParameters(0.7*(nf+1)*9e+09, 3.5e+04, 4.0);

    TFitResultPtr result_rho_avr = g_Rho_avr->Fit(f_Rho_avr, "SQ");  
    Rho_avr_p0_lst[nf] = result_rho_avr->Parameter(0);
    Rho_avr_p1_lst[nf] = result_rho_avr->Parameter(1);
    Rho_avr_p2_lst[nf] = result_rho_avr->Parameter(2);

    param_avr_file << Rho_avr_p0_lst[nf] << " ";
    param_avr_file << Rho_avr_p1_lst[nf] << " ";
    param_avr_file << Rho_avr_p2_lst[nf] << "\n";

  //-------------------------------------------------------------

  TCanvas F3("F"); 
  F3.SetLogy(); 

  string file_name = ent->d_name;
  int size_name = file_name.size();
  file_name.resize(size_name - 5);

  cout<<endl;
  string path_Rho_avr_fit = plots_dir_Rho_avr;
  path_Rho_avr_fit.append(file_name);
  path_Rho_avr_fit.append("_fit.png");

  const char* Rho_avr_fit_path = path_Rho_avr_fit.c_str();

  g_Rho_avr->Draw("AP");
  F3.SaveAs(Rho_avr_fit_path);
  g_Rho_avr->Delete();

  nf++;

  //-------------------------------------------------------------

  }
  }
  }										//End FILES LOOP

  param_rE_file.close();

  //-------------------------------------------------------------

  auto stop = high_resolution_clock::now(); 
  auto duration = duration_cast<microseconds>(stop - start); 			//Stoping counting time of computations

  cout<<"\n" <<endl;
  cout<< "Time of the computations: " << duration.count() << " * 10^{-6} s = " <<  duration.count()/pow(10, 6) << " s = "  << duration.count()/(pow(10, 6)*60) << " min" <<endl;

  //-------------------------------------------------------------

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~				//End MACRO


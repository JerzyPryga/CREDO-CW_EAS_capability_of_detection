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

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/* 	Used functions */

string get_current_dir_fit_macro();						//Function returning path to current directory

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/* 	The second macro fits functions
	to graphs created later for later use.

	Its starting parameter is a path to directory where the plots are stored.
*/

void NR_E_fit_macro(string plots_dir = get_current_dir_fit_macro())				//Start MACRO
{

  //-------------------------------------------------------------

  /*	In this part the functions of meanN_part(E),
	sigmaN_part(E) and cN_part(E) are defined.
  */

  TF1 * f_meanN_part_E = new TF1("f_meanN_part_E" , " [0]*pow(x, [1]) + [2]");
  f_meanN_part_E->SetParameters(15.5473, 0.91507, 0.0322778);

  TF1 * f_sigmaN_part_E = new TF1("f_sigmaN_part_E" , " [0]*pow(x, [1]) + [2]");
  f_sigmaN_part_E->SetParameters(4.71603, 0.859976, 5.07782);

  TF1 * f_cN_part_E = new TF1("f_cN_part_E" , " [0]*pow(x, [1]) + [2]");
  f_cN_part_E->SetParameters(0.00736672, 0.884137, 105.498);

  TF1 * f_fracN_part_E = new TF1("f_fracN_part_E" , " [0]*pow(x, [1]) + [2]");
  f_fracN_part_E->SetParameters(0.00736672, 0.884137, 105.498);

  //-------------------------------------------------------------

  string plots_dir_N_part = plots_dir;
  plots_dir_N_part.append("/N_part(E)/");

  int nf = 0;

  //-------------------------------------------------------------

  /*	Here the graphs of meanN_part(E), sigmaN_part(E),
	cN_part(E) and sigmaN_part/meanN_part(E) are loaded
	from the directory in which it were previously saved
	as .root files.
  */

  ofstream param_file ("N(E)_fit_params.txt");

  DIR *dir;
  struct dirent *ent;

  if((dir = opendir(plots_dir_N_part.c_str())) != NULL) {				
  while((ent = readdir(dir)) != NULL) {						//A loop over all files in directory - Start FILES LOOP
  if(opendir(ent->d_name) == NULL) {

								
  cout<< "\nFile name: " <<endl;
  cout<< ent->d_name <<endl;

  string file_graph = plots_dir_N_part;
  file_graph.append(ent->d_name);						//Adding name of the file into directory path


  if(file_graph.find(".root") > file_graph.size()) continue;			//Checking if its a root file

  cout<< "\nOpening file nr " << nf+1 << " : " <<endl;
  cout<< file_graph <<endl;
  
  TFile *file = new TFile(file_graph.c_str(), "READ");				//Opening each file
  TGraph * g_N_part_E = (TGraph *) file->Get("Graph");				//Getting a graph from the file

  //-------------------------------------------------------------

  /*	In this part the functions are fitted
	and their parameters are saved into
	a N(E)_params.txt file.
  */

  //param_file.open("N(E)_params.txt");

  if(file_graph.find("cN_part(E)") < file_graph.size()) {

    TFitResultPtr result_N_part_E = g_N_part_E->Fit(f_cN_part_E, "S");  
    param_file << result_N_part_E->Parameter(0) << "\n";
    param_file << result_N_part_E->Parameter(1) << "\n";
    param_file << result_N_part_E->Parameter(2) << "\n";

  }

  if(file_graph.find("sigmaN_part(E)") < file_graph.size()) {

    TFitResultPtr result_N_part_E = g_N_part_E->Fit(f_sigmaN_part_E, "S");  
    param_file << result_N_part_E->Parameter(0) << "\n";
    param_file << result_N_part_E->Parameter(1) << "\n";
    param_file << result_N_part_E->Parameter(2) << "\n";

  }

  if(file_graph.find("meanN_part(E)") < file_graph.size()) {

    TFitResultPtr result_N_part_E = g_N_part_E->Fit(f_meanN_part_E, "S");  
    param_file << result_N_part_E->Parameter(0) << "\n";
    param_file << result_N_part_E->Parameter(1) << "\n";
    param_file << result_N_part_E->Parameter(2) << "\n";

  }

  if(file_graph.find("N_part_frac(E)") < file_graph.size()) {

    TFitResultPtr result_N_part_E = g_N_part_E->Fit(f_fracN_part_E, "S");  

  }

  //-------------------------------------------------------------

  TCanvas F1("F"); 
  F1.SetLogy(); 
  F1.SetLogx(); 

  string file_name = ent->d_name;
  int size_name = file_name.size();
  file_name.resize(size_name - 5);

  cout<<endl;
  string path_N_part_fit = plots_dir_N_part;
  path_N_part_fit.append(file_name);
  path_N_part_fit.append("_fit.png");

  const char* N_part_fit_path = path_N_part_fit.c_str();

  g_N_part_E->Draw("AP");
  F1.SaveAs(N_part_fit_path);
  g_N_part_E->Delete();

  nf++;

  //-------------------------------------------------------------

  }
  }
  }										//End FILES LOOP

  param_file.close();

  //-------------------------------------------------------------


  /*	In this part the functions of meanR_prc(E),
	sigmaR_prc(E) and cR_prc(E) are defined.
  */

  TF1 * f_meanR_prc_E = new TF1("f_meanR_prc_E" , " [0]*pow(x, [1]) + [2]");
  f_meanR_prc_E->SetParameters(-3.32045e+07, 0.000957474, 4.77929e+06);

  TF1 * f_sigmaR_prc_E = new TF1("f_sigmaR_prc_E" , " [0]*pow(x, [1]) + [2]");
  f_sigmaR_prc_E->SetParameters(-3.32045e+07, 0.000137065, 3.32486e+07);

  TF1 * f_cR_prc_E = new TF1("f_cR_prc_E" , " [0]*pow(x, [1]) + [2]");
  f_cR_prc_E->SetParameters(-3.32045e+07, 0.000137065, 3.32486e+07);

  //-------------------------------------------------------------

  string plots_dir_R_prc = plots_dir;
  plots_dir_R_prc.append("/R_prc(E)/");

  nf = 0;

  //-------------------------------------------------------------

  /*	Here the graphs of meanR_prc(E), sigmaR_prc(E) and
	cR_prc(E) are loaded from the directory in
	which it were previously saved as .root files.
  */

  if((dir = opendir(plots_dir_R_prc.c_str())) != NULL) {				
  while((ent = readdir(dir)) != NULL) {						//A loop over all files in directory - Start FILES LOOP
  if(opendir(ent->d_name) == NULL) {

								
  cout<< "\nFile name: " <<endl;
  cout<< ent->d_name <<endl;

  string file_graph = plots_dir_R_prc;
  file_graph.append(ent->d_name);						//Adding name of the file into directory path


  if(file_graph.find(".root") > file_graph.size()) continue;			//Checking if its a root file

  cout<< "\nOpening file nr " << nf+1 << " : " <<endl;
  cout<< file_graph <<endl;
  
  TFile *file = new TFile(file_graph.c_str(), "READ");				//Opening each file
  TGraph * g_R_prc_E = (TGraph *) file->Get("Graph");				//Getting a graph from the file

  //-------------------------------------------------------------

  /*	In this part the functions are fitted.
  */


  if(file_graph.find("cR_prc(E)") < file_graph.size()) {

    TFitResultPtr result_R_prc_E = g_R_prc_E->Fit(f_cR_prc_E, "S");  

  }

  if(file_graph.find("sigmaR_prc(E)") < file_graph.size()) {

    TFitResultPtr result_R_prc_E = g_R_prc_E->Fit(f_sigmaR_prc_E, "S");  

  }

  if(file_graph.find("meanR_prc(E)") < file_graph.size()) {

    TFitResultPtr result_R_prc_E = g_R_prc_E->Fit(f_meanR_prc_E, "S");  

  }

  //-------------------------------------------------------------

  TCanvas F2("F"); 
  F2.SetLogy(); 
  F2.SetLogx(); 

  string file_name = ent->d_name;
  int size_name = file_name.size();
  file_name.resize(size_name - 5);

  cout<<endl;
  string path_R_prc_fit = plots_dir_R_prc;
  path_R_prc_fit.append(file_name);
  path_R_prc_fit.append("_fit.png");

  const char* R_prc_fit_path = path_R_prc_fit.c_str();

  g_R_prc_E->Draw("AP");
  F2.SaveAs(R_prc_fit_path);
  g_R_prc_E->Delete();

  nf++;

  //-------------------------------------------------------------

  }
  }
  }										//End FILES LOOP

  //-------------------------------------------------------------

  /*	In this part the functions of meanR_rho(E),
	sigmaR_rho(E) and cR_rho(E) are defined.
  */

  TF1 * f_meanR_rho_E = new TF1("f_meanR_rho_E" , " [0]*pow(x, [1]) + [2]");
  f_meanR_rho_E->SetParameters(-3.32045e+07, 0.000957474, 4.77929e+06);

  TF1 * f_sigmaR_rho_E = new TF1("f_sigmaR_rho_E" , " [0]*pow(x, [1]) + [2]");
  f_sigmaR_rho_E->SetParameters(-3.32045e+07, 0.000137065, 3.32486e+07);

  TF1 * f_cR_rho_E = new TF1("f_cR_rho_E" , " [0]*pow(x, [1]) + [2]");
  f_cR_rho_E->SetParameters(-3.32045e+07, 0.000137065, 3.32486e+07);

  //-------------------------------------------------------------

  string plots_dir_R_rho = plots_dir;
  plots_dir_R_rho.append("/R_rho(E)/");

  nf = 0;

  //-------------------------------------------------------------

  /*	Here the graphs of meanR_rho(E), sigmaR_rho(E) and
	cR_rho(E) are loaded from the directory in
	which it were previously saved as .root files.
  */

  if((dir = opendir(plots_dir_R_rho.c_str())) != NULL) {				
  while((ent = readdir(dir)) != NULL) {						//A loop over all files in directory - Start FILES LOOP
  if(opendir(ent->d_name) == NULL) {

								
  cout<< "\nFile name: " <<endl;
  cout<< ent->d_name <<endl;

  string file_graph = plots_dir_R_rho;
  file_graph.append(ent->d_name);						//Adding name of the file into directory path


  if(file_graph.find(".root") > file_graph.size()) continue;			//Checking if its a root file

  cout<< "\nOpening file nr " << nf+1 << " : " <<endl;
  cout<< file_graph <<endl;
  
  TFile *file = new TFile(file_graph.c_str(), "READ");				//Opening each file
  TGraph * g_R_rho_E = (TGraph *) file->Get("Graph");				//Getting a graph from the file

  //-------------------------------------------------------------

  /*	In this part the functions are fitted.
  */


  if(file_graph.find("cR_rho(E)") < file_graph.size()) {

    TFitResultPtr result_R_rho_E = g_R_rho_E->Fit(f_cR_rho_E, "S");  

  }

  if(file_graph.find("sigmaR_rho(E)") < file_graph.size()) {

    TFitResultPtr result_R_rho_E = g_R_rho_E->Fit(f_sigmaR_rho_E, "S");  

  }

  if(file_graph.find("meanR_rho(E)") < file_graph.size()) {

    TFitResultPtr result_R_rho_E = g_R_rho_E->Fit(f_meanR_rho_E, "S");  

  }

  //-------------------------------------------------------------

  TCanvas F3("F"); 
  F3.SetLogy(); 
  F3.SetLogx(); 

  string file_name = ent->d_name;
  int size_name = file_name.size();
  file_name.resize(size_name - 5);

  cout<<endl;
  string path_R_rho_fit = plots_dir_R_rho;
  path_R_rho_fit.append(file_name);
  path_R_rho_fit.append("_fit.png");

  const char* R_rho_fit_path = path_R_rho_fit.c_str();

  g_R_rho_E->Draw("AP");
  F3.SaveAs(R_rho_fit_path);
  g_R_rho_E->Delete();

  nf++;

  //-------------------------------------------------------------

  }
  }
  }										//End FILES LOOP

  //-------------------------------------------------------------

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~				//End MACRO

/* Used functions definitions
*/

string get_current_dir_fit_macro() {						//Function returning path to current directory

   char buff[FILENAME_MAX]; 							//Create string buffer to hold path

   GetCurrentDir( buff, FILENAME_MAX );
   string current_working_dir(buff);

   return current_working_dir;
}


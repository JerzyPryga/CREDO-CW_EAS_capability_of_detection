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

std::string get_current_dir_fit_macro();						//Function returning path to current directory

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/* 	This macro fits functions
	to graphs created later for later use.

	Its starting parameters are:
	- plots_dir - path to directory where the plots are stored,
	- r_number - number of distance points analysed.
*/

void Rho_rNE_fit_macro(string plots_dir = get_current_dir_fit_macro(), int r_number = 10, int n_datas = 18)	//Start MACRO
{

  //-------------------------------------------------------------

  /*	In this part the functions of rho(N)/rho_avr(E)
	for different distances r are defined.
  */

  TF1 * f_Rho_rN = new TF1("f_Rho_rN" , " [0]*pow(x, [1]) + [2]");
  f_Rho_rN->SetParameters(15.5473, 0.91507, 0.0322778);

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

  double Rho_rN_p0_lst[r_number], Rho_rN_p1_lst[r_number], Rho_rN_p2_lst[r_number];	//Arrays of parameters of all functions
  int fit_status = 0;

  cout<< "\nOpening file nr " << nr+1 << " : " <<endl;
  cout<< file_graph <<endl;
  
  TFile *file = new TFile(file_graph.c_str(), "READ");				//Opening each file
  TGraph * g_Rho_rN = (TGraph *) file->Get("Graph");				//Getting a graph from the file

  //-------------------------------------------------------------

  /*	In this part the functions are fitted
	and their parameters are saved into
	a N(E)_params.txt file.
  */

    f_Rho_rN->SetParameters(3.54041e-01, 2.81199e+00, -1.90981e-01);
    if((nr > 0) && (fit_status > 0.0) && (fit_status < 1.0)) f_Rho_rN->SetParameters(Rho_rN_p0_lst[nr-1], Rho_rN_p1_lst[nr-1], Rho_rN_p2_lst[nr-1]);

    TFitResultPtr result_Rho_rN = g_Rho_rN->Fit(f_Rho_rN, "S");  
    Rho_rN_p0_lst[nr] = result_Rho_rN->Parameter(0);
    Rho_rN_p1_lst[nr] = result_Rho_rN->Parameter(1);
    Rho_rN_p2_lst[nr] = result_Rho_rN->Parameter(2);
    fit_status = result_Rho_rN->Error(0);

    param_rN_file << Rho_rN_p0_lst[nr] << " ";
    param_rN_file << Rho_rN_p1_lst[nr] << " ";
    param_rN_file << Rho_rN_p2_lst[nr] << "\n";


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
  TGraph * g_Rho_E = (TGraph *) file->Get("Graph");				//Getting a graph from the file

  //-------------------------------------------------------------

  /*	In this part the functions are fitted.
  */

    f_Rho_E->SetParameters(2.80690e-04, 9.99054e-01, -9.79221e-05);
    if(nr > 0) f_Rho_E->SetParameters(Rho_E_p0_lst[nr-1], Rho_E_p1_lst[nr-1], Rho_E_p2_lst[nr-1]);

    TFitResultPtr result_Rho_E = g_Rho_E->Fit(f_Rho_E, "S");  
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

  TF1 * f_Rho_avr = new TF1("f_Rho_avr" , " [0]*exp([1]*(x + [2])) + [3]*exp([4]*(x + [5]))");
  int n_params = 6;

  //-------------------------------------------------------------

  string plots_dir_Rho_avr = plots_dir;
  plots_dir_Rho_avr.append("/Rho_avr(r)/");

  nf = 0;

  long double** Rho_avr_p_lst = 0;							//Arrays of parameters of all functions		
  Rho_avr_p_lst = new long double*[n_datas];

  for(int i = 0; i < n_datas; i++) {
    Rho_avr_p_lst[i] = new long double[n_params];

    for(int j = 0; j < n_params; j++) {
      Rho_avr_p_lst[i][j] = 0;
    }
  }

  //-------------------------------------------------------------

  /*	Here the histograms of Rho_avr are
	loaded from the directory in
	which it were previously saved as .root files.
  */

  ofstream param_avr_file("Rho(r)_fit_params.txt");

  if((dir = opendir(plots_dir_Rho_avr.c_str())) != NULL) {				
  while((ent = readdir(dir)) != NULL) {						//A loop over all files in directory - Start FILES LOOP
  if(opendir(ent->d_name) == NULL) {

								
  cout<< "\nFile name: " <<endl;
  cout<< ent->d_name <<endl;

  string file_hist = plots_dir_Rho_avr;
  file_hist.append(ent->d_name);						//Adding name of the file into directory path


  if(file_hist.find(".root") > file_hist.size()) continue;			//Checking if its a root file

  cout<< "\nOpening file nr " << nf+1 << " : " <<endl;
  cout<< file_hist <<endl;
  
  TFile * file = new TFile(file_hist.c_str(), "READ");				//Opening each file
  TH1F * h_Rho_avr = (TH1F *) file->Get("h_Rho_avr");				//Getting a histogram from the file

  //-------------------------------------------------------------

  /*	In this part the function is fitted.
  */

    f_Rho_avr->SetParameters(3.09561e-07, -6.68758e-05, 9.88480e+03, 7.22502e-10, -8.41600e-06, 7.99379e+04);
    //if(nf > 0) f_Rho_avr->SetParameters(Rho_avr_p_lst[nf-1][0], Rho_avr_p_lst[nf-1][1], Rho_avr_p_lst[nf-1][2], Rho_avr_p_lst[nf-1][3], Rho_avr_p_lst[nf-1][4], Rho_avr_p_lst[nf-1][5]);

    TFitResultPtr result_rho_avr = h_Rho_avr->Fit(f_Rho_avr, "S");  
    Rho_avr_p_lst[nf][0] = result_rho_avr->Parameter(0);
    Rho_avr_p_lst[nf][1] = result_rho_avr->Parameter(1);
    Rho_avr_p_lst[nf][2] = result_rho_avr->Parameter(2);
    Rho_avr_p_lst[nf][3] = result_rho_avr->Parameter(3);
    Rho_avr_p_lst[nf][4] = result_rho_avr->Parameter(4);
    Rho_avr_p_lst[nf][5] = result_rho_avr->Parameter(5);
    param_avr_file << Rho_avr_p_lst[nf][0] << " ";
    param_avr_file << Rho_avr_p_lst[nf][1] << " ";
    param_avr_file << Rho_avr_p_lst[nf][2] << " ";
    param_avr_file << Rho_avr_p_lst[nf][3] << " ";
    param_avr_file << Rho_avr_p_lst[nf][4] << " ";
    param_avr_file << Rho_avr_p_lst[nf][5] << "\n";

  //-------------------------------------------------------------

  TCanvas F3("F"); 
  F3.SetLogy(); 
  //F3.SetLogx(); 

  string file_name = ent->d_name;
  int size_name = file_name.size();
  file_name.resize(size_name - 5);

  cout<<endl;
  string path_Rho_avr_fit = plots_dir_Rho_avr;
  path_Rho_avr_fit.append(file_name);
  path_Rho_avr_fit.append("_fit.png");

  const char* Rho_avr_fit_path = path_Rho_avr_fit.c_str();

  h_Rho_avr->Draw();
  f_Rho_avr->Draw("same");
  F3.SaveAs(Rho_avr_fit_path);
  h_Rho_avr->Delete();

  nf++;

  //-------------------------------------------------------------

  }
  }
  }										//End FILES LOOP

  param_rE_file.close();

  //-------------------------------------------------------------

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~				//End MACRO

/*	Used functions definitions
*/

std::string get_current_dir_fit_macro() {					//Function returning path to current directory

   char buff[FILENAME_MAX]; 							
   GetCurrentDir( buff, FILENAME_MAX );
   std::string current_working_dir(buff);
   return current_working_dir;

}

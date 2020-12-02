#include <string>
#include <iostream>
#include <vector>
#include <experimental/filesystem>
#include <dirent.h>
#include <codecvt>
#include <chrono> 
using namespace std::chrono; 
namespace fs = std::experimental::filesystem;
using namespace std;

//int lst_dir(const char* path);
void quicksort(double tablica[],int tabid[], int x, int y);

//void fit_Rho_rNE(char *katalog, int n_sim = 1000, int bins_N = 20)
void fit_Rho_rNE(const int n_datas = 18, int prc = 95, int bins_N = 30, int bins_rho = 50, int bins_rN = 40, const int n_min = 0, int r_fit_prime = 0, int r_fit_last = 20)
{

  auto start = high_resolution_clock::now(); 		//Pomiar czasu obliczeń

  //-------------------------------------------------------------

  //vector<string> sim_lst;

  //int ll = lst_dir("/home/jerzy/CREDO/Analiza/new_sim");

  /*string path = "/home/jerzy/CREDO/Analiza/new_sim";
  for(const auto & entry : fs::directory_iterator(path)) {
    cout << entry.path() << endl;
    //sim_list.push_back(entry.path());
    //n_datas++;
  }*/

  const int numb_sim = 18;
  char* casc_lst[numb_sim] = {"/home/jerzy/CREDO/Analiza/More_simulations/p_1_e3.root", "/home/jerzy/CREDO/Analiza/More_simulations/p_1_460_e3.root", "/home/jerzy/CREDO/Analiza/More_simulations/p_2_132_e3.root", "/home/jerzy/CREDO/Analiza/More_simulations/p_3_112_e3.root", "/home/jerzy/CREDO/Analiza/More_simulations/p_4_544_e3.root", "/home/jerzy/CREDO/Analiza/More_simulations/p_6_634_e3.root", "/home/jerzy/CREDO/Analiza/More_simulations/p_1_e4.root", "/home/jerzy/CREDO/Analiza/More_simulations/p_1_4141_e4.root", "/home/jerzy/CREDO/Analiza/More_simulations/p_2_2908_e4.root", "/home/jerzy/CREDO/Analiza/More_simulations/p_3_7111_e4.root", "/home/jerzy/CREDO/Analiza/More_simulations/p_6_012_e4.root", "/home/jerzy/CREDO/Analiza/More_simulations/p_1_e5.root", "/home/jerzy/CREDO/Analiza/More_simulations/p_1_57777_e5.root", "/home/jerzy/CREDO/Analiza/More_simulations/p_2_91888_e5.root", "/home/jerzy/CREDO/Analiza/More_simulations/p_5_39993_e5.root", "/home/jerzy/CREDO/Analiza/More_simulations/p_1_e6.root", "/home/jerzy/CREDO/Analiza/More_simulations/p_1_84813_e6.root", "/home/jerzy/CREDO/Analiza/More_simulations/p_3_69626_e6.root"};
  //double E_lst[n_datas]; 

  double E_lst[n_datas], N_lst[n_datas], sigma_lst[n_datas], c_lst[n_datas];				//Rozkład ilości kaskad
  double E_error[n_datas];										//Rozkład gęstości od E
  double N_mi_lst[bins_rN], rho_N_lst[bins_rN], N_mi_error[bins_rN], rho_N_error[bins_rN];		//Rozkład gęstości od Nmi
  double n_N_mi_counter[bins_rN], n_N_mi[bins_rN], sum_n_N_mi[bins_rN];

  double Rho_N_p0_lst[n_datas], Rho_N_p1_lst[n_datas], Rho_N_p2_lst[n_datas];			//Rozkład gęstości od N - parametry

  const int numb_r = 20;

  long double** rho_E_lst = 0;			//Rozkład gęstości od E,r
  rho_E_lst = new long double*[n_datas];

  for(int i = 0; i < n_datas; i++) {
    rho_E_lst[i] = new long double[numb_r];

    for(int j = 0; j < numb_r; j++) {
      rho_E_lst[i][j] = 0;
    }

  }

  long double** rho_E_error = 0;			//Rozkład gęstości - errors - od E,r
  rho_E_error = new long double*[n_datas];

  for(int i = 0; i < n_datas; i++) {
    rho_E_error[i] = new long double[numb_r];

    for(int j = 0; j < numb_r; j++) {
      rho_E_error[i][j] = 0;
    }

  }


  //-------------------------------------------------------------  

  for(int n = 0; n < n_datas; n++) {
    E_lst[n] = N_lst[n_datas] = sigma_lst[n_datas] = c_lst[n_datas] = 0;
    rho_E_error[n_datas] = 0;
    Rho_N_p0_lst[n] = Rho_N_p1_lst[n] = Rho_N_p2_lst[n] = 0;
  }


  //-------------------------------------------------------------  

  int rho_fit_R_lst[numb_r] = {1000, 2000, 4000, 6000, 8000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 110000, 120000, 130000, 140000, 145000}; 

  //n_datas = sizeof(casc_lst)/sizeof(cacs_lst[0]);


  //-------------------------------------------------------------

  //Rho(r, E) - 2D
  TGraph2D * g_Rho_rE = new TGraph2D();

  double r_rE_max = 150000;
  //double E_rE_max = 4*pow(10, 6);
  //double E_rE_min = 1000;
  double E_rE_max = 4*pow(10, 3);
  double E_rE_min = 1;
  int rE_point = 0;



  //Rho(r, E) - 2D
  //bins_rN = 1000;
  TH2 * h_rho_rE = new TH2F("h_rho_rE", " #rho(r, file) distribution; r [cm]; Number of file ; #rho [cm^{-2}] ", bins_rho, E_rE_min, r_rE_max, n_datas, 0, n_datas);


  //-------------------------------------------------------------  

  for(int nf = n_min; nf < n_datas; nf++) {		//pętla dla każdego pliku


  for(int n = 0; n < bins_rN; n++) {
    N_mi_lst[n] = rho_N_lst[n] = N_mi_error[n] = rho_N_error[n] = 0;
    n_N_mi_counter[n] = n_N_mi[n] = sum_n_N_mi[n] = 0;
  }

  //-------------------------------------------------------------  

  TFile *file = new TFile(casc_lst[nf], "READ");
  TTree *sim = (TTree*) file->Get("sim");
  
  TLeafI* l_id = (TLeafI*)sim->GetLeaf("particle..ParticleID");
  TLeaf* l_x = (TLeaf*)sim->GetLeaf("particle..x");
  TLeaf* l_y = (TLeaf*)sim->GetLeaf("particle..y");
  TLeaf* l_pz = (TLeaf*)sim->GetLeaf("particle..Pz");
  TLeaf* l_px = (TLeaf*)sim->GetLeaf("particle..Px");
  TLeaf* l_py = (TLeaf*)sim->GetLeaf("particle..Py");

  TLeaf* E_sim = (TLeaf*)sim->GetLeaf("shower.Energy");

  sim->GetEntry(0);
  //E_lst[nf] = E_sim->GetValue(0);
  E_lst[nf] = E_sim->GetValue(0)/1000;
  E_error[nf] = 0.05*E_lst[nf];

  //cout<< "\nEnergia cząstki pierwotnej: " << E_lst[nf] << " GeV" <<endl;
  cout<< "\nEnergia cząstki pierwotnej: " << E_lst[nf] << " TeV" <<endl;

  int n_sim = 0;
  while(sim->GetEntry(n_sim) != NULL) {
    n_sim++;
  }

  cout<< "Ilość symulacji w pliku: " << n_sim <<endl;
  cout<< "\n" <<endl;

  long double R_prc_lst[n_sim];			//Rozkład promieni - procenty
  double r_prc_max = 0;

  long double R_rho_lst[n_sim];				//Rozkład primieni - gęstość > tło
  double r_rho_max = 0;


  int tab_N[n_sim];				//Rozkład ilości
  int N_max = 0;


  //-------------------------------------------------------------

  //Rho(r, N)	- 2D
  string title_rho_rN = " #rho_{avr}(r, N_{mi}) distribution for E_{casc} = ";
  title_rho_rN.append(to_string(E_lst[nf]));
  //title_rho_rN.append(" GeV ; r [cm]; N_{mi} ; #rho [particles/cm^{2}]");
  title_rho_rN.append(" TeV ; r [cm]; N_{mi} ; #rho [particles/cm^{2}]");

  const char* rho_rN_title = title_rho_rN.c_str();


  TF1 * f_NE = new TF1("f_NE" , " [0]*pow(x, [1]) + [2] " , 0, E_lst[n_datas]);
  f_NE->SetParameters(15.5244, 0.915374, 0.45552);						//Znaleziona funkcja - N(E)

  TF1 * f_sigmaNE = new TF1("f_sigmaNE" , " [0]*pow(x, [1]) + [2] " , 0, E_lst[n_datas]);
  f_sigmaNE->SetParameters(4.83389, 0.85693, 4.33344);						//Znaleziona funkcja - sigmaN(E)

  TF1 * f_CNE = new TF1("f_CNE" , " [0]*pow(x, [1]) + [2] " , 0, E_lst[n_datas]);
  f_CNE->SetParameters(0.00646735, 0.902696, 105.379);						//Znaleziona funkcja - CN(E)

  double r_rN_max = 150000;

  double N_rN_max = f_NE->Eval(E_lst[nf]) + 4*f_sigmaNE->Eval(E_lst[nf]);
  double N_rN_min = f_NE->Eval(E_lst[nf]) - 4*f_sigmaNE->Eval(E_lst[nf]);
  if(N_rN_min <= 0) N_rN_min = 1;

  TH2 * h_rho_rN = new TH2F("h_rho_rN", rho_rN_title, bins_rho, 0, r_rN_max, bins_rN, N_rN_min, N_rN_max);

  int rho_stop = 0;

  //-------------------------------------------------------------


  //Rho(r, N)	- 1D
  string title_rho_Nr = " #rho_{avr}(r, N_{mi}) distribution for r = ";
  title_rho_Nr.append(to_string(rho_fit_R_lst[r_fit_prime]));
  //title_rho_Nr.append(" GeV ; N_{mi} ; #rho [particles/cm^{2}]");
  title_rho_Nr.append(" TeV ; N_{mi} ; #rho [particles/cm^{2}]");

  const char* rho_Nr_title = title_rho_Nr.c_str();

  //-------------------------------------------------------------

  //Rho(r)
  string title_rho_avr = " #rho_{avr}(r) distribution for E_{casc} = ";
  title_rho_avr.append(to_string(E_lst[nf]));
  //title_rho_avr.append(" GeV ; r [cm] ; #rho [particles/cm^{2}]");
  title_rho_avr.append(" TeV ; r [cm] ; #rho [particles/cm^{2}]");

  const char* rho_avr_title = title_rho_avr.c_str();

  TH1 * h_nr_temp = new TH1F(" ", " ; r [cm]; Number of particles", bins_rho, 0, r_rE_max);
  TH1 * h_nr_sim = new TH1F(" ", " ; r [cm]; Number of particles", bins_rho, 0, r_rE_max);
  TH1 * h_rho_avr = new TH1F(" ", rho_avr_title , bins_rho, 0, r_rE_max);



  //-------------------------------------------------------------

  for(int n = 0; n < n_sim; n++) {			//pętla dla różnych symulacji

  sim->GetEntry(n);
  int len0 = l_pz->GetLen();

  double r_lst_temp[len0];
  for(int i = 0; i < len0; i++) r_lst_temp[i] = 0;
  int ID_lst_temp[len0];
  for(int i = 0; i < len0; i++) ID_lst_temp[i] = 0;

  int ID, N_mi;
  double p_mi, x, y;

  double p_x, p_y, p_z;

  N_mi = 0;

  for(int i = 0; i < len0; i++)  {
    ID = l_id->GetValue(i);
    //p_mi = l_pz->GetValue(i);

    p_x = l_px->GetValue(i);
    p_y = l_py->GetValue(i);    
    p_z = l_pz->GetValue(i);       
    p_mi = sqrt(pow(p_x, 2) + pow(p_y, 2) + pow(p_z, 2));

    //if((ID > 4) && (ID < 7) && (p_mi >= 0.3)) {
    if(ID <= 9) {
      x = l_x->GetValue(i);
      y = l_y->GetValue(i);

      r_lst_temp[N_mi] = sqrt((x*x) + (y*y));
      ID_lst_temp[N_mi] = ID;
      N_mi++;
    }
  }

  tab_N[n] = N_mi;
 
  //-------------------------------------------------------------


  double r_lst[N_mi];
  int ID_lst[N_mi];

  int i_mi = 0;

  for(int i = 0; i < len0; i++)  {
    if(r_lst_temp[i] != 0) {
      r_lst[i_mi] = r_lst_temp[i];
      ID_lst[i_mi] = ID_lst_temp[i];
      i_mi++;
    }
  }

  quicksort(r_lst, ID_lst, 0, N_mi-1);

  double maxR[100]; 
 
  for(int j = 0; j < 100; j++) {
    int sum = 0;
    int n_prc = (j+1)*N_mi/100;
    maxR[j] = 0;

    for(int i = N_mi-1; i > 0; i--) {

      if (sum >= n_prc) break;
      sum = sum + 1;
      maxR[j] = r_lst[i];

    }
  }


  //-------------------------------------------------------------

  for(int ir = 0; ir < N_mi; ir++) {
    h_nr_temp->Fill(r_lst[ir]);
    h_nr_sim->Fill(r_lst[ir]);
  }


  //-------------------------------------------------------------


  int b_i = h_nr_temp->GetXaxis()->FindBin(rho_fit_R_lst[r_fit_prime]);
  int n_i = 0;

  for(int bn = 0; bn < bins_rN; bn++) {
    double n_min = (double(bn)/double(bins_rN))*N_rN_max;
    double n_max = (double(bn+1)/double(bins_rN))*N_rN_max;

    if((tab_N[n] > n_min) && (tab_N[n] <= n_max)) {
      n_i = bn;

      double rho_temp = h_nr_sim->GetBinContent(b_i)/((2*M_PI)*h_nr_sim->GetBinWidth(b_i)*h_nr_sim->GetBinCenter(b_i));

      n_N_mi[n_i] = n_N_mi[n_i] + h_nr_sim->GetBinContent(b_i);
      n_N_mi_counter[n_i]++;

      sum_n_N_mi[n_i] = sum_n_N_mi[n_i] + tab_N[n];
      N_mi_error[n_i] = ((n_max - n_min)/2.0) / f_NE->Eval(E_lst[nf]);

      break;
    }
  }

  for(int b = 0; b < bins_rho; b++) {

    double rho_temp = h_nr_sim->GetBinContent(b)/((2*M_PI)*h_nr_sim->GetBinWidth(b)*h_nr_sim->GetBinCenter(b));
    h_rho_rN->SetBinContent(b, n_i, rho_temp);

    h_nr_sim->SetBinContent(b, 0);

  }


  //-------------------------------------------------------------

  }			//pętla dla różnych symulacji


  //Rho_avr(r)
  int b_i = h_nr_temp->GetXaxis()->FindBin(rho_fit_R_lst[r_fit_prime]);
  for(int b = 0; b < bins_rho; b++) {

    double rho_avr = (1.0/double(n_sim)) * (h_nr_temp->GetBinContent(b)/((2*M_PI)*h_nr_temp->GetBinWidth(b)*h_nr_temp->GetBinCenter(b)));
    h_rho_avr->SetBinContent(b, rho_avr);
    h_rho_rE->SetBinContent(b, nf, rho_avr);

    g_Rho_rE->SetPoint(rE_point, h_nr_temp->GetBinCenter(b), E_lst[nf], rho_avr);
    rE_point++;
  }

  TCanvas A("A");
  A.SetLogy();

  string path_rho_avr = "/home/jerzy/CREDO/Analiza/Fitting/Rho_arv(r)/Rho_avr_";
  path_rho_avr.append(to_string(E_lst[nf]));
  path_rho_avr.append(".png");

  const char* rho_avr_path = path_rho_avr.c_str();

  h_rho_avr->Draw();
  A.SaveAs(rho_avr_path);


  //-------------------------------------------------------------

  //Rho_(E)/Rho_max

  for(int r_fit = r_fit_prime; r_fit < r_fit_last; r_fit++) {					//pętla dla różnych konkretnych R

    int b_temp = h_nr_temp->GetXaxis()->FindBin(rho_fit_R_lst[r_fit]);

    rho_E_lst[nf][r_fit] = h_rho_avr->GetBinContent(b_temp); 
    //rho_E_error[nf] = h_rho_avr->GetBinError(b);
    rho_E_error[nf][r_fit] = 0.1*rho_E_lst[nf][r_fit];

  }				//pętla dla różnych konkretnych R

  //-------------------------------------------------------------

  //Rho_(N)/Rho_avr

  for(int bn = 0; bn < bins_rN; bn++) {

    long double N_mi_ring = 0;
    long double rho_temp = 0;
    if(n_N_mi_counter[bn] > 0) {

      N_mi_ring = (n_N_mi[bn]/n_N_mi_counter[bn]);
      rho_temp = N_mi_ring / ((2*M_PI)*h_nr_temp->GetBinWidth(b_i)*h_nr_temp->GetBinCenter(b_i));
      rho_N_lst[bn] = rho_temp / h_rho_avr->GetBinContent(b_i);

    }

    if((rho_temp == 0) && (bn > 0)) rho_N_lst[bn] = 1.1 * rho_N_lst[bn-1];
    if((rho_temp == 0) && (bn == 0)) rho_N_lst[bn] = 0.1;

    if(n_N_mi_counter[bn] > 0) N_mi_lst[bn] = (sum_n_N_mi[bn]/n_N_mi_counter[bn]) / f_NE->Eval(E_lst[nf]);
    if(n_N_mi_counter[bn] == 0) N_mi_lst[bn] = N_mi_lst[bn-1] + (1.0/double(bins_rN));

    if(n_N_mi_counter[bn] > 0) rho_N_error[bn] = (rho_N_lst[bn] / sqrt(n_N_mi_counter[bn]));
    if((n_N_mi_counter[bn] == 0) && (bn > 0)) rho_N_error[bn] = rho_N_error[bn-1];
    if((n_N_mi_counter[bn] == 0) && (bn == 0)) rho_N_error[bn] = (rho_N_lst[bn] / sqrt(double(n_sim) / double(bins_rN)) );


  }

  h_rho_avr->Delete();
  h_nr_temp->Delete();

  //-------------------------------------------------------------

  //Rho(N_mi, r) - 2D
  TCanvas J("J");
  //J.SetLogx();
  //J.SetLogy();
  J.SetLogz();

  string path_rho_rN = "/home/jerzy/CREDO/Analiza/Fitting/Rho_2D_(r,N)/Rho_avr_rN_";
  path_rho_rN.append(to_string(E_lst[nf]));
  path_rho_rN.append(".png");

  const char* rho_rN_path = path_rho_rN.c_str();

  h_rho_rN->SetMinimum(pow(10, -10));
  h_rho_rN->SetMaximum(pow(10, -3));

  h_rho_rN->Draw("lego");
  J.SaveAs(rho_rN_path);

  h_rho_rN->Delete();


  //-------------------------------------------------------------


  //Rho(N_mi, r) - 1D - funkcja
  TF1 * f_rho_rN = new TF1("f_rho_rN" , " [0]*pow([1], x) + [2]", N_rN_min, N_rN_max);
  if(nf == 0) f_rho_rN->SetParameters(3.54041e-01, 2.81199e+00, -1.90981e-01);
  if(nf > 0) f_rho_rN->SetParameters(Rho_N_p0_lst[nf-1], Rho_N_p1_lst[nf-1], Rho_N_p2_lst[nf-1]);

  //-------------------------------------------------------------


  //Rho(N_mi, r) - 1D
  TCanvas L("L"); 
  L.SetLogy(); 
  //L.SetLogx();  

  cout<<endl;
  for(int l = 0; l < bins_rN; l++) {
    if(l == 0) cout<< "Rho_avr( " << rho_fit_R_lst[r_fit_prime] << ", N ) = ";
    cout<< rho_N_lst[l] << ", ";
  }
  cout<<endl;


  //TGraph * g_Rho_N = new TGraph(bins_rN, N_mi_lst, rho_N_lst);
  TGraphErrors * g_Rho_N = new TGraphErrors(bins_rN, N_mi_lst, rho_N_lst, N_mi_error, rho_N_error);

  cout<< "\n Dopasowanie Rho(N)" <<endl;
  //g_Rho_N->Fit(f_rho_rN);

  TFitResultPtr resultRho_N = g_Rho_N->Fit(f_rho_rN, "S");  
  Rho_N_p0_lst[nf] = resultRho_N->Parameter(0);
  Rho_N_p1_lst[nf] = resultRho_N->Parameter(1);
  Rho_N_p2_lst[nf] = resultRho_N->Parameter(2);


  string title_rhoN_fit = "#rho(N)/#rho_{avr} distribution for r = ";
  title_rhoN_fit.append(to_string(rho_fit_R_lst[r_fit_prime]));
  title_rhoN_fit.append(" cm and E = ");
  title_rhoN_fit.append(to_string(E_lst[nf]));
  title_rhoN_fit.append(" TeV");

  const char* rho_fitN_title = title_rhoN_fit.c_str();

  g_Rho_N->SetTitle(rho_fitN_title);
  g_Rho_N->SetMinimum(pow(10, -2));
  g_Rho_N->SetMaximum(pow(10, 2));

  //g_Rho_N->GetXaxis()->SetTitle("E_{casc} [GeV]");
  g_Rho_N->GetXaxis()->SetTitle("N_{mi}/N_{mi}_{avr}");
  g_Rho_N->GetXaxis()->SetLimits(0.0, 2.0);

  g_Rho_N->GetYaxis()->SetTitle("#rho(N)/#rho_{avr} [particles/cm^{2}]");

  g_Rho_N->SetMarkerStyle(8);

  g_Rho_N->SetMarkerSize(1);



  string name_rhoN_fit = "/home/jerzy/CREDO/Analiza/Fitting/Rho_1D_(r,N)/Rho(N)_";
  name_rhoN_fit.append(to_string(E_lst[nf]));
  name_rhoN_fit.append(".root");

  const char* rho_fitN_name = name_rhoN_fit.c_str();
  g_Rho_N->SaveAs(rho_fitN_name);				//zapisywanie do .root

  string path_rhoN_fit = "/home/jerzy/CREDO/Analiza/Fitting/Rho_1D_(r,N)/Rho(N)_";
  path_rhoN_fit.append(to_string(E_lst[nf]));
  path_rhoN_fit.append(".png");

  const char* rho_fitN_path = path_rhoN_fit.c_str();

  g_Rho_N->Draw("AP");
  L.SaveAs(rho_fitN_path);
  g_Rho_N->Delete();


  //-------------------------------------------------------------



  }		//pętla dla każdego pliku


  //-------------------------------------------------------------



  //Rho(r, E) - r_const
  TCanvas N("N"); 
  N.SetLogy(); 
  N.SetLogx();  

  int E_norm = n_datas - 1;

  cout<<endl;
  for(int r_fit = r_fit_prime; r_fit < r_fit_last; r_fit++) {					//pętla dla różnych konkretnych R
    cout<<endl;
    for(int l = 0; l < n_datas; l++) {
      if(l == 0) cout<< "Rho_avr/Rho_max( " << rho_fit_R_lst[r_fit] << ", E ) = ";
      rho_E_lst[l][r_fit] = rho_E_lst[l][r_fit]/rho_E_lst[E_norm][r_fit];
      cout<< rho_E_lst[l][r_fit] << ", ";
    }
  }					//pętla dla różnych konkretnych R
  cout<<endl;

  for(int r_fit = r_fit_prime; r_fit < r_fit_last; r_fit++) {					//pętla dla różnych konkretnych R

  double rho_E_temp[n_datas], rho_E_error_temp[n_datas];;
  for(int e_temp = 0; e_temp < n_datas; e_temp++) {

    rho_E_temp[e_temp] = rho_E_lst[e_temp][r_fit];
    rho_E_error_temp[e_temp] = rho_E_error[e_temp][r_fit];

  }

  //TGraph * g_Rho_E = new TGraph(n_datas, E_lst, rho_E_lst);
  TGraphErrors * g_Rho_E = new TGraphErrors(n_datas, E_lst, rho_E_temp, E_error, rho_E_error_temp);


  string title_rho_fit = " #rho_{avr}(E) distribution for r = ";
  title_rho_fit.append(to_string(rho_fit_R_lst[r_fit]));
  title_rho_fit.append(" cm");

  const char* rho_fit_title = title_rho_fit.c_str();

  g_Rho_E->SetTitle(rho_fit_title);
  g_Rho_E->SetMinimum(pow(10, -4));
  g_Rho_E->SetMaximum(2);


  //g_Rho_E->GetXaxis()->SetTitle("E_{casc} [GeV]");
  g_Rho_E->GetXaxis()->SetTitle("E_{casc} [TeV]");

  g_Rho_E->GetYaxis()->SetTitle("#rho_{avr}/#rho_max");

  g_Rho_E->SetMarkerStyle(8);

  g_Rho_E->SetMarkerSize(1);


  string name_rho_fit = "/home/jerzy/CREDO/Analiza/Fitting/Rho_fit(E)/Rho(E)_";
  name_rho_fit.append(to_string(rho_fit_R_lst[r_fit]));
  name_rho_fit.append(".root");

  const char* rho_fit_name = name_rho_fit.c_str();
  g_Rho_E->SaveAs(rho_fit_name);				//zapisywanie do .root

  string path_rho_fit = "/home/jerzy/CREDO/Analiza/Fitting/Rho_fit(E)/Rho(E)_";
  path_rho_fit.append(to_string(rho_fit_R_lst[r_fit]));
  path_rho_fit.append(".png");

  const char* rho_fit_path = path_rho_fit.c_str();

  g_Rho_E->Draw("AP");
  N.SaveAs(rho_fit_path);
  g_Rho_E->Delete();

  }					//pętla dla różnych konkretnych R


  //-------------------------------------------------------------



  //Rho(r, E)
  TCanvas O("O");
  //O.SetLogx(); 
  O.SetLogy(); 
  O.SetLogz(); 

  g_Rho_rE->SetTitle("#rho_{avr}(r, E) distribution; r [cm]; E_{casc} [TeV]; #rho_{avr} [particles/cm^{2}]");

  g_Rho_rE->GetXaxis()->SetTitle("r [cm]");
  //g_Rho_rE->GetXaxis()->SetTitle("E_{casc} [GeV]");
  g_Rho_rE->GetYaxis()->SetTitle("E_{casc} [TeV]");
  g_Rho_rE->GetZaxis()->SetTitle("#rho_{avr} [particles/cm^{2}]");

  g_Rho_rE->Draw();
  //g_Rho_rE->Draw("pcol");
  //g_Rho_rE->Draw("colz");
  //g_Rho_rE->Draw("surf1");
  //g_Rho_rE->Draw("lego");

  //g_Rho_rE->GetXaxis()->SetRangeUser(0, r_rE_max);
  //g_Rho_rE->GetYaxis()->SetRangeUser(E_rE_min, E_rE_max);
  g_Rho_rE->SetMinimum(pow(10, -11));
  g_Rho_rE->SetMaximum(pow(10, -3));


  O.Update();
  O.SaveAs("Rho_rE.png");
  g_Rho_rE->Delete();


  //-------------------------------------------------------------


  auto stop = high_resolution_clock::now(); 
  auto duration = duration_cast<microseconds>(stop - start); 

  cout<<"\n" <<endl;
  cout<< "Czas obliczeń: " << duration.count() << " * 10^{-6} s = " <<  duration.count()/pow(10, 6) << " s = "  << duration.count()/(pow(10, 6)*60) << " min" <<endl;		//mierzenie czasu

}

/*int lst_dir(const char* path)
{
  struct dirent* entry;
  DIR* dp;

  dp = opendir(path);
  if(dp == NULL) {
    perror("ERROR");
    return -1;

  }

  while(entry == readdir(dp)) puts(entry->d_name);

  closedir(dp);
  return 0;
}*/

void quicksort(double tablica[], int tabid[], int x, int y)  
{

  int i, j, t, idt;
  double v;
  i = x;
  j = y;
  v = tablica[ x ];
  do {

    while( v < tablica[ i ] ) i++;
       
    while( v > tablica[ j ] ) j--;
       
    if( i <= j ) {

      t = tablica[ i ];
      tablica[ i ] = tablica[ j ];
      tablica[ j ] = t;
      idt = tabid[i];
      tabid[i] = tabid[j];
      tabid[j] = idt;
      i++;
      j--;
    }
  }
  while( i <= j ); {
  
    if( x < j ) quicksort( tablica, tabid, x, j );
   
    if( i < y ) quicksort( tablica, tabid, i, y );
  }
}

#include <string>
using namespace std;

void fit_Rho_E(const int ir = 20)
{

  const int numb_r = 20;

  char* R_graphs_lst[numb_r] = {"/home/jerzy/CREDO/Analiza/Fitting/Rho_fit(E)/Rho(E)_1000.root", "/home/jerzy/CREDO/Analiza/Fitting/Rho_fit(E)/Rho(E)_2000.root", "/home/jerzy/CREDO/Analiza/Fitting/Rho_fit(E)/Rho(E)_4000.root", "/home/jerzy/CREDO/Analiza/Fitting/Rho_fit(E)/Rho(E)_6000.root", "/home/jerzy/CREDO/Analiza/Fitting/Rho_fit(E)/Rho(E)_8000.root", "/home/jerzy/CREDO/Analiza/Fitting/Rho_fit(E)/Rho(E)_10000.root", "/home/jerzy/CREDO/Analiza/Fitting/Rho_fit(E)/Rho(E)_20000.root", "/home/jerzy/CREDO/Analiza/Fitting/Rho_fit(E)/Rho(E)_30000.root", "/home/jerzy/CREDO/Analiza/Fitting/Rho_fit(E)/Rho(E)_40000.root", "/home/jerzy/CREDO/Analiza/Fitting/Rho_fit(E)/Rho(E)_50000.root", "/home/jerzy/CREDO/Analiza/Fitting/Rho_fit(E)/Rho(E)_60000.root", "/home/jerzy/CREDO/Analiza/Fitting/Rho_fit(E)/Rho(E)_70000.root", "/home/jerzy/CREDO/Analiza/Fitting/Rho_fit(E)/Rho(E)_80000.root", "/home/jerzy/CREDO/Analiza/Fitting/Rho_fit(E)/Rho(E)_90000.root", "/home/jerzy/CREDO/Analiza/Fitting/Rho_fit(E)/Rho(E)_100000.root", "/home/jerzy/CREDO/Analiza/Fitting/Rho_fit(E)/Rho(E)_110000.root","/home/jerzy/CREDO/Analiza/Fitting/Rho_fit(E)/Rho(E)_120000.root","/home/jerzy/CREDO/Analiza/Fitting/Rho_fit(E)/Rho(E)_130000.root", "/home/jerzy/CREDO/Analiza/Fitting/Rho_fit(E)/Rho(E)_140000.root", "/home/jerzy/CREDO/Analiza/Fitting/Rho_fit(E)/Rho(E)_145000.root"};


  double R_lst[numb_r] = {1000, 2000, 4000, 6000, 8000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 110000, 120000, 130000, 140000, 145000}; 
  double Rho_E_p0_lst[ir], Rho_E_p1_lst[ir], Rho_E_p2_lst[ir]; 

  //Rho(r, E) - 1D
  TF1 * f_rho_rE = new TF1("f_rho_rE", " [0]*pow(x, [1]) + [2]", 0, 10000);


  for(int nr = 0; nr < ir; nr++) {

    TFile * graph_E = new TFile(R_graphs_lst[nr], "READ");
    TGraphErrors * g_Rho_E = (TGraphErrors *) graph_E->Get("Graph");

    if(nr == 0) f_rho_rE->SetParameters(2.80690e-04, 9.99054e-01, -9.79221e-05);
    if(nr > 0) f_rho_rE->SetParameters(Rho_E_p0_lst[nr-1], Rho_E_p1_lst[nr-1], Rho_E_p2_lst[nr-1]);

    TFitResultPtr resultRho_E = g_Rho_E->Fit(f_rho_rE, "S");  
    Rho_E_p0_lst[nr] = resultRho_E->Parameter(0);
    Rho_E_p1_lst[nr] = resultRho_E->Parameter(1);
    Rho_E_p2_lst[nr] = resultRho_E->Parameter(2);

    TCanvas D("D"); 
    D.SetLogy(); 
    D.SetLogx(); 

    string path_rho_fit = "/home/jerzy/CREDO/Analiza/Fitting/E_rho_fit/Rho(E)_fit_";
    path_rho_fit.append(to_string(R_lst[nr]));
    path_rho_fit.append(".png");

    const char* rho_fit_path = path_rho_fit.c_str();

    g_Rho_E->Draw("AP");
    D.SaveAs(rho_fit_path);
    g_Rho_E->Delete();

  }

  //----------------------------------------------


  TF1 * f_Rho_E_param0 = new TF1("f_Rho_E_param0" , " [0]*pow(x, [1]) + [2] " , 0, 10000);
  f_Rho_E_param0->SetParameters(3.72033e-06, 0.452302, 0.00017608);						//Funkcja param0


  TCanvas A("A"); 
  //A.SetLogy(); 
  //A.SetLogx(); 

  TGraph * g_Rho_E_param0 = new TGraph(ir, R_lst, Rho_E_p0_lst);
  g_Rho_E_param0->Fit(f_Rho_E_param0, "S");

  g_Rho_E_param0->SetTitle("#rho_{p0}(E)");
  //g_Rho_E_param0->SetMinimum(pow(10, -4));
  //g_Rho_E_param0->SetMaximum(2);

  g_Rho_E_param0->GetXaxis()->SetTitle("r [cm]");

  g_Rho_E_param0->GetYaxis()->SetTitle("p0");

  g_Rho_E_param0->SetMarkerStyle(8);

  g_Rho_E_param0->SetMarkerSize(1);


  g_Rho_E_param0->Draw("AP");
  A.SaveAs("/home/jerzy/CREDO/Analiza/Fitting/E_rho_fit/Rho_rE_param1.png");

  //----------------------------------------------

  TF1 * f_Rho_E_param1 = new TF1("f_Rho_E_param1" , " [0] * x + [2] " , 0, 10000);
  f_Rho_E_param1->SetParameters(0.969034, -1.07468e-06);						//Funkcja param1


  TCanvas B("B"); 
  //B.SetLogy(); 
  //B.SetLogx(); 

  TGraph * g_Rho_E_param1 = new TGraph(ir, R_lst, Rho_E_p1_lst);
  g_Rho_E_param1->Fit("pol3", "S");
  //g_Rho_E_param1->Fit(f_Rho_E_param1, "S");

  g_Rho_E_param1->SetTitle("#rho_{p1}(E)");
  //g_Rho_E_param0->SetMinimum(pow(10, -4));
  //g_Rho_E_param0->SetMaximum(2);

  g_Rho_E_param1->GetXaxis()->SetTitle("r [cm]");

  g_Rho_E_param1->GetYaxis()->SetTitle("p1");

  g_Rho_E_param1->SetMarkerStyle(8);

  g_Rho_E_param1->SetMarkerSize(1);

  g_Rho_E_param1->Draw("AP");
  B.SaveAs("/home/jerzy/CREDO/Analiza/Fitting/E_rho_fit/Rho_rE_param2.png");

  //----------------------------------------------

  TF1 * f_Rho_E_param2 = new TF1("f_Rho_E_param2" , " [0]*pow(x, [1]) + [2] " , 0, 10000);
  f_Rho_E_param2->SetParameters(1.90287e-06, 0.434794, -8.08257e-05);						//Funkcja param2


  TCanvas C("C"); 
  //C.SetLogy(); 
  //C.SetLogx(); 

  TGraph * g_Rho_E_param2 = new TGraph(ir, R_lst, Rho_E_p2_lst);
  g_Rho_E_param2->Fit(f_Rho_E_param2, "S");

  g_Rho_E_param2->SetTitle("#rho_{p2}(E)");
  //g_Rho_E_param0->SetMinimum(pow(10, -4));
  //g_Rho_E_param0->SetMaximum(2);

  g_Rho_E_param2->GetXaxis()->SetTitle("r [cm]");

  g_Rho_E_param2->GetYaxis()->SetTitle("p2");

  g_Rho_E_param2->SetMarkerStyle(8);

  g_Rho_E_param2->SetMarkerSize(1);

  g_Rho_E_param2->Draw("AP");
  C.SaveAs("/home/jerzy/CREDO/Analiza/Fitting/E_rho_fit/Rho_rE_param3.png");


}

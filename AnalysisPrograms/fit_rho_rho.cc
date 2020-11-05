#include <string>
using namespace std;

void fit_rho_rho()
{

  TF1 * f_rho_rho1 = new TF1("f_rho_rho1" , " [0]" , 0, 2*pow(10, -6));
  f_rho_rho1->SetParameter(0, 1.1);

  TF1 * f_rho_rho2 = new TF1("f_rho_rho2" , " [0]*pow(x, [1]) " , 2*pow(10, -6), pow(10, -4));
  f_rho_rho2->SetParameters(5.73886, 0.12589);

  TF1 * f_rho_rho3 = new TF1("f_rho_rho3" , " [0]*pow(x, [1]) " , pow(10, -4), 7*pow(10, -4));
  f_rho_rho3->SetParameters(0.174958, -0.25308);

  TF1 * f_rho_rho4 = new TF1("f_rho_rho4" , " [0] " , 7*pow(10, -4), 1.0);
  f_rho_rho4->SetParameter(0, 1.1);

  TCanvas N("N"); 
  //N.SetLogy(); 
  N.SetLogx(); 

  double rho1[2] = {pow(10, -8), 2*pow(10, -6)};
  double rho_rho1[2] = {1.1, 1.1};

  double rho2[2] = {2*pow(10, -6), pow(10, -4)};
  double rho_rho2[2] = {1.1, 1.8};

  double rho3[2] = {pow(10, -4), 7*pow(10, -4)};
  double rho_rho3[2] = {1.8, 1.1};

  double rho4[2] = {7*pow(10, -4), 1.0};
  double rho_rho4[2] = {1.1, 1.1};

  TGraph * g_Rho_E1 = new TGraph(2, rho1, rho_rho1);
  TGraph * g_Rho_E2 = new TGraph(2, rho2, rho_rho2);
  TGraph * g_Rho_E3 = new TGraph(2, rho3, rho_rho3);
  TGraph * g_Rho_E4 = new TGraph(2, rho4, rho_rho4);

  cout<< "\n Dopasowanie #rho_2/#rho_{avr}" <<endl;
  g_Rho_E1->Fit(f_rho_rho1);
  g_Rho_E2->Fit(f_rho_rho2);
  g_Rho_E3->Fit(f_rho_rho3);
  g_Rho_E4->Fit(f_rho_rho4);

  g_Rho_E1->SetTitle("#rho_2/#rho_{avr}(#rho) distribution ");


  g_Rho_E1->GetXaxis()->SetTitle("#rho [particles/m^2]");
  g_Rho_E1->GetXaxis()->SetLimits(pow(10, -8), pow(10, -2));
  g_Rho_E2->GetXaxis()->SetLimits(pow(10, -8), pow(10, -2));
  g_Rho_E3->GetXaxis()->SetLimits(pow(10, -8), pow(10, -2));
  g_Rho_E4->GetXaxis()->SetLimits(pow(10, -8), pow(10, -2));

  g_Rho_E1->GetYaxis()->SetTitle("#rho_2/#rho_{avr}");

  g_Rho_E1->SetMarkerStyle(8);

  g_Rho_E1->SetMarkerSize(1);


  g_Rho_E1->Draw("AP");
  g_Rho_E2->Draw("same");
  g_Rho_E3->Draw("same");
  g_Rho_E4->Draw("same");
  N.SaveAs("/home/jerzy/CREDO/Analiza/fit_2/E_rho_fit/Rho_rho.png");
}

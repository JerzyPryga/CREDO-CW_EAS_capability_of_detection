#include <string>
#include <chrono> 

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#include "TROOT.h"
#include "TCanvas.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussLegendreIntegrator.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TF2.h"
#include "TRandom3.h"

using namespace std;
using namespace std::chrono; 

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

int symbol_newtonaB(int a, int b);
long double QnkB(int n, int k, long double P);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/* 	This macro evaluates number of signals caused
	by the background.

	The meaning of each starting parameters of macro:

	- I_criterium - list of particles IDs,

	- Aprime - surface area of detectors [cm^2] (default = 25),
        - eff - efficiency of the detector [%] (default = 95),
	- n - number of detectors (default = 4),
	- Tdet - detector's registration time [s] (default = 10^(-7)),
	- f_fake - frequency of fake signals from the device [1/s] (default = 0.1),

	- Tp - time of measurement [s] (default = 3600*24*7 = 1 week),
	- fi - maximum vertical angle of arriving particles [degrees] (default = 90).

*/

void BgSignalEstimation(int I_criterium[], double Aprime = 25.0, double eff = 0.95, int n = 4, double Tdet = 2*pow(10, -7), double f_fake = 0.1, int Tp = 3600*24*7, double fi = 90.0)

{
  //-----------------------------------------------------------------

  auto start = high_resolution_clock::now(); 			//Time of computations

  //-----------------------------------------------------------------

  /*	This part of code calculates flux of
	cosmic rays background.
  */

  double phi = (fi/180.0)*M_PI;					//Angle into radians
  double omega = 2*M_PI*(1 - cos(phi));				//Spherical angle

  TF1 f_bg("f_tlo","pow(cos(x), 3)*sin(x)", 0, phi);		//Relationship between effective 

  ROOT::Math::WrappedTF1 wf1(f_bg);
  ROOT::Math::GaussLegendreIntegrator ig;
  ig.SetFunction(wf1);
  ig.SetNumberPoints(1000);
  ig.SetRelTolerance(0.0001);
  double integral1 = ig.Integral(0, phi);				//Integral of previous function

  double I[3] = {70, 32, 2};						//Primary fluxes for different particle types
  string I_names[3] = {"muons", "electrons & positons", "photons"};	//Names of different particle types

  int I_len = sizeof(I)/sizeof(I[0]);					//Lengths of previous arrays			
  int Ic_len = sizeof(I_criterium)/sizeof(I_criterium[0]);

  double Ibg = 0;						//Background flux

  cout<<endl;
  cout<< "INCLUDED PARTICLES TYPES" <<endl;
  cout<<"________________________________________________" <<endl;
  cout<<endl;

  for(int j = 0; j < I_len; j++) {
  bool if_I = false;
	
    for(int jc = 0; jc <= Ic_len; jc++) {			//Checking if this type is in the list
      //cout<< I_criterium[jc] <<endl;
      if(j == I_criterium[jc]) {
        if_I = true;
        break;
      }
    }

    if(if_I == true) {						//Adding each type contribution to the flux
      Ibg = Ibg + I[j]*2*M_PI*integral1;
      cout<< "Including: " << I_names[j] <<endl;
    }
  }

  cout<<endl;
  cout<< "BACKGROUND LEVEL" <<endl;
  cout<<"________________________________________________" <<endl;
  cout<<endl;
  cout<< "Background flux					I_bg = "<< Ibg << " [m^(-2)s^(-1)]"<<endl;
  cout<< "Spherical angle of arriving particles	 	Omega = " << omega << " [rad]" <<endl;

  //-----------------------------------------------------------------

  /*	Here informations about detectors
	are printed and probability of detection
	is calculated.
  */

  long double Pbg;

  Pbg = 1 - exp(- Tdet * eff * ((Aprime/10000) * Ibg + f_fake));

  cout<<endl;
  cout<<"PARAMETERS OF THE SYSTEM" <<endl;
  cout<<"________________________________________________" <<endl;
  cout<<endl;

  cout<<"Number of detectors	 			n = " << n <<endl;
  cout<<"Detector's surface	 			Adet = " << Aprime << " [cm^2]" <<endl;
  cout<<"Detector's registration time 			Tdet = " << Tdet << " [s]" <<endl;
  cout<<"Detector's efficiency 				e_det = " << 100*eff << " %" <<endl;
  cout<<"Fake signals frequency 				f_fake = " << f_fake << " [s^(-1)]" <<endl;
  cout<<"\nProbability of detection in Tdet	 	Pbg = " << Pbg <<endl;
  cout<<"\nTime of measurement	 			Tp = " << Tp << " [s]" <<endl;
  cout<<"						Tp = " << Tp/60 << " [min]" <<endl;
  cout<<"						Tp = " << (Tp/60)/60 << " [h]"  <<endl;
  cout<<"						Tp = " << ((Tp/60)/60)/24 << " [days]" <<endl;
  cout<<"						Tp = " << (((Tp/60)/60)/24)/7 << " [weeks]" <<endl;

  //-----------------------------------------------------------------

  /*	Here some arrays are created.
  */

  long double Qk_back[n];				//An array of coincidance probabilities for Background
  long double NiK_back[n];				//An array of number of coincidance signals

  //-----------------------------------------------------------------

  /*	This part of code computes probability
	of coincidance and expected number
	oc such signals.
  */

  for(int k = (n - 1); k >= 0; k--) {

    Qk_back[k] = QnkB(n, k+1, Pbg);			//Probability of coincidance
    NiK_back[k] = Qk_back[k]*(Tp/Tdet);			//Expected number of signals

  }

  //-----------------------------------------------------------------

  /*	This part just prints the results.
  */

  cout<<endl;
  cout<<"EXPECTED NUMBER OF COINCIDANCE SIGNALS" <<endl;
  cout<<"________________________________________________" <<endl;
  cout<<endl;

  cout<<"k	|	N_bg(k)		|	f_bg(k)" <<endl;
  cout<<"________________________________________________" <<endl;

  for(int i = 0; i < n; i++) {

    cout<< i+1 << "	|	" << NiK_back[i];
    cout<< "	|	" << NiK_back[i]/Tp << " [s^(-1)]";
    cout<< "	|	" << 60*NiK_back[i]/Tp << " [min^(-1)]";
    cout<< "	|	" << 60*60*NiK_back[i]/Tp << " [h^(-1)]";
    cout<< "	|	" << 60*60*24*NiK_back[i]/Tp << " [day^(-1)]";
    cout<<endl;
  }

  //-----------------------------------------------------------------

  cout<<endl;

  //-------------------------------------------------------------

  auto stop = high_resolution_clock::now(); 
  auto duration = duration_cast<microseconds>(stop - start); 			//Stoping counting time of computations

  cout<<"\n" <<endl;
  cout<< "Time of the computations: " << duration.count() << " * 10^{-6} s = " <<  duration.count()/pow(10, 6) << " s = "  << duration.count()/(pow(10, 6)*60) << " min" <<endl;

  //-------------------------------------------------------------

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*	Here used functions are defined
*/

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

int symbol_newtonaB(int a, int b)
{

  return (TMath::Factorial(a)/(TMath::Factorial(b)*TMath::Factorial(a - b)));

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

long double QnkB(int n, int k, long double P)
{
   return (symbol_newtonaB(n, k) * pow(P, k) * pow((1 - P), (n - k)));
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

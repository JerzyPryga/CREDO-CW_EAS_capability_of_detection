using namespace std;

void fit_ro(char* plik_in, int rmax = 150000, int en = 10, int bins = 50, int n_it = 100)

{
  //specify the filename

  double bin_width = rmax/bins;

  TH1 * h_nr = new TH1F(" ", " ; r [cm]; Number of particles", bins, 1000, rmax);
  TH1 * h_ro = new TH1F(" ", " ; r [cm]; #rho [cm^{-2}]", bins, 0, rmax);

  TGraph * gr_rho_r = new TGraph();	//Rysowanie gęstości
  TF1 * fg = new TF1("fg", " [0]*exp([1]*x) + [2]*exp([3]*x) + [4]*exp([5]*x) + [6]*exp([7]*x)", 0, rmax);
  //TF1 * fg = new TF1("fg", " [0]/(pow(x + [2], [3])) ", 0, rmax);

  const int n_params = 8;
  double rho_param_lst[n_params] = {8.97103e-06, -5.48389e-05, 5.59388e-05, -2.54412e-04, 8.20991e-03, -5.46281e-03, -3.88348e+06, -8.43938e+05};

  for(int np = 0; np < n_params; np++) { 

    fg->SetParameter(np, rho_param_lst[np]);

  }

  

  TFile *file = new TFile(plik_in, "READ");
  TTree *tree = (TTree*) file->Get("sim");
  

  TLeafI* l_id = (TLeafI*)tree->GetLeaf("particle..ParticleID");
  TLeaf* l_x = (TLeaf*)tree->GetLeaf("particle..x");
  TLeaf* l_y = (TLeaf*)tree->GetLeaf("particle..y");
  TLeaf* l_px = (TLeaf*)tree->GetLeaf("particle..Px");
  TLeaf* l_py = (TLeaf*)tree->GetLeaf("particle..Py");
  TLeaf* l_pz = (TLeaf*)tree->GetLeaf("particle..Pz");
  TLeaf* l_t = (TLeaf*)tree->GetLeaf("particle..Time");


  int Nl = 0;

  for(int ec = 0; ec < en; ec++) {			//pętla dla różnych symulacji

  tree->GetEntry(ec);

  
  int len0 = l_x->GetLen();
  if((ec % 10) == 0) cout << "len0: " << len0 <<endl;

  //len0 = 600000;

  int ID;
  double tabr[len0];
  int tabID[len0];
  double rr;
  double xx, yy, r;

  double p_mi;
  double p_x, p_y, p_z;

  int len = 0;
  for(int i=0; i<len0; i++)  {
    ID = l_id->GetValue(i);
    //p_mi = l_pz->GetValue(i);

    p_x = l_px->GetValue(i);
    p_y = l_py->GetValue(i);    
    p_z = l_pz->GetValue(i);       
    p_mi = sqrt(pow(p_x, 2) + pow(p_y, 2) + pow(p_z, 2));

    //if((ID > 4) && (ID < 7) && (p_mi >= 0.3)) {
    if(ID <= 9) {
      xx = l_x->GetValue(i);
      yy = l_y->GetValue(i);

      rr = sqrt((xx*xx)+(yy*yy));
      tabr[len] = rr;
      len++;
   
    }   
  }

  Nl = Nl + len;

  if((ec % 10) == 0) {
    cout << " " <<endl;
    cout << "Symulacja nr						: " << ec <<endl;
    cout << "Suma N cząstek (mionów): 					" << Nl <<endl;
    cout << "N cząstek (mionów): 						" << len <<endl;
    cout << " " <<endl;
  }

  for(int ir = 0; ir < len; ir++) {

    h_nr->Fill(tabr[ir]);

  }

  }			//pętla dla różnych symulacji

  for(int rr = 0; rr < bins; rr++) {

    h_ro->SetBinContent(rr, (1/(2*M_PI))*(1/double(en))*h_nr->GetBinContent(rr)/(h_nr->GetBinWidth(rr)*h_nr->GetBinCenter(rr)));

  }


  for(int it = 0; it < n_it; it++) {

    TFitResultPtr resultRho_N = h_ro->Fit(fg, "S"); 

    for(int np = 0; np < n_params; np++) { 

      rho_param_lst[np] = resultRho_N->Parameter(np);

      fg->SetParameter(np, rho_param_lst[np]);

    }
    //h_ro->Fit(fg);

  }

  double suma = 0;

  for(int j = 1; j <= bins; j++) {

    if((j % 10) == 0) {
      cout<<endl;
      cout << "Odległość r			: " << h_nr->GetBinCenter(j) << " [cm]" <<endl;
      cout << "Gestość				: " << h_ro->GetBinContent(j) << " [1/cm^2] =  " << 10000*h_ro->GetBinContent(j) << " [1/m^2]" <<endl;
      cout << "Gestość F(r)			: " << fg->Eval(h_nr->GetBinCenter(j)) << " [1/cm^2] =  " << 10000*fg->Eval(h_nr->GetBinCenter(j)) << " [1/m^2]"  <<endl;
    }
  }

  cout<<endl;

  cout << "Suma N cząstek (mionów): 					" << Nl <<endl;


  TCanvas C("C");
  C.SetLogy();
  //C.SetLogx();

  h_ro->SetLineColor(4);
  h_ro->SetLineStyle(1);
  h_ro->SetLineWidth(2);

  h_ro->Draw();

  C.SaveAs("RoHist(R)2.png");
 
}

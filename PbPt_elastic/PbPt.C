#include <dirent.h>
#include <boost/filesystem.hpp>
//If you don't have boost, the code runs fine
//with std::filesystem instead
//#include <filesystem>
#include "../utils/utils.C"
#include "../utils/cosmetics.C"
double Ebeam = 10.6;
double Pmass = 938.272081/1000.;
  
//-------------------------------------------------------------------------------Define cuts---------------------------------------------------------------------------------------//
TCut exclu_cuts(string target)
{
  //This function reads exclusivity cuts
  //Input: string for target type
  //Ouput: TCut holding the exclusivity cuts
  cout << "Exclusivity cuts" << endl;                                                                                                                                                                    
  TCut cut_all_topo;
  // Open file streams for CD and FD cuts files
  string path_CD = target+"_CD_cuts.txt";
  string path_FD = target+"_FD_cuts.txt";
  fstream file_CD, file_FD;
  file_CD.open(path_CD, ios::in);
  file_FD.open(path_FD, ios::in);
  //Check files sanity
  if (file_CD.bad())
    {
      std::cerr << "Unable to open " + path_CD + "\n";
      exit(8);
    }
  if (file_FD.bad())
    {
      std::cerr << "Unable to open " + path_FD + "\n";
      exit(8);
    }
  cout << ">>> Cut files opened " << path_CD << " and " << path_FD << endl;
  string CD,FD;
  TCut cut_CD, cut_FD;
  // Define topological cuts for CD and FD
  TCut topo_CD = "Nuc_status >= 4000 && Nuc_status < 8000";
  TCut topo_FD = "Nuc_status < 4000 && Nuc_status >=2000";
  //Read exclu cuts
  getline(file_CD, CD);
  cut_CD = CD.c_str();
  file_CD.close();
  getline(file_FD, FD);
  cut_FD = FD.c_str();
  file_FD.close();
  //Build final cut to be applied
  cut_all_topo = (cut_CD && topo_CD) || (cut_FD && topo_FD);
  cut_all_topo =  cut_all_topo + "El_fiducial && Nuc_chi2pid<3.4 && Nuc_chi2pid>-2.6 && TMath::Abs(El_chi2pid)<3";
  cout << ">>> Cuts read" << endl;
  return cut_all_topo;
}

//-----------------------------------------------------------------------------Propagating stat errors-------------------------------------------------------------------//
double error_asym(double Np, double Nm, double dNp, double dNm, double Df=1, double dDf=0)
{
  //This function computes the stat error including the error on the dilution factor A = (Np-Nm)/(Df*(Np+Nm))
  //Input: (double) FCup Normalized yields and their errors, Dilution Factor and its error
  //Ouput: (double) stat error
  double den = Df*pow((Np + Nm),2); 
  double term_Np = 2*Np / den;        
  double term_Nm = 2*Nm / den;
  double asym = (Np-Nm)/(Df*(Np+Nm));
  double term_Df = asym/Df;
  double err_squared = term_Np*term_Np*dNp*dNp + term_Nm*term_Nm*dNm*dNm + term_Df*term_Df*dDf*dDf;
  return TMath::Sqrt(err_squared);
}

//-----------------------------------------------------------------------------Dilution factor-------------------------------------------------------------------//

vector<double> dilution_factor(TTree *REC_C, TTree* REC_signal, TCut exclusivity, TCut Q2_bin, double FCup_run_C, double FCup_run_signal, string target)
{
  //Computing the dilution factor
  //Input: (TTree*) C and Signal data, (TCut) exclusivity and binning, (double) FCup norms
  //Output:(vector<double>) Dilution Factor [0] and its stat error [1]
  //Dilution factor
  //Get Yields
  double yield_C = REC_C->GetEntries(exclusivity+Q2_bin);
  double yield_signal = REC_signal->GetEntries(exclusivity+Q2_bin);
  //Errors
  double dC = TMath::Sqrt(yield_C);
  double dS = TMath::Sqrt(yield_signal);
  //Normalized yields
  double normalised_C =  yield_C/FCup_run_C;
  double normalised_signal =  yield_signal/FCup_run_signal;
  //Dilution Factor
  double constant = 0.7;
  if(target=="NH3") constant = 0.82;
  double Df = 1 - constant*(normalised_C/normalised_signal);
  //Stat Error
  double dDf1 = dC/yield_signal;
  double dDf2 = dS*yield_C/(yield_signal*yield_signal);
  double dDf = constant * (FCup_run_signal/FCup_run_C) * TMath::Sqrt(dDf1*dDf1+dDf2*dDf2);
  cout << ">>> Dilution factor " << Df << " pm "<<  dDf << endl;
  return {Df,dDf};
}

//-----------------------------------------------------------------------------------------MAIN-------------------------------------------------------------------------------------//
int PbPt(string analysis_folder, string filename_C, string filename_signal, string target)
{
  //Main function to compute PbPt
  //Input: (string) filenames of the files containing the list of runs, (string) target type
  //Output: (int) sanity check
  //Results from the analysis are saved in the /results folder
  //It contains: - PbPt plot with the Q^2 dependence of the asymmetries and their ratio
  //             - Dilution factor plot
  //             - Text file with the PbPt result
  
  //Reading the input runs
  vector <int> run_C = get_runs_from_file(filename_C);     
  vector <int> run_signal = get_runs_from_file(filename_signal);      
  //Run range
  int run_min_C = *min_element(run_C.begin(), run_C.end());        
  int run_max_C = *max_element(run_C.begin(), run_C.end());      
  int run_min_signal = *min_element(run_signal.begin(), run_signal.end());  
  int run_max_signal = *max_element(run_signal.begin(), run_signal.end());   
     
  //-------------------------------------------------------------------------------Output plots-------------------------------------------------------------------------------------//
  //The output directory name contains an int that is increased by 1 each time
  //So that previous data is not erased.
  cout << "Creating output dir" << endl;
  
  int i = 0;
  string path = "./results/PbPt_elastic_" + to_string(i) + "/";
  bool condition = DirectoryExists( path.c_str());
  while (condition == true)
    {
      i++;
      path = "./results/PbPt_elastic_" + to_string(i) + "/";
      condition = DirectoryExists( path.c_str());
    }
  
  boost::filesystem::create_directories(path);
  cout << ">>>" << path << endl;
  cout << endl;
  
  //-------------------------------------------------------------------------------Open files and get trees-------------------------------------------------------------------------//
  cout << "Opening input files" << endl;

  vector<int> valid_runs_ND3, valid_runs_C;
  //Data trees for elastic events and fcup norms
  //CARBON
  TChain *REC_C = new TChain("Elastic_loose");
  TChain *scaler_C = new TChain("Scaler info");
  for (int i_files = 0; i_files < run_C.size(); i_files++)        
    {        
      int run_number = run_C.at(i_files);
      string filename_C = analysis_folder + "/pass1_verytight_cut_"+target+"_analysis_elastic_" + to_string(run_number) + ".root";                                                     
      //Check if the file exists        
      if (boost::filesystem::exists(filename_C))  
        {                                  
          REC_C->Add(filename_C.c_str()); 
          scaler_C->Add(filename_C.c_str());
	  valid_runs_C.push_back(run_number);
        }   
      else {cout << filename_C << " does not exist, continuing" << endl;}       
    }
  //SIGNAL
  TChain *REC_signal = new TChain("Elastic_loose");
  TChain *scaler_signal = new TChain("Scaler info");
  for (int i_files = 0; i_files < run_signal.size(); i_files++)       
    {                     
      int run_number = run_signal.at(i_files);
      string filename_signal = analysis_folder + "/pass1_verytight_cut_"+target+"_analysis_elastic_" + to_string(run_number) + ".root";                                                     
      //Check if the file exists
      if (boost::filesystem::exists(filename_signal))  
        {                         
          REC_signal->Add(filename_signal.c_str());    
          scaler_signal->Add(filename_signal.c_str());
	  valid_runs_ND3.push_back(run_number);
        }       
      else {cout << filename_signal << " does not exist, continuing" << endl;}     
    }

  //---------------------------------------------------------------------------Get FCup normalizations-------------------------------------------------------------------------------//
  cout << "FCup normalization" << endl;
  //Adding normalizations between runs
  double FCup_p_signal = 0;
  double FCup_n_signal = 0;
  double FCup_run_signal = 0;
  double FCup_hel_p_signal, FCup_hel_n_signal, FCup_read_run_signal;
  scaler_signal->SetBranchAddress("FCup_hel_p", &FCup_hel_p_signal);
  scaler_signal->SetBranchAddress("FCup_hel_n", &FCup_hel_n_signal);
  scaler_signal->SetBranchAddress("FCup_run", &FCup_read_run_signal);
  
  for (int iEntry = 0; scaler_signal->LoadTree(iEntry) >= 0; ++iEntry)
    {
      scaler_signal->GetEntry(iEntry);
      if(FCup_read_run_signal==0) {cerr << "At least one signal file has 0 FCup normalization." << endl; return 0;}
      FCup_p_signal+=FCup_hel_p_signal;
      FCup_n_signal+=FCup_hel_n_signal;
      FCup_run_signal+=FCup_read_run_signal;
    }
  
  double FCup_p_C = 0;
  double FCup_n_C = 0;
  double FCup_run_C = 0;
  double FCup_hel_p_C, FCup_hel_n_C, FCup_read_run_C;
  scaler_C->SetBranchAddress("FCup_hel_p", &FCup_hel_p_C);
  scaler_C->SetBranchAddress("FCup_hel_n", &FCup_hel_n_C);
  scaler_C->SetBranchAddress("FCup_run", &FCup_read_run_C);
  
  for (int iEntry = 0; scaler_C->LoadTree(iEntry) >= 0; ++iEntry)
    {
      scaler_C->GetEntry(iEntry);
      if(FCup_read_run_C==0) {cerr << "At least one C file has 0 FCup normalization." << endl; return 0;}
      FCup_p_C+=FCup_hel_p_C;
      FCup_n_C+=FCup_hel_n_C;
      FCup_run_C+=FCup_read_run_C;
    }                                                                                                                                                                              
     
   cout << ">>>FCup results" << endl;
   cout << "Target | FC+    | FC-     | FC_run " << endl;
   cout << "Carbon |" << FCup_p_C << " | " << FCup_n_C << " | " << FCup_run_C << endl;
   cout << "Signal |" << FCup_p_signal << " | " << FCup_n_signal << " | " << FCup_run_signal << endl;
   cout << "FC Asymmetry for C runs " << (FCup_p_C - FCup_n_C)/(FCup_p_C + FCup_n_C) << endl;
   cout << "FC Asymmetry for signal runs " << (FCup_p_signal - FCup_n_signal)/(FCup_p_signal + FCup_n_signal) << endl;
   cout << endl;
  
  //-------------------------------------------------------------------------Get exclu cuts------------------------------------------------------------------------------------------//
  TCut total = exclu_cuts(target);
  TCut hel_plus("(Helicity == 1)");
  TCut hel_minus("(Helicity == -1)");
  TCut total_hel_plus  = total + hel_plus ;
  TCut total_hel_minus = total + hel_minus;
    
  //------------------------------------------------------------------------Q2 binning----------------------------------------------------------------------------------------//
  //To be changed with whatever is needed
  //Q2 binning for: - max likelihood estimator
  //                - checking TPol indep of it
  cout << endl;
  cout << "Q2 binning" << endl;
  double Q2_low = 1.8;//1.8;//1.2;
  double Q2_high = 3.5;//3.5;//5;
  int nbins = 10;

  // Create logarithmically spaced bin edges
  double logQ2_low = TMath::Log10(Q2_low);
  double logQ2_high = TMath::Log10(Q2_high);
  double step = (logQ2_high - logQ2_low) / nbins;
    
  cout << ">>> Q2 from " << Q2_low << " GeV^2 to " << Q2_high << " GeV^2 in " << nbins << " bins" << endl;
  cout << endl;
  
  //Binned values for asym computation
  double centers_Q2[nbins], centers_theta[nbins];
  double A_in_bins[nbins];
  double err_Q2[nbins];
  //Binned asym and PbPt are a sanity check
  double N_plus[nbins], N_minus[nbins], asym[nbins], PbPt[nbins];
  double N_plus_C[nbins], N_minus_C[nbins], asym_C[nbins];
  double err_N_plus[nbins], err_N_minus[nbins], err_asym[nbins], err_PbPt[nbins];   
  double err_N_plus_C[nbins], err_N_minus_C[nbins], err_asym_C[nbins];
  double Df_vect[nbins], err_Df[nbins];
  int count_bins = -1;
  double Np_C=0, Nm_C=0; //Integrated yields for C
  
  //Max likelihood binned value
  double numerator = 0;
  double numerator_for_error = 0;
  double denominator = 0;
  double denominator_for_error = 0;
  double numerator_error_DF = 0;
  
  for(double i_log = logQ2_low; TMath::Abs(i_log-logQ2_high) > 0.000001; i_log+=step){    
    
    count_bins++;
    
    double i = TMath::Power(10, i_log);
    double i_plus_step = TMath::Power(10, (i_log+step));

    //Binning cut and display
    TCut Q2_bin = ("Q2>=" + to_string(i) + "&&Q2<" + to_string(i_plus_step)).c_str();
    auto Q2_bin_fortext = ("Q^{2} >= " + to_string(i) + " Q^{2} < " + to_string(i_plus_step));

    cout << "Q2>=" + to_string(i) + "&&Q2<" + to_string(i_plus_step) << endl;

    //Mean Q2 of the bin
    string name_Q2_distrib = ("h_Q2_distrib"+to_string(count_bins));
    auto h_Q2_distrib = new TH1F(name_Q2_distrib.c_str(),"",80,Q2_low,Q2_high);
    string todo_Q2_distrib = "Q2>>" + name_Q2_distrib;
    REC_signal->Draw(todo_Q2_distrib.c_str(),total+Q2_bin,"goff");
    centers_Q2[count_bins] = h_Q2_distrib->GetMean();
    err_Q2[count_bins] = 0;

    //Mean theta_e of the bin
    string name_theta_distrib = ("h_theta_distrib"+to_string(i));         
    auto h_theta_distrib = new TH1F(name_theta_distrib.c_str(),"",80,0,50);                                                                                                                            
    string todo_theta_distrib = "El_Theta>>" + name_theta_distrib;
    REC_signal->Draw(todo_theta_distrib.c_str(),total+Q2_bin,"goff");                                                                                                                             
    centers_theta[count_bins] = h_theta_distrib->GetMean();

    //Theoretical asymmetry. arXiv:0707.1861
    double tau =centers_Q2[count_bins]/(4*Pmass*Pmass);
    double theta_rad = centers_theta[count_bins]*M_PI/180.;
    double tan2 = tan(theta_rad/2)*tan(theta_rad/2);
    double epsilon = 1/(1+2*(1+tau)*tan2);
    double G_E = (1+3.439*tau-1.602*pow(tau,2)+0.068*pow(tau,3))/(1+15.055*tau+48.061*pow(tau,2)+99.304*pow(tau,3)+0.012*pow(tau,4)+8.65*pow(tau,5));
    double G_M = 2.7928473*(1-1.465*tau+1.26*pow(tau,2)+0.262*pow(tau,3))/(1+9.627*tau+11.179*pow(tau,4)+13.245*pow(tau,5));
    double G = G_M/G_E;
    double A_num = 2*tau*G*(Pmass/Ebeam + G*(tau*Pmass/Ebeam + (1+tau)*tan2));
    double A_den = 1 + G*G*tau/epsilon;
    A_in_bins[count_bins] = A_num/A_den;
    //Dilution Factor
    vector<double> Df_and_errors = dilution_factor(REC_C, REC_signal, total, Q2_bin, FCup_run_C, FCup_run_signal, target);
    Df_vect[count_bins] = Df_and_errors[0];
    err_Df[count_bins] = Df_and_errors[1]; 
    
    //Signal
    N_plus[count_bins] = REC_signal->GetEntries(total_hel_plus+Q2_bin); 
    err_N_plus[count_bins] = TMath::Sqrt(N_plus[count_bins]);
    N_minus[count_bins] = REC_signal->GetEntries(total_hel_minus+Q2_bin);
    err_N_minus[count_bins] = TMath::Sqrt(N_minus[count_bins]);
    //FCup normalisation
    N_plus[count_bins] = N_plus[count_bins]/FCup_p_signal;
    err_N_plus[count_bins] = err_N_plus[count_bins]/FCup_p_signal;
    N_minus[count_bins] = N_minus[count_bins]/FCup_n_signal;
    err_N_minus[count_bins] = err_N_minus[count_bins]/FCup_n_signal;

    //Check Carbon Asymmetry
    N_plus_C[count_bins] = REC_C->GetEntries(total_hel_plus+Q2_bin); 
    err_N_plus_C[count_bins] = TMath::Sqrt(N_plus_C[count_bins]);
    Np_C+=N_plus_C[count_bins];
    N_minus_C[count_bins] = REC_C->GetEntries(total_hel_minus+Q2_bin); 
    err_N_minus_C[count_bins] = TMath::Sqrt(N_minus_C[count_bins]);
    Nm_C+=N_minus_C[count_bins];
    N_plus_C[count_bins] = N_plus_C[count_bins]/FCup_p_C;
    err_N_plus_C[count_bins] = err_N_plus_C[count_bins]/FCup_p_C;
    N_minus_C[count_bins] = N_minus_C[count_bins]/FCup_n_C;
    err_N_minus_C[count_bins] = err_N_minus_C[count_bins]/FCup_n_C;
    
    //-------------------------------------------------------------------------Asymmetries----------------------------------------------------------------------------------//
    asym[count_bins] = (N_plus[count_bins] - N_minus[count_bins])/(Df_vect[count_bins]*(N_plus[count_bins] + N_minus[count_bins])); 
    err_asym[count_bins] = error_asym(N_plus[count_bins], N_minus[count_bins], err_N_plus[count_bins], err_N_minus[count_bins], Df_vect[count_bins], err_Df[count_bins]);

    asym_C[count_bins] = (N_plus_C[count_bins] - N_minus_C[count_bins])/((N_plus_C[count_bins] + N_minus_C[count_bins]));
    err_asym_C[count_bins] = error_asym(N_plus_C[count_bins], N_minus_C[count_bins], err_N_plus_C[count_bins], err_N_minus_C[count_bins]);
    
    //Variables needed for the max likelihood method
    double DA = Df_vect[count_bins] * A_in_bins[count_bins];
    numerator += (N_plus[count_bins] - N_minus[count_bins]) * DA;
    denominator += (N_plus[count_bins] + N_minus[count_bins]) * DA * DA;
    numerator_for_error += DA*DA*(N_plus[count_bins]/FCup_p_signal + N_minus[count_bins]/FCup_n_signal);
    numerator_error_DF +=  err_Df[count_bins]*err_Df[count_bins]*A_in_bins[count_bins]*A_in_bins[count_bins]*(N_plus[count_bins] - N_minus[count_bins])*(N_plus[count_bins] - N_minus[count_bins]);
    PbPt[count_bins] = asym[count_bins]/A_in_bins[count_bins];
    err_PbPt[count_bins] = err_asym[count_bins]/A_in_bins[count_bins];
    
    cout << ">>> Asym: " << asym[count_bins] << " pm " <<  err_asym[count_bins] << endl;   
    cout << ">>> PbPt: " << PbPt[count_bins] << " pm " <<  err_PbPt[count_bins] << endl;
    cout << ">>> Asym for Carbon: " << asym_C[count_bins] << " pm " <<  err_asym_C[count_bins] << endl;   
  }

  //Plotting
  style_PbPt_plots();
  TCanvas *c3 = new TCanvas("Asym", "Asym", 1200,1200);
  c3->SetBorderMode(0);
  c3->Divide(1,2);
  c3->cd(1);
  auto gr_asym = new TGraphErrors(nbins,centers_Q2,asym, err_Q2, err_asym);
  gr_asym->SetMarkerStyle(21);
  gr_asym->SetMarkerColor(color_accent);
  gr_asym->SetLineColor(color_accent);
  gr_asym->Draw("ap");
  gr_asym->SetTitle("");
  gr_asym->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
  gr_asym->GetYaxis()->SetTitle("#frac{N^{+} - N^{-}}{N^{+} + N^{-}}");
  gr_asym->GetYaxis()->SetRangeUser(-0.3,0.3);
  gr_asym->GetXaxis()->SetLimits(Q2_low,Q2_high);
  auto gr_asym_th = new TGraph(nbins,centers_Q2,A_in_bins);                                                                                         
  gr_asym_th->SetMarkerStyle(22);
  gr_asym_th->SetMarkerSize(2);
  gr_asym_th->SetMarkerColor(color_accent2);                                                                                                                                                             
  gr_asym_th->SetLineColor(color_accent2);                                                                                                                                     
  gr_asym_th->SetTitle("");
  gr_asym_th->Draw("p && same");
  auto gr_asym_C = new TGraphErrors(nbins,centers_Q2,asym_C, err_Q2, err_asym_C);
  gr_asym_C->SetMarkerStyle(23);
  gr_asym_C->SetMarkerSize(2);
  gr_asym_C->SetMarkerColor(color_C);                                                                                                                                                             
  gr_asym_C->SetLineColor(color_C);                                                                                                                                     
  gr_asym_C->SetTitle("");
  gr_asym_C->Draw("p && same");

  gPad->SetBottomMargin(0.2);
  auto legend = new TLegend(0.12,0.02,0.82,0.07);
  legend->SetFillStyle(0);//Transparent
  legend->SetNColumns(3);
  legend->SetColumnSeparation(0.3);
  legend->AddEntry(gr_asym_th, "Theoretical asymmetry","p");
  legend->AddEntry(gr_asym, "Measured asymmetry","p");
  legend->AddEntry(gr_asym_C, "Carbon asymmetry","p");
  legend->SetTextSize(0.05);
  legend->Draw();
  
  c3->cd(2);
  auto gr_PbPt = new TGraphErrors(nbins,centers_Q2,PbPt, err_Q2, err_PbPt);
  gr_PbPt->SetMarkerStyle(21);
  gr_PbPt->SetMarkerColor(color_accent);
  gr_PbPt->SetLineColor(color_accent);
  gr_PbPt->Draw("ap");
  gr_PbPt->GetYaxis()->SetRangeUser(-1,1);
  TF1 * mean_PbPt = new TF1("mean_PbPt","pol0",Q2_low,Q2_high);
  mean_PbPt->SetLineColor(color_accent);
  mean_PbPt->SetLineStyle(2);
  cout << ">>> Linear fit for mean PbPt" << endl;
  gr_PbPt->Fit("mean_PbPt");
  gr_PbPt->SetTitle("");
  gr_PbPt->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
  gr_PbPt->GetYaxis()->SetTitle("#frac{A}{A_{th}}");
  style_large_plots();
  TCanvas *cdf = new TCanvas("Dilution factor","Dilution factor", 1200,1200);
  cdf->SetBorderMode(0);
  auto gr_df = new TGraphErrors(nbins,centers_Q2,Df_vect,err_Q2, err_Df);
  gr_df->SetMarkerColor(color_C);
  gr_df->SetTitle("");
  gr_df->Draw("ap");
  gr_df->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
  gr_df->GetYaxis()->SetTitle("Dilution factor");
  gr_df->GetYaxis()->SetRangeUser(0,1);
  gr_df->GetXaxis()->SetLimits(Q2_low,Q2_high);
  cdf->Print((path + "Df.pdf").c_str());
  cdf->Print((path + "Df.png").c_str());
  
  //Final PbPt results
  //Asymmetry results for max likelihood and 1 bin integrated methods
  double PbPt_weighted = numerator/denominator;
  double error_DF_2 = numerator_error_DF/(denominator*denominator);
  double error_PbPt_weighted_2 = numerator_for_error/(denominator*denominator);
  double error_PbPt_weighted = TMath::Sqrt(error_PbPt_weighted_2 + error_DF_2);
    
  cout << "Computed PbPt from max likelihood " << PbPt_weighted << " pm " << error_PbPt_weighted << endl;
  cout << "Computed PbPt from the mean of the binned values " << mean_PbPt->GetParameter(0) << " pm " << mean_PbPt->GetParError(0) << endl;

  //Asymmetry results for Carbon, integrated
  double asym_C_integrated = (Np_C/FCup_p_C - Nm_C/FCup_n_C)/(Np_C/FCup_p_C + Nm_C/FCup_n_C);
  double err_asym_C_integrated = error_asym(Np_C/FCup_p_C, Nm_C/FCup_n_C, TMath::Sqrt(Np_C)/FCup_p_C, TMath::Sqrt(Nm_C)/FCup_n_C);   
  cout << "Carbon asymmetry = " << asym_C_integrated << " pm " << err_asym_C_integrated << endl;
  c3->cd(2);
  // Create a text box
  TString text = Form("P_{b} #times P_{t} = %.2f +/- %.2f", PbPt_weighted, error_PbPt_weighted);
  TLatex *textbox = new TLatex(0.6, 0.85, text);
  textbox->SetNDC();
  textbox->SetTextSize(0.05);
  textbox->Draw();
  c3->Print((path + "PbPt.pdf").c_str());
  c3->Print((path + "PbPt.png").c_str());
 
  //Writting results to a txt file
  ofstream outFile(path+"PbPt.txt");
  if (!outFile)
    {
      cerr << "Error: Unable to open the file." << endl;
      return 1;
    }
  outFile << "Signal runs: ";
  for(int i=0; i<valid_runs_ND3.size()-1;i++) outFile << valid_runs_ND3.at(i) << ",";
  outFile << valid_runs_ND3.at(valid_runs_ND3.size()-1) << endl;
  outFile << "C runs: ";
  for(int i=0; i<valid_runs_C.size()-1;i++) outFile << valid_runs_C.at(i) << ",";
  outFile << valid_runs_C.at(valid_runs_C.size()-1) << endl;
  outFile << "PbPt = " << PbPt_weighted << " pm " << error_PbPt_weighted << endl;
  outFile << "Integrated Carbon asymmetry " << asym_C_integrated << " pm " << err_asym_C_integrated << endl;
  outFile << "FC Asymmetry for C runs " << (FCup_p_C - FCup_n_C)/(FCup_p_C + FCup_n_C) << endl;
  outFile << "FC Asymmetry for signal runs " << (FCup_p_signal - FCup_n_signal)/(FCup_p_signal + FCup_n_signal) << endl;
  outFile << "Mean_Q2" << ","<< "Df" << ","<< "dDf" << ","<< "Ath" << ","<< "A_Signal" << ","<< "dA_Signal" << ","<< "A_Carbon" << ","<< "dA_Carbon" << endl;
  for (int i = 0; i < nbins; ++i) {
    outFile << centers_Q2[i] << ","<< Df_vect[i] << ","<< err_Df[i] << ","<< A_in_bins[i]  << ","<< asym[i] << ","<< err_asym[i] << ","<< asym_C[i] << ","<< err_asym_C[i] << endl; 
  }
  outFile.close();
  return 0;
}




#include "Particle.h"
#include "Particle.cpp"
#include "hipo4/reader.h"

int analysis_elastic(int run_number, string target_type)
{
  //Input file
  clas12root::HipoChain hipochain;
  TChain *chain = new TChain("", "");
  string filename = "/cache/clas12/rg-c/production/summer22/pass1/10.5gev/"+target_type+"/dst/train/gmn/gmn_0"+ to_string(run_number) +".hipo"; 
  if(run_number>=16843 && run_number<=17408) filename = "/cache/clas12/rg-c/production/fall22/pass1/"+target_type+"/dst/train/gmn/gmn_0"+ to_string(run_number) +".hipo";
  ifstream my_file(filename);
  if (my_file)
    {
      chain->AddFile(filename.c_str());
      hipochain.Add(filename.c_str());
    } 
  else 
    {
      cerr << "file " + filename + " does not exist" << endl;
      return 0;
    }

  //Output file                                                                                                                                                                                          
  TFile *outfile = new TFile(std::string("/volatile/clas12/pilleux/analysis_results/pass1/elastic/analysis_elastic_"+ to_string(run_number) +  ".root").c_str(), "RECREATE"); 

  //Define hipo variables
  hipo::reader reader;
  hipo::dictionary factory;
  hipo::structure particles;
  hipo::event event;
  
  //Output tree
  TTree *tree = new TTree("Elastic", "Elastic");                                                                                                                                                    
  TTree *tree_scaler = new TTree("Scaler info", "Scaler info");
  int RunNumber, Helicity, Event_ID;
  double El_px, El_py, El_pz;
  double El_E, El_P;
  double El_Theta, El_Phi;
  double El_vx, El_vy,El_vz;
  double El_status, El_chi2pid, El_beta;
  double El_sampling, El_edep_pcal, El_edep_ecin, El_edep_ecout;
  bool El_fiducial;
  int El_sector_pcal;
  double Nuc_px, Nuc_py, Nuc_pz;
  double Nuc_E, Nuc_P;
  double Nuc_Theta, Nuc_Phi;
  double Nuc_vx, Nuc_vy,Nuc_vz;
  double Nuc_status, Nuc_chi2pid, Nuc_beta;
  double mm2_eN_N, Q2, delta_Q2, W2, Emiss,p_perp,p_long;
  double pe_calc,Ee_calc,Ee_calc_purelastic, pp_calc, Q2_calc, dPhi;
  double delta_pe, delta_pp, mm2_total;
  double Ath;
  double FCup_hel_p, FCup_hel_n, FCup_run;
  double nucleon_mass = 938.272081/1000.;
  
  tree->Branch("RunNumber", &RunNumber);
  tree_scaler->Branch("RunNumber", &RunNumber);
  tree->Branch("Helicity", &Helicity);
  tree->Branch("EventID", &Event_ID);
  tree->Branch("El_px", &El_px);
  tree->Branch("El_py", &El_py);
  tree->Branch("El_pz", &El_pz);
  tree->Branch("El_E", &El_E);
  tree->Branch("El_P", &El_P);
  tree->Branch("El_Theta", &El_Theta);
  tree->Branch("El_Phi", &El_Phi);
  tree->Branch("El_vx", &El_vx);
  tree->Branch("El_vy", &El_vy);
  tree->Branch("El_vz", &El_vz);
  tree->Branch("El_status", &El_status);
  tree->Branch("El_chi2pid", &El_chi2pid);
  tree->Branch("El_beta", &El_beta);
  tree->Branch("El_fiducial", &El_fiducial);
  tree->Branch("El_sector_pcal", &El_sector_pcal);
  tree->Branch("El_edep_pcal",&El_edep_pcal);
  tree->Branch("El_edep_ecin",&El_edep_ecin);
  tree->Branch("El_edep_ecout",&El_edep_ecout);
  tree->Branch("El_sampling",&El_sampling);
  tree->Branch("Nuc_px", &Nuc_px);
  tree->Branch("Nuc_py", &Nuc_py);
  tree->Branch("Nuc_pz", &Nuc_pz);
  tree->Branch("Nuc_E", &Nuc_E);
  tree->Branch("Nuc_P", &Nuc_P);
  tree->Branch("Nuc_Theta", &Nuc_Theta);
  tree->Branch("Nuc_Phi", &Nuc_Phi);
  tree->Branch("Nuc_vx", &Nuc_vx);
  tree->Branch("Nuc_vy", &Nuc_vy);
  tree->Branch("Nuc_vz", &Nuc_vz);
  tree->Branch("Nuc_status", &Nuc_status);
  tree->Branch("Nuc_chi2pid", &Nuc_chi2pid);
  tree->Branch("Nuc_beta", &Nuc_beta);
  tree->Branch("mm2_eN_N", &mm2_eN_N);
  tree->Branch("mm2_total", &mm2_total);
  tree->Branch("Q2", &Q2);
  tree->Branch("delta_Q2", &delta_Q2);
  tree->Branch("W2", &W2); 
  tree->Branch("Emiss", &Emiss);
  tree->Branch("p_perp", &p_perp);
  tree->Branch("p_long", &p_long); 
  tree->Branch("pe_calc", &pe_calc);
  tree->Branch("Ee_calc", &Ee_calc);
  tree->Branch("Ee_calc_purelastic", &Ee_calc_purelastic);
  tree->Branch("pp_calc", &pp_calc);
  tree->Branch("Q2_calc", &Q2_calc);
  tree->Branch("dPhi", &dPhi);
  tree->Branch("delta_pe", &delta_pe);
  tree->Branch("delta_pp", &delta_pp);
  tree->Branch("Ath", &Ath);
  tree_scaler->Branch("FCup_hel_p", &FCup_hel_p);
  tree_scaler->Branch("FCup_hel_n", &FCup_hel_n);
  tree_scaler->Branch("FCup_run", &FCup_run);

  //Incident electron beam                                                                                                                                                    
  double Ebeam = 10.5473;
  TVector3 beam_energy(0,0,Ebeam);
  Particle electron_beam(beam_energy, 11);
  electron_beam.construct_all_properties();
  
  //Proton target                                                                                                                                                                                        
  TVector3 target(0,0,0);
  Particle proton_target(target, 2212);
  proton_target.construct_all_properties();

  //Reconstructed particles
  vector<Particle> electrons;
  vector<Particle> protons;
  
  //Handling hipo files
  reader.open(filename.c_str());
  
  reader.readDictionary(factory);
  hipo::bank CONF(factory.getSchema("RUN::config"));
  hipo::bank HEL(factory.getSchema("REC::Event"));
  hipo::bank PART(factory.getSchema("REC::Particle"));
  hipo::bank HEL_SCALER(factory.getSchema("HEL::scaler"));
  hipo::bank RUN_SCALER(factory.getSchema("RUN::scaler"));
  hipo::bank CALO(factory.getSchema("REC::Calorimeter"));
  hipo::bank TRAJ(factory.getSchema("REC::Traj"));

  //FCup counts from scalers
  FCup_hel_n = 0;
  FCup_hel_p = 0;      
  FCup_run = 0;

  //Start loop on events
  while (reader.next() == true)
    {
      
      //Get hipo info
      reader.read(event);  
      event.getStructure(PART);
      event.getStructure(HEL);
      event.getStructure(CONF);
      event.getStructure(HEL_SCALER);
      event.getStructure(RUN_SCALER);
      event.getStructure(TRAJ);
      event.getStructure(CALO);
      //Event info                                                                                                                                                                                   
      RunNumber = CONF.getInt("run", 0);
      Helicity = HEL.getInt("helicity", 0);
      
      //Scalers information
      for(int row=0; row<HEL_SCALER.getRows(); row++)
	{
	  int hel = HEL_SCALER.getByte("helicity",row);
	  if(hel == -1) FCup_hel_n+=HEL_SCALER.getFloat("fcupgated",row);
	  if(hel == 1) FCup_hel_p+=HEL_SCALER.getFloat("fcupgated",row);
	}

      if(RUN_SCALER.getRows()>0)
	{
	  float entry = RUN_SCALER.getFloat("fcupgated",0);
	  if(entry>FCup_run) FCup_run=entry;
	  //FCup_run=RUN_SCALER.getFloat("fcupgated",0);
	}
      
      //Number of particles in each event
      int size = PART.getRows();	  
      
      //Start loop on REC particles
      for (int i = 0; i < PART.getRows(); i++)
	{
	  //Get particle info
	  int    pid = PART.getInt("pid", i);
	  if (pid!= 11 && pid!=2212) continue;
	  double px  = PART.getFloat("px", i);
	  double py  = PART.getFloat("py", i);
	  double pz  = PART.getFloat("pz", i);
	  TVector3 momentum(px,py,pz);
	  double vx  = PART.getFloat("vx", i);
	  double vy  = PART.getFloat("vy", i);
	  double vz  = PART.getFloat("vz", i);
	  TVector3 position(vx,vy,vz);
	  int status = PART.getInt("status", i);
	  double chi2pid = PART.getFloat("chi2pid", i);
	  int charge = PART.getInt("charge", i);      	    
	  double beta=PART.getFloat("beta",i);
	  
	  //Electron
	  if(pid==11)
	    {
	      Particle electron(momentum, position, pid, chi2pid, status, charge);
	      electron.construct_all_properties();
	      
	      TLorentzVector virtual_ph = electron_beam.LorentzVector- electron.LorentzVector;
	      Q2 = - virtual_ph.Mag2();
	      W2 = (proton_target.LorentzVector + virtual_ph).Mag2();
	      mm2_eN_N = (electron_beam.LorentzVector + proton_target.LorentzVector - electron.LorentzVector).M2();
	      El_px = px;
	      El_py = py;
	      El_pz = pz;
	      El_E = electron.energy;
	      El_P = electron.getmomentum().Mag();
	      El_Theta = electron.theta;
	      El_Phi = electron.phi;
	      El_vx = vx;
	      El_vy = vy;
	      El_vz = vz;
	      El_status = status;
	      El_chi2pid = chi2pid;
	      El_beta = beta;
	      
	      pe_calc = Ebeam / (1 + (Ebeam/nucleon_mass)*(1-cos(El_Theta*TMath::Pi()/180.)));
	      Ee_calc_purelastic = Ebeam / (1 + 2*(Ebeam/nucleon_mass) * pow(sin((El_Theta*TMath::Pi()/180.)/2.),2));
	      pp_calc = pe_calc * TMath::Sqrt(pow((1+Ebeam/nucleon_mass),2)*pow((1-cos(El_Theta*TMath::Pi()/180.)),2)+pow(sin(El_Theta*TMath::Pi()/180.),2));
	      Q2_calc = 4*Ebeam*El_E*pow(sin(El_Theta*TMath::Pi()/180.)/2.,2);		  
	      
	      delta_pe = pe_calc - El_P;
	      delta_Q2 = Q2_calc - Q2;

	      double sampling=0, edep_pcal=0, edep_ecin=0, edep_ecout=0, edep=0;
	      int sector=0;
	      bool fiducial=true;

	      //fiducial DC                                                                                                                                                                            
	      for (int i_traj = 0; i_traj < TRAJ.getRows(); i_traj++)
		{
		  int pindex = TRAJ.getInt("pindex", i_traj);
		  if(pindex!=i) continue; //looking at current particle			  
		  int detector = TRAJ.getInt("detector" , i_traj);
		  if(detector!=6) continue;//LOOKING AT DRIFT CHAMBERS                                                                                                                                 
		  if(TRAJ.getFloat("edge" ,i_traj)<4) fiducial=false;
		}

	      
	      for (int i_calo = 0; i_calo < CALO.getRows(); i_calo++)
		{
		  if(CALO.getInt("pindex", i_calo)!=i) continue;//Only current particle                                                                                                                
		  if(CALO.getInt("detector", i_calo)!=7) continue;//ECAL                                                                                                                               
		  double edep_calo = CALO.getFloat("energy", i_calo);
		  if(CALO.getInt("layer", i_calo)==4) {edep_ecin = edep_calo; edep+=edep_calo;}//ECIN                                                                                                  
		  if(CALO.getInt("layer", i_calo)==7) {edep_ecout = edep_calo; edep+=edep_calo;}//ECOUT                                                                                                
		  if(CALO.getInt("layer", i_calo)==1)
		    {//PCAL                                                                                                                                                                            
		      edep_pcal = edep_calo;
		      sector = CALO.getInt("sector", i_calo);
		      edep+=edep_calo;
		      if(CALO.getFloat("lv", i_calo)<8 || CALO.getFloat("lw", i_calo)<8) fiducial=false;
		    }
		}
	      sampling = edep/momentum.Mag();
	      El_edep_pcal = edep_pcal;
	      El_edep_ecin = edep_ecin;
	      El_edep_ecout = edep_ecout;
	      El_sampling=sampling;
	      El_fiducial=fiducial;
	      El_sector_pcal=sector;
	      electrons.push_back(electron);
	    }
	  
	  //Proton
	  if(pid==2212)
	    {
	      Particle proton(momentum, position, pid, chi2pid, status, charge);
	      proton.construct_all_properties();
	      Nuc_px = px;
	      Nuc_py = py;
	      Nuc_pz = pz;
	      Nuc_E = proton.energy;
	      Nuc_P = proton.getmomentum().Mag();
	      Nuc_Theta = proton.theta;
	      Nuc_Phi = proton.phi;
	      Nuc_vx = vx;
	      Nuc_vy = vy;
	      Nuc_vz = vz;
	      Nuc_status = status;
	      Nuc_chi2pid = chi2pid;
	      Nuc_beta=beta;
	      
	      protons.push_back(proton);
	    }
	}//END LOOP ON PARTICLES
      
      if(protons.size() == 1 && electrons.size()==1)
	{
	  Ee_calc = nucleon_mass*(1/(tan(electrons[0].theta*(TMath::Pi()/180.)/2)*tan(protons[0].theta*TMath::Pi()/180.))-1.);
	  dPhi = electrons[0].phi - protons[0].phi;
	  delta_pp = pp_calc - protons[0].getmomentum().Mag();   
	  Emiss = electron_beam.energy + proton_target.energy - electrons[0].energy - protons[0].energy;
	  mm2_total = (electron_beam.LorentzVector + proton_target.LorentzVector - electrons[0].LorentzVector - protons[0].LorentzVector).M2();
	  p_long = electron_beam.getmomentum().z() + proton_target.getmomentum().z() - electrons[0].getmomentum().z()- protons[0].getmomentum().z();
	  double missing_px = electron_beam.getmomentum().x() + proton_target.getmomentum().x() - electrons[0].getmomentum().x()- protons[0].getmomentum().x();
	  double missing_py = electron_beam.getmomentum().y() + proton_target.getmomentum().y() - electrons[0].getmomentum().y()- protons[0].getmomentum().y();
	  p_perp = sqrt(missing_px*missing_px + missing_py*missing_py);
	  
	  double tau = Q2/(4*nucleon_mass*nucleon_mass);
	  double Ge = (1+3.439*tau-1.602*pow(tau,2)+0.068*pow(tau,3))/(1+15.055*tau+48.061*pow(tau,2)+99.304*pow(tau,3)+0.012*pow(tau,4)+8.65*pow(tau,5));
	  double Gm = 2.7928473*(1-1.465*tau+1.26*pow(tau,2)+0.262*pow(tau,3))/(1+9.627*tau+11.179*pow(tau,4)+13.245*pow(tau,5));
	  double G = Gm/Ge;
	  double ME = nucleon_mass/electron_beam.energy;
	  double tanterm = (1+tau)*pow(tan(El_Theta*M_PI/(180*2.)),2);
	  double num_Ath = 2*tau*G*(ME + G*(tau*ME + tanterm));
	  double epsilon = 1./(1+2*tanterm);
	  double den_Ath = 1+G*G*tau/epsilon;
	  Ath = num_Ath/den_Ath;
	  tree->Fill();
	}
      
      electrons.clear();
      protons.clear();
      
    }//END LOOP ON EVENTS
  
  tree_scaler->Fill();
  outfile->cd(); 
  tree->Write();
  tree_scaler->Write();
  outfile->Close();
  return 1;
}

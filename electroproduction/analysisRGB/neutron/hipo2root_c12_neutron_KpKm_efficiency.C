#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
//#include <HipoChain.h>
#include <filesystem> // Pour parcourir un dossier
#include <vector>
#include <algorithm>  // pour std::min_element
#include <iostream>

namespace fs = std::filesystem;

using namespace clas12;
using namespace clas12root;
using namespace QA;


void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	      rp->par()->getPz(),p4.M());

}

void hipo2root_c12_neutron_KpKm_efficiency(){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();

  bool isMC = true;

  double Lum;
  double Q;
  double f2; //facteur correctif si le nombre de simu qui a run n'est pas 1000
  double N_gen = 2000000;

  bool fall2019 = true;
  bool spring2020 = false;
  bool spring2019 = false;

  if(fall2019 == true){

      f2 = 1.28; // facteur car toute les simu n'ont pas reussi sur osg
      Q = 1.2157*1e7;
      Lum = 1.51434*Q*f2;

  }

  if(spring2020 == true){

       Q = 2.83024*1e7;
       Lum = 2.83024*1e7;
        
  }

  if(spring2019 == true){

       Lum = 0;
        
  }

  /////////////////////////////////////
  //LECTURE FICHIER HIPO

  HipoChain chain;

  // Dossier contenant tes fichiers
  TString folder = "/lustre24/expphy/volatile/clas12/osg/mronay/10255"; //

  ///lustre24/expphy/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/nSidis

  int fileCount = 0;
  const int maxFiles = 1200;

  // Parcourir tous les fichiers du dossier
  for (const auto &entry : fs::directory_iterator(folder.Data())) {
      TString filename = entry.path().c_str();
      if (filename.EndsWith(".hipo")) {
          chain.Add(filename.Data());
          fileCount++;
          std::cout << "Ajout du fichier : " << filename << std::endl;

          if (fileCount >= maxFiles) {
            std::cout << "Limite de " << maxFiles << " fichiers atteinte, arrêt du chargement." << std::endl;
            break;
          }
          
      }
  }

  auto config_c12 = chain.GetC12Reader();
  //config_c12->db()->turnOffQADB();
  auto &c12 = chain.C12ref();
  

  //uint EVENTD = config_c12->addBank("DECAYS::Event");
  //uint PARTD = config_c12->addBank("DECAYS::Particle");

  uint EVENTMC = config_c12->addBank("MC::Event");
  uint PartMC = config_c12->addBank("MC::Particle");



  /////CREER LE FICHIER ROOT ET LES TREEs et variables associées

  cout<<"Variables initi ..."<<endl;

  TFile *outFile = new TFile("fall2019_outbending_neutronKpKm_mc_efficiency.root", "recreate");   // nom a choisir rgA_outbendingfall2018_data.root or outbending_mc.root for exemple


  TTree *outT = new TTree("tree", "tree");


  TTree *outT_gen = new TTree("tree_gen", "tree_gen");

  //VARIABLE REC (MC ou DATA)
  
  auto db=TDatabasePDG::Instance();
  double Eb = 10.410;
  double m_target = db->GetParticle(2112)->Mass(); // cas e n -> e' n' K+ K-


  TLorentzVector Missing_nucleon(0,0,0,0);
  TLorentzVector beam(0,0, Eb, Eb);
  TLorentzVector target(0,0,0,m_target); //cible neutron
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector neutron(0,0,0,m_target); //neutron
  TLorentzVector Kp(0,0,0,db->GetParticle(321)->Mass());
  TLorentzVector Km(0,0,0,db->GetParticle(-321)->Mass());
  TLorentzVector q(0, 0, 0, 0);
  TLorentzVector miss(0, 0, 0, 0);


  double Q2, t, W;

  double pid;

  double MinvKpKm, MM;

  double status_Kp, status_Km, status_n, status_el;

  double e_vx, e_vy, e_vz;
  double Kp_vx, Kp_vy, Kp_vz;
  double Km_vx, Km_vy, Km_vz;
  double n_vx, n_vy, n_vz;

  double E_PCAL_n;
  double Lu_PCAL_n, Lv_PCAL_n, Lw_PCAL_n;
  double x_PCAL_n, y_PCAL_n;

  double E_ECIN_n;
  double Lu_ECIN_n, Lv_ECIN_n, Lw_ECIN_n;
  double x_ECIN_n, y_ECIN_n;

  double E_ECOUT_n;
  double Lu_ECOUT_n, Lv_ECOUT_n, Lw_ECOUT_n;
  double x_ECOUT_n, y_ECOUT_n;

  double weight;
  double real_weight;

  double t_missing_nucleon;
  double angle_neutron_missnucl;


  //VARIABLES GEN

  TLorentzVector el_gen(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector neutron_gen(0,0,0,m_target); //neutron
  TLorentzVector Kp_gen(0,0,0,db->GetParticle(321)->Mass());
  TLorentzVector Km_gen(0,0,0,db->GetParticle(-321)->Mass());
  TLorentzVector q_gen(0, 0, 0, 0);
  TLorentzVector miss_gen(0, 0, 0, 0);

  double Q2_gen, t_gen, W_gen;

  double MinvKpKm_gen, MM_gen;

  double e_vz_gen, Kp_vz_gen, Km_vz_gen, n_vz_gen;

  double weight_gen, real_weight_gen;




  std::vector<double> L;


	outT->Branch("Electron", "TLorentzVector", &el);
	outT->Branch("Neutron", "TLorentzVector", &neutron);
  outT->Branch("Kp", "TLorentzVector", &Kp);
  outT->Branch("Km", "TLorentzVector", &Km);
  outT->Branch("Missing", "TLorentzVector", &miss);
  outT->Branch("Missing_nucleon", "TLorentzVector" , &Missing_nucleon);

  outT->Branch("Q2", &Q2, "Q2/D");
  outT->Branch("t", &t, "t/D");
  outT->Branch("t_missing_nucleon", &t_missing_nucleon, "t_missing_nucleon/D");
  outT->Branch("W", &W, "W/D");

  outT->Branch("angle_neutron_missnucl", &angle_neutron_missnucl, "angle_neutron_missnucl/D");

  outT->Branch("e_vx", &e_vx, "e_vx/D");
  outT->Branch("e_vy", &e_vy, "e_vy/D");
  outT->Branch("e_vz", &e_vz, "e_vz/D");

  outT->Branch("Kp_vx", &Kp_vx, "Kp_vx/D");
  outT->Branch("Kp_vy", &Kp_vy, "Kp_vy/D");
  outT->Branch("Kp_vz", &Kp_vz, "Kp_vz/D");

  outT->Branch("Km_vx", &Km_vx, "Km_vx/D");
  outT->Branch("Km_vy", &Km_vy, "Km_vy/D");
  outT->Branch("Km_vz", &Km_vz, "Km_vz/D");

  outT->Branch("n_vx", &n_vx, "n_vx/D");
  outT->Branch("n_vy", &n_vy, "n_vy/D");
  outT->Branch("n_vz", &n_vz, "n_vz/D");

  outT->Branch("MinvKpKm", &MinvKpKm, "MinvKpKm/D");
  outT->Branch("MissingMass", &MM, "MM/D");

  outT->Branch("Status_Kp", &status_Kp, "status_Kp/D");
  outT->Branch("Status_Km", &status_Km, "status_Km/D");
  outT->Branch("Status_n", &status_n, "status_n/D");
  outT->Branch("Status_el", &status_el, "status_el/D");

  //calo neutron
  outT->Branch("E_PCAL_n", &E_PCAL_n, "E_PCAL_n/D");
  outT->Branch("Lu_PCAL_n", &Lu_PCAL_n, "Lu_PCAL_n/D");
  outT->Branch("Lv_PCAL_n", &Lv_PCAL_n, "Lv_PCAL_n/D");
  outT->Branch("Lw_PCAL_n", &Lw_PCAL_n, "Lw_PCAL_n/D");
  outT->Branch("x_PCAL_n", &x_PCAL_n, "x_PCAL_n/D");
  outT->Branch("y_PCAL_n", &y_PCAL_n, "y_PCAL_n/D");

  outT->Branch("E_ECIN_n", &E_ECIN_n, "E_ECIN_n/D");
  outT->Branch("Lu_ECIN_n", &Lu_ECIN_n, "Lu_ECIN_n/D");
  outT->Branch("Lv_ECIN_n", &Lv_ECIN_n, "Lv_ECIN_n/D");
  outT->Branch("Lw_ECIN_n", &Lw_ECIN_n, "Lw_ECIN_n/D");
  outT->Branch("x_ECIN_n", &x_ECIN_n, "x_ECIN_n/D");
  outT->Branch("y_ECIN_n", &y_ECIN_n, "y_ECIN_n/D");

  outT->Branch("E_ECOUT_n", &E_ECOUT_n, "E_ECOUT_n/D");
  outT->Branch("Lu_ECOUT_n", &Lu_ECOUT_n, "Lu_ECOUT_n/D");
  outT->Branch("Lv_ECOUT_n", &Lv_ECOUT_n, "Lv_ECOUT_n/D");
  outT->Branch("Lw_ECOUT_n", &Lw_ECOUT_n, "Lw_ECOUT_n/D");
  outT->Branch("x_ECOUT_n", &x_ECOUT_n, "x_ECOUT_n/D");
  outT->Branch("y_ECOUT_n", &y_ECOUT_n, "y_ECOUT_n/D");


  

  if (isMC){

    outT_gen->Branch("Electron_gen", "TLorentzVector", &el_gen);
	  outT_gen->Branch("Neutron_gen", "TLorentzVector", &neutron_gen);
    outT_gen->Branch("Kp_gen", "TLorentzVector", &Kp_gen);
    outT_gen->Branch("Km_gen", "TLorentzVector", &Km_gen);
    outT_gen->Branch("Missing_gen", "TLorentzVector", &miss_gen);
    

    outT->Branch("weight", &weight, "weight/D");
    outT->Branch("real_weight", &real_weight, "real_weight/D");

    outT_gen->Branch("Q2_gen", &Q2_gen, "Q_gen/D");
    outT_gen->Branch("t_gen", &t_gen, "t_gen/D");
    outT_gen->Branch("W_gen", &W_gen, "W_gen/D");

    outT_gen->Branch("MinvKpKm_gen", &MinvKpKm_gen, "MinvKpKm_gen/D");
    outT_gen->Branch("MissingMass_gen", &MM_gen, "MM_gen/D");

    outT_gen->Branch("weight_gen", &weight_gen, "weight_gen/D");
    outT_gen->Branch("real_weight_gen", &real_weight_gen, "real_weight_gen/D");

    outT_gen->Branch("e_vz_gen", &e_vz_gen, "e_vz_gen/D");
    outT_gen->Branch("Kp_vz_gen", &Kp_vz_gen, "Kp_vz_gen/D");
    outT_gen->Branch("Km_vz_gen", &Km_vz_gen, "Km_vz_gen/D");
    outT_gen->Branch("n_vz_gen", &n_vz_gen, "n_vz_gen/D");


  }

  
  gBenchmark->Start("timer");
  int counter=0;

  //QADB

  double total_charge = 0.0;
  long total_good_events = 0;
  int counter2 = 0;
  double charge_tempo = 0;

  QADB *qa;
  qa = new QADB("latest");
  qa->CheckForDefect("TotalOutlier", true);
  qa->CheckForDefect("TerminalOutlier", true);
  qa->CheckForDefect("MarginalOutlier", true);
  qa->CheckForDefect("SectorLoss", true);
  qa->CheckForDefect("LowLiveTime", true);
  qa->CheckForDefect("ChargeHigh", true);
  qa->CheckForDefect("ChargeNegative", true);
  qa->CheckForDefect("ChargeUnknown", true);
  qa->CheckForDefect("PossiblyNoBeam", true);
  qa->CheckForDefect("Misc", true);
  std::vector<int> runs_allow_misc = {5046, 5047, 5051, 5128, 5129, 5130, 5158, 5159, 5160, 5163, 5165, 5166, 5167, 5168, 5169, 5180, 5181, 5182, 5183, 5400,5448, 5495, 5496, 5505, 5567, 5610, 5617, 5621, 5623, 6736, 6737, 6738, 6739, 6740, 6741, 6742, 6743, 6744, 6746, 6747, 6748, 6749, 6750, 6751, 6753, 6754, 6755, 6756, 6757};
          
  for (int run_id : runs_allow_misc) {
      qa->AllowMiscBit(run_id);
  }

  bool pass_qadb = true;

  //Parcourir fichier et evenements
  while(chain.Next()){

   int run   = c12->runconfig()->getRun();
   int event = c12->runconfig()->getEvent();
   counter2++;

   L.clear();
    
   // get particles by type
   auto electrons=c12->getByID(11);
   auto neutrons=c12->getByID(2112);
   auto Kps=c12->getByID(321);
   auto Kms=c12->getByID(-321);
   auto protons = c12->getByID(2212);

   if(isMC == false){ // si c'est des datas on test le run et l'event 

    pass_qadb = qa->Pass(run,event);

   }

   if(isMC == true){ // si c'est un MC on mets true direct -> pas de qadb

    pass_qadb = true;

   }


   if(pass_qadb == true){

    total_good_events++;
    qa->AccumulateCharge();
    charge_tempo = qa->GetAccumulatedCharge();

    if(counter2 % 1000000 == 0){
		    cout <<" valeur de qa :" << charge_tempo << endl;
	  }

    //SELECTION: ev avec exactement 1 electron >=1 neutron et pi+ pi-


    if(isMC){ // REMPLIR LE TREE GEN A CHAQUE EVENT


          e_vz_gen = c12->getBank(PartMC)->getFloat("vz", 1); // electron tjr en colone 1
          Kp_vz_gen = c12->getBank(PartMC)->getFloat("vz", 2); // kp tjr en colone 2
          Km_vz_gen = c12->getBank(PartMC)->getFloat("vz", 3); // km tjr en colone 3
          n_vz_gen = c12->getBank(PartMC)->getFloat("vz", 0); // neutron tjr en colone 0

          neutron_gen.SetXYZM(c12->getBank(PartMC)->getFloat("px", 0), c12->getBank(PartMC)->getFloat("py", 0), c12->getBank(PartMC)->getFloat("pz", 0), m_target);
          Kp_gen.SetXYZM(c12->getBank(PartMC)->getFloat("px", 2), c12->getBank(PartMC)->getFloat("py", 2), c12->getBank(PartMC)->getFloat("pz", 2), db->GetParticle(321)->Mass());
          Km_gen.SetXYZM(c12->getBank(PartMC)->getFloat("px", 3), c12->getBank(PartMC)->getFloat("py", 3), c12->getBank(PartMC)->getFloat("pz", 3), db->GetParticle(-321)->Mass());
          el_gen.SetXYZM(c12->getBank(PartMC)->getFloat("px", 1), c12->getBank(PartMC)->getFloat("py", 1), c12->getBank(PartMC)->getFloat("pz", 1), db->GetParticle(11)->Mass());

          weight_gen = c12->getBank(EVENTMC)->getFloat("weight", 0);
          real_weight_gen = weight_gen*(Lum/N_gen);

          q_gen = beam - el_gen;
          Q2_gen = -q_gen.M2();
          t_gen = (target - neutron_gen).M2();
          W_gen = (target + q_gen).M();

          miss_gen = beam+target-el_gen-neutron_gen-Kp_gen-Km_gen;

          MinvKpKm_gen = (Kp_gen + Km_gen).M();
          MM_gen = miss_gen.M();

          outT_gen->Fill();

    }


    
    if(neutrons.size()>=1 && electrons.size()==1){

      // ELECTRON :
      SetLorentzVector(el,electrons[0]);
      status_el = electrons[0]->par()->getStatus();

      e_vz = electrons[0]->par()->getVz();
      e_vy = electrons[0]->par()->getVy();
      e_vx = electrons[0]->par()->getVx();


      // NEUTRON :
      for(int i = 0; i < neutrons.size(); i++){

          TLorentzVector neutron_test(0,0,0,m_target);

          SetLorentzVector(neutron_test,neutrons[i]);

          double angle = neutron_test.Vect().Angle(neutron_gen.Vect());

          if(neutron_test.Theta()*180./3.141592 > 4 && neutron_test.P()>0.25){

            L.push_back(angle);

          } else {

            L.push_back(100000000);

          }
        
      }

      auto it = std::min_element(L.begin(), L.end());
      double min_x = *it;
      int idx_min = std::distance(L.begin(), it);

      SetLorentzVector(neutron,neutrons[idx_min]);

      status_n = neutrons[idx_min]->par()->getStatus();

      n_vx = neutrons[idx_min]->par()->getVx();
      n_vy = neutrons[idx_min]->par()->getVy();
      n_vz = neutrons[idx_min]->par()->getVz();

      E_PCAL_n = neutrons[idx_min]->cal(PCAL)->getEnergy(); 
      Lu_PCAL_n = neutrons[idx_min]->cal(PCAL)->getLu(); 
      Lv_PCAL_n = neutrons[idx_min]->cal(PCAL)->getLv(); 
      Lw_PCAL_n = neutrons[idx_min]->cal(PCAL)->getLw(); 
      x_PCAL_n = neutrons[idx_min]->cal(PCAL)->getX(); 
      y_PCAL_n = neutrons[idx_min]->cal(PCAL)->getY(); 

      E_ECIN_n = neutrons[idx_min]->cal(ECIN)->getEnergy(); 
      Lu_ECIN_n = neutrons[idx_min]->cal(ECIN)->getLu(); 
      Lv_ECIN_n = neutrons[idx_min]->cal(ECIN)->getLv(); 
      Lw_ECIN_n = neutrons[idx_min]->cal(ECIN)->getLw(); 
      x_ECIN_n = neutrons[idx_min]->cal(ECIN)->getX(); 
      y_ECIN_n = neutrons[idx_min]->cal(ECIN)->getY(); 

      E_ECOUT_n = neutrons[idx_min]->cal(ECOUT)->getEnergy(); 
      Lu_ECOUT_n = neutrons[idx_min]->cal(ECOUT)->getLu(); 
      Lv_ECOUT_n = neutrons[idx_min]->cal(ECOUT)->getLv(); 
      Lw_ECOUT_n = neutrons[idx_min]->cal(ECOUT)->getLw(); 
      x_ECOUT_n = neutrons[idx_min]->cal(ECOUT)->getX(); 
      y_ECOUT_n = neutrons[idx_min]->cal(ECOUT)->getY(); 

      q = beam - el;

      Q2 = -q.M2();
      t = (target - neutron).M2();
      W = (target + q).M();

      if (isMC){

           weight = c12->getBank(EVENTMC)->getFloat("weight", 0);
           real_weight = weight*(Lum/N_gen);

      }


      outT->Fill();
      counter += 1;

      }

   }
  }



   total_charge = qa->GetAccumulatedCharge();

   std::cout << "==========================" << std::endl;
   std::cout << "Total good events: " << total_good_events << std::endl;
   std::cout << "Total accumulated charge (nC): " << total_charge << std::endl;

   cout << "nombre d'evenements retenus : " << counter <<endl;



   outFile->cd();
   outT->Write("", TObject::kOverwrite);

   if (isMC){

    outT_gen->Write("", TObject::kOverwrite);
    
   }


   outFile->Close();

   gBenchmark->Stop("timer");
   gBenchmark->Print("timer");
  
   auto finish = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed = finish - start;
   std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";

}

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


void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	      rp->par()->getPz(),p4.M());

}

void hipo2root_c12_KpKm(){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();

  bool isMC = false;
  bool isDATA = true;

  /////////////////////////////////////
  //LECTURE FICHIER HIPO

  HipoChain chain;

  // Dossier contenant tes fichiers
  TString folder = "/lustre24/expphy/cache/clas12/rg-b/production/recon/fall2019/torus+1/pass2/v1/dst/train/sidisdvcs"; //

  ///lustre24/expphy/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/nSidis

  int fileCount = 0;
  const int maxFiles = 12;

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
  config_c12->db()->turnOffQADB();

  auto &c12 = chain.C12ref();

  //uint EVENTD = config_c12->addBank("DECAYS::Event");
  //uint PARTD = config_c12->addBank("DECAYS::Particle");

  uint EVENTMC = config_c12->addBank("MC::Event");
  uint PartMC = config_c12->addBank("MC::Particle");



  /////CREER LE FICHIER ROOT ET LES TREEs et variables associées

  cout<<"Variables initi ..."<<endl;

  TFile *outFile = new TFile("fall2019_outbending_neutronKpKm_mc.root", "recreate");   // nom a choisir rgA_outbendingfall2018_data.root or outbending_mc.root for exemple


  TTree *outT = new TTree("tree", "tree");


  TTree *outT_gen = new TTree("tree_gen", "tree_gen");

  double Eb = 10.410;
  double m_target = db->GetParticle(2112)->Mass(); // cas e n -> e' n' K+ K-


  auto db=TDatabasePDG::Instance();
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

  double weight;

  double e_vx_gen, e_vy_gen, e_vz_gen;
  double Kp_vx_gen, Kp_vy_gen, Kp_vz_gen;
  double Km_vx_gen, Km_vy_gen, Km_vz_gen;
  double n_vx_gen, n_vy_gen, n_vz_gen;

  std::vector<double> L;


	outT->Branch("Electron", "TLorentzVector", &el);
	outT->Branch("Neutron", "TLorentzVector", &neutron);
  //outT->Branch("PiPlusD", "TLorentzVector", &pipD);
	//outT->Branch("PiMinusD", "TLorentzVector", &pimD);
  outT->Branch("Kp", "TLorentzVector", &Kp);
  outT->Branch("Km", "TLorentzVector", &Km);
  outT->Branch("Missing", "TLorentzVector", &miss);

  outT->Branch("Q2", &Q2, "Q2/D");
  outT->Branch("t", &t, "t/D");
  outT->Branch("W", &W, "W/D");

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

  if (isMC){

    outT->Branch("weight", &weight, "weight/D");

    outT_gen->Branch("e_vx_gen", &e_vx_gen, "e_vx_gen/D");
    outT_gen->Branch("e_vy_gen", &e_vy_gen, "e_vy_gen/D");
    outT_gen->Branch("e_vz_gen", &e_vz_gen, "e_vz_gen/D");

    outT_gen->Branch("Kp_vx_gen", &Kp_vx_gen, "Kp_vx_gen/D");
    outT_gen->Branch("Kp_vy_gen", &Kp_vy_gen, "Kp_vy_gen/D");
    outT_gen->Branch("Kp_vz_gen", &Kp_vz_gen, "Kp_vz_gen/D");

    outT_gen->Branch("Km_vx_gen", &Km_vx_gen, "Km_vx_gen/D");
    outT_gen->Branch("Km_vy_gen", &Km_vy_gen, "Km_vy_gen/D");
    outT_gen->Branch("Km_vz_gen", &Km_vz_gen, "Km_vz_gen/D");

    outT_gen->Branch("n_vx_gen", &n_vx_gen, "n_vx_gen/D");
    outT_gen->Branch("n_vy_gen", &n_vy_gen, "n_vy_gen/D");
    outT_gen->Branch("n_vz_gen", &n_vz_gen, "n_vz_gen/D");

  }

  
  gBenchmark->Start("timer");
  int counter=0;

  //Parcourir fichier et evenements
  while(chain.Next()){

   L.clear();
    
   // get particles by type
   auto electrons=c12->getByID(11);
   auto neutrons=c12->getByID(2112);
   auto Kps=c12->getByID(321);
   auto Kms=c12->getByID(-321);

   //SELECTION: ev avec exactement 1 electron >=1 neutron et pi+ pi-
    
   if(neutrons.size()>=1 && electrons.size()==1 && Kps.size()==1 && Kms.size()==1){

     // ELECTRON :
     SetLorentzVector(el,electrons[0]);
     status_el = electrons[0]->par()->getStatus();

     e_vz = electrons[0]->par()->getVz();
     e_vy = electrons[0]->par()->getVy();
     e_vx = electrons[0]->par()->getVx();

     if (isMC) {

        e_vx_gen = c12->getBank(PartMC)->getFloat("vx", 1); // e- tjr en colone 1
        e_vy_gen = c12->getBank(PartMC)->getFloat("vy", 1);
        e_vz_gen = c12->getBank(PartMC)->getFloat("vz", 1);

     }

     //CHOIX DES K+ K-

     SetLorentzVector(Kp,Kps[0]);
     SetLorentzVector(Km,Kms[0]);

     status_Kp = Kps[0]->par()->getStatus();
     status_Km = Kms[0]->par()->getStatus();

     Kp_vz = Kps[0]->par()->getVz();
     Kp_vy = Kps[0]->par()->getVy();
     Kp_vx = Kps[0]->par()->getVx();

     Km_vz = Kms[0]->par()->getVz();
     Km_vy = Kms[0]->par()->getVy();
     Km_vx = Kms[0]->par()->getVx();



     if (isMC) {

        Kp_vx_gen = c12->getBank(PartMC)->getFloat("vx", 1); // e- tjr en colone 1
        Kp_vy_gen = c12->getBank(PartMC)->getFloat("vy", 1);
        Kp_vz_gen = c12->getBank(PartMC)->getFloat("vz", 1);

     }


     if (isMC) {

        Km_vx_gen = c12->getBank(PartMC)->getFloat("vx", 1); // e- tjr en colone 1
        Km_vy_gen = c12->getBank(PartMC)->getFloat("vy", 1);
        Km_vz_gen = c12->getBank(PartMC)->getFloat("vz", 1);

     }

     // NEUTRON :
     for(int i = 0; i < neutrons.size(); i++){

        TLorentzVector neutron_test(0,0,0,m_target);

        SetLorentzVector(neutron_test,neutrons[i]);

        miss = beam+target-el-neutron_test-Kp-Km;;
        double delta = miss.M2(); 
        
        L.push_back(delta);
        
     }

     auto it = std::min_element(L.begin(), L.end());
     double min_x = *it;
     int idx_min = std::distance(L.begin(), it);

     SetLorentzVector(neutron,neutrons[idx_min]);

     status_n = neutrons[idx_min]->par()->getStatus();

     n_vx = neutrons[idx_min]->par()->getVx();
     n_vy = neutrons[idx_min]->par()->getVy();
     n_vz = neutrons[idx_min]->par()->getVz();


     if (isMC) {

        p_vx_gen = c12->getBank(PartMC)->getFloat("vx", 0); //neutron tjr en colone 0 ?? (a verifier)
        p_vy_gen = c12->getBank(PartMC)->getFloat("vy", 0);
        p_vz_gen = c12->getBank(PartMC)->getFloat("vz", 0);

     }

     q = beam - el;

     Q2 = -q.M2();
     t = (target - neutron).M2();
     W = (target + q).M();

     miss = beam+target-el-neutron-Kp-Km;

     MinvKpKm = (Kp + Km).M();
     MM = miss.M();

     if (isMC){

          weight = c12->getBank(EVENTMC)->getFloat("weight", 0);

     }

    // cuts larges :

     if((Kp + Km).M2() < 3.0 && (Kp + Km).M2() > 0.0 && miss.M2() < 3.5 && miss.M2() > -3.5) {

       outT->Fill();
       outT_gen->Fill();
       counter += 1;

     }
     

    }
  }

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
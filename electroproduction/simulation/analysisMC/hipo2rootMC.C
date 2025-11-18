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

void hipo2rootMC(){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();


  /////////////////////////////////////
  //LECTURE FICHIER HIPO

  HipoChain chain;

  //TString fname = "nSidis_005169.hipo";

  //chain.Add(fname.Data());
  // Dossier contenant tes fichiers
  TString folder = "/w/hallb-scshelf2102/clas12/mronay/simu/9925"; //

  ///lustre24/expphy/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/nSidis

  int fileCount = 0;
  const int maxFiles = 1000;

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

  uint EVENTMC = config_c12->addBank("MC::Event");
  //uint PARTD = config_c12->addBank("DECAYS::Particle");



  /////CREER LE FICHIER ROOT ET LES TREEs et variables associées

  TFile *outFile = new TFile("inbendingMC.root", "recreate");
  TTree *outT = new TTree("tree", "tree");

  auto db=TDatabasePDG::Instance();
  TLorentzVector beam(0,0, 10.6, 10.6);
  TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector pr(0,0,0,db->GetParticle(2212)->Mass());


  TLorentzVector q(0, 0, 0, 0);

  TLorentzVector pip, pim, miss;

  double Q2, t, W;

  double Minv_pip_pim, MM, Minv_Ks_Kl;

  double status_pim;
  double status_pip;
  double status_pr;
  double status_el;

  double e_vx, e_vy, e_vz;
  double pip_vx, pip_vy, pip_vz;
  double pim_vx, pim_vy, pim_vz;

  double weight;


	outT->Branch("Electron", "TLorentzVector", &el);
	outT->Branch("Proton", "TLorentzVector", &pr);
  outT->Branch("PiPlus", "TLorentzVector", &pip);
	outT->Branch("PiMinus", "TLorentzVector", &pim);
  outT->Branch("Missing", "TLorentzVector", &miss);

  outT->Branch("Q2", &Q2, "Q2/D");
  outT->Branch("t", &t, "t/D");
  outT->Branch("W", &W, "W/D");

  outT->Branch("e_vx", &e_vx, "e_vx/D");
  outT->Branch("e_vy", &e_vy, "e_vy/D");
  outT->Branch("e_vz", &e_vz, "e_vz/D");

  outT->Branch("pip_vx", &pip_vx, "pip_vx/D");
  outT->Branch("pip_vy", &pip_vy, "pip_vy/D");
  outT->Branch("pip_vz", &pip_vz, "pip_vz/D");

  outT->Branch("pim_vx", &pim_vx, "pim_vx/D");
  outT->Branch("pim_vy", &pim_vy, "pim_vy/D");
  outT->Branch("pim_vz", &pim_vz, "pim_vz/D");

  outT->Branch("MinvPiPlusMoins", &Minv_pip_pim, "Minv_pip_pim/D");
  outT->Branch("MissingMass", &MM, "MM/D");
  outT->Branch("MinvKsKl", &Minv_Ks_Kl, "Minv_Ks_Kl/D");

  outT->Branch("Status_el", &status_el, "status_el/D");
  outT->Branch("Status_pim", &status_pim, "status_pim/D");
  outT->Branch("Status_pip", &status_pip, "status_pip/D");
  outT->Branch("Status_pr", &status_pr, "status_pr/D");

  outT->Branch("weight", &weight, "weight/D");

  
   gBenchmark->Start("timer");

   int counter=0;
   int counter2=0;
   int counter3=0;

   //Parcourir fichier et evenements

   cout<<"Variables initi "<<endl;

   while(chain.Next()){

    counter3 += 1;
    
    // get particles by type
    auto electrons=c12->getByID(11);
    auto protons=c12->getByID(2212);
    auto pips=c12->getByID(211);
    auto pims=c12->getByID(-211);

    //SELECTION: ev avec exactement 1 electron 1 proton et pi+ pi-

    if(protons.size()==1 && electrons.size()>=1 && pims.size()>=1 && pips.size()>=1){

      counter2 += 1;

      }
    
    if(protons.size()==1 && electrons.size()==1 && pims.size()==1 && pips.size()==1){

     // PROTON :

     SetLorentzVector(pr,protons[0]);
     status_pr = protons[0]->par()->getStatus();
     
     // ELECTRON :

     SetLorentzVector(el,electrons[0]);
     status_el = electrons[0]->par()->getStatus();

     e_vz = electrons[0]->par()->getVz();
     e_vy = electrons[0]->par()->getVy();
     e_vx = electrons[0]->par()->getVx();

     

     //PI+ PI-

     SetLorentzVector(pip,pips[0]);
     SetLorentzVector(pim,pims[0]);

     status_pim = pims[0]->par()->getStatus();
     status_pip = pips[0]->par()->getStatus();

     pim_vx = pims[0]->par()->getVx();
     pim_vy = pims[0]->par()->getVy();
     pim_vz = pims[0]->par()->getVz();

     pip_vx = pips[0]->par()->getVx();
     pip_vy = pips[0]->par()->getVy();
     pip_vz = pips[0]->par()->getVz();

     q = beam - el;

     Q2 = -q.M2();
     t = (target - pr).M2();
     W = (target + q).M();

     miss = beam+target-el-pr-pip-pim;

     Minv_pip_pim = (pim + pip).M();
     MM = miss.M();
     Minv_Ks_Kl = (miss + pip + pim).M();


     //Cut sur masse inv pi+ pi-
     //Sur missing mass 
     //Sur masse inv Ks + Kl

     weight = c12->getBank(EVENTMC)->getFloat("weight", 0);

     if((pim + pip).M() < 3.0 && (pim + pip).M() > 0.0 && miss.M() < 3.0 && miss.M() > 0.0 && (miss + pip + pim).M() < 3.0 && (miss + pip + pim).M() > 0.0 ) {

       outT->Fill();

       counter += 1;

     }
     

    }
  }

  cout << "nombre d'evenements total reconstruit est : " << counter3 <<endl;
  cout << "nombre d'evenements retenus avec exactement 1 pi+-: " << counter <<endl;
  cout << "nombre d'evenements retenus >=1 pi+-: " << counter2 <<endl;
  

   outFile->cd();
   outT->Write("", TObject::kOverwrite);
   outFile->Close();

   gBenchmark->Stop("timer");
   gBenchmark->Print("timer");
  
   auto finish = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed = finish - start;
   std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";

}
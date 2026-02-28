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

void hipo2root_c12(){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();


  /////////////////////////////////////
  //LECTURE FICHIER HIPO

  HipoChain chain;

  //TString fname = "nSidis_005169.hipo";

  //chain.Add(fname.Data());
  // Dossier contenant tes fichiers
  TString folder = "/w/hallb-scshelf2102/clas12/mronay/Displaced_Vertex/output_displace_vertex"; //

  ///lustre24/expphy/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/nSidis

  int fileCount = 0;
  const int maxFiles = 50;

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

  uint EVENTD = config_c12->addBank("DECAYS::Event");
  uint PARTD = config_c12->addBank("DECAYS::Particle");



  /////CREER LE FICHIER ROOT ET LES TREEs et variables associées

  TFile *outFile = new TFile("rgA_outbending.root", "recreate");
  TTree *outT = new TTree("tree", "tree");

  auto db=TDatabasePDG::Instance();
  TLorentzVector beam(0,0, 10.6, 10.6);
  TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector pr(0,0,0,db->GetParticle(2212)->Mass());


  TLorentzVector q(0, 0, 0, 0);
  TLorentzVector Ks(0,0,0,0);

  TLorentzVector pipD, pimD, missD;

  double Q2, t, W;

  double pid;
  double index_pr, index_pip, index_pim;
  double Minv_pipD_pimD, MMD, Minv_Ks_KlD;

  double status_pim;
  double status_pip;
  double status_pr;
  double status_el;

  double number_pr;
  double number_pip;
  double number_pim;
  double number_el;
  
  double number_event;

  double e_vx, e_vy, e_vz;
  double pip_vx, pip_vy, pip_vz;
  double pim_vx, pim_vy, pim_vz;
  double Ks_vx, Ks_vy, Ks_vz;

  std::vector<double> L;
  std::vector<double> L2;
  

	outT->Branch("Electron", "TLorentzVector", &el);
	outT->Branch("Proton", "TLorentzVector", &pr);
  outT->Branch("PiPlusD", "TLorentzVector", &pipD);
	outT->Branch("PiMinusD", "TLorentzVector", &pimD);
  outT->Branch("MissingD", "TLorentzVector", &missD);

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

  outT->Branch("Ks_vx", &Ks_vx, "Ks_vx/D");
  outT->Branch("Ks_vy", &Ks_vy, "Ks_vy/D");
  outT->Branch("Ks_vz", &Ks_vz, "Ks_vz/D");

  outT->Branch("MinvPiPlusMoinsD", &Minv_pipD_pimD, "Minv_pipD_pimD/D");
  outT->Branch("MissingMassD", &MMD, "MMD/D");
  outT->Branch("MinvKsKlD", &Minv_Ks_KlD, "Minv_Ks_KlD/D");

  outT->Branch("Status_pim", &status_pim, "status_pim/D");
  outT->Branch("Status_pip", &status_pip, "status_pip/D");
  outT->Branch("Status_pr", &status_pr, "status_pr/D");

  outT->Branch("number_pim", &number_pim, "number_pim/D");
  outT->Branch("number_pip", &number_pip, "number_pip/D");
  outT->Branch("number_pr", &number_pr, "number_pr/D");
  outT->Branch("number_el", &number_el, "number_el/D");
  outT->Branch("number_event", &number_event, "number_event/D");
  
   gBenchmark->Start("timer");
   int counter=0;

   //Parcourir fichier et evenements

   cout<<"Variables initi "<<endl;

   while(chain.Next()){
    L.clear();
    L2.clear();
    
    // get particles by type
    auto electrons=c12->getByID(11);
    auto protons=c12->getByID(2212);
    auto pips=c12->getByID(211);
    auto pims=c12->getByID(-211);

    //SELECTION: ev avec exactement 1 electron 1 proton et pi+ pi-
    
    if(protons.size()==1 && electrons.size()>=1 && pims.size()>=1 && pips.size()>=1){

     number_el = electrons.size();
     number_pr = protons.size();
     number_pip = pips.size();
     number_pim = pims.size();

     number_event = number_pr*1000 + number_el*100 + number_pip*10 + number_pim; // normalement bcp de 1111 et 1122


     // CHOIX DU PROTON :
     SetLorentzVector(pr,protons[0]);
     status_pr = protons[0]->par()->getStatus();
     
     // CHOIX DE L'ELECTRON :
     for(int i = 0; i < electrons.size(); i++){

        e_vz = electrons[i]->par()->getVz();
        e_vy = electrons[i]->par()->getVy();
        e_vx = electrons[i]->par()->getVx();

        double delta = sqrt((e_vz + 2.5)*(e_vz + 2.5) + (e_vy - 0)*(e_vy - 0) + (e_vx - 0)*(e_vx - 0)); // abs(evz - (-2.5))

        status_el = electrons[i]->par()->getStatus();

        if (abs(status_el) <= 2999 && abs(status_el) >= 2000){ // lectron dans le FD

          L.push_back(delta);

        } else {

          L.push_back(1000); // on add une distance qui sera jamais choisi 

        }

        if (counter % 100000 ==0){

              // Affiche chaque valeur ajoutée
              std::cout << "Electron " << i 
              << ", delta = " << delta << std::endl;

        }
        
     }

     auto it = std::min_element(L.begin(), L.end());
     double min_x = *it;
     int idx_min = std::distance(L.begin(), it);

     SetLorentzVector(el,electrons[idx_min]);

     if (counter % 100000 == 0){

          std::cout << "\n--- Résumé ---" << std::endl;
          std::cout << "Liste des deltas : ";
          for (size_t j = 0; j < L.size(); j++) {
             std::cout << L[j] << " ";
          }
          std::cout << std::endl;

          std::cout << "Valeur delta minimale = " << min_x << std::endl;
          std::cout << "Indice du delta  minimum = " << idx_min << std::endl;

          std::cout << "----------------\n" << std::endl;


     }

     //CHOIX DES PI+ PI-

      auto bank = c12->getBank(PARTD);
      auto schema = bank->getSchema();
      int nb_Rows = bank->getRows();
      int nbKs = nb_Rows/3;

      for(int i = 0; i < nbKs; i++){

        for (int j = 0; j < 3; j++) {

          pid = c12->getBank(PARTD)->getInt("pid", i*3 + j);

          if(counter % 100000 == 0){

             cout << "Le pid avant le if vaut " << pid << " pour le lidentifiant boucle : " << i*3 + j <<endl;

          }

          if (pid == 310){

              e_vx = c12->getBank(PARTD)->getFloat("ovx", i*3 + j);
              e_vy = c12->getBank(PARTD)->getFloat("ovy", i*3 + j);
              e_vz = c12->getBank(PARTD)->getFloat("ovz", i*3 + j);

              Ks_vx = c12->getBank(PARTD)->getFloat("vx", i*3 + j);
              Ks_vy = c12->getBank(PARTD)->getFloat("vy", i*3 + j);
              Ks_vz = c12->getBank(PARTD)->getFloat("vx", i*3 + j);

              if (counter % 100000 == 0){

                  cout << "dans la boucle Ks "  <<endl;
                  cout << "le pid vaut  : " << pid << " et l'identifiant vaut : " << i*3 + j << endl;
              }

          }



          if (pid == 211){

               pipD.SetXYZM(c12->getBank(PARTD)->getFloat("px", i*3 + j), c12->getBank(PARTD)->getFloat("py", i*3 + j), c12->getBank(PARTD)->getFloat("pz", i*3 + j), db->GetParticle(211)->Mass());

               pip_vx = c12->getBank(PARTD)->getFloat("vx", i*3 + j);
               pip_vy = c12->getBank(PARTD)->getFloat("vy", i*3 + j);
               pip_vz = c12->getBank(PARTD)->getFloat("vx", i*3 + j);

              if(counter % 100000 == 0){

                  cout << "dans la boucle pi+ "  <<endl;
                  cout << "le pid vaut  : " << pid << " et l'identifiant vaut : " << i*3 + j << endl;


              }

          }

          if (pid == -211){

                pimD.SetXYZM(c12->getBank(PARTD)->getFloat("px", i*3 + j), c12->getBank(PARTD)->getFloat("py", i*3 + j), c12->getBank(PARTD)->getFloat("pz", i*3 + j), db->GetParticle(-211)->Mass());

                pim_vx = c12->getBank(PARTD)->getFloat("vx", i*3 + j);
                pim_vy = c12->getBank(PARTD)->getFloat("vy", i*3 + j);
                pim_vz = c12->getBank(PARTD)->getFloat("vx", i*3 + j);

                if(counter % 100000 == 0){

                    cout << "dans la boucle pi- " <<endl;
                    cout << "le pid vaut  : " << pid << " et l'identifiant vaut : " << i*3 + j <<endl;

                }
           }
          

        }

        double Mks = (pimD + pipD).M();
        double deltaMks = abs(Mks - 0.497); // verif proche masse Ks
        L2.push_back(deltaMks);

      }

        auto it2 = std::min_element(L2.begin(), L2.end());
        double min_x2 = *it2;
        int idx_min2 = std::distance(L2.begin(), it2);

        if (counter % 100000 == 0){

          std::cout << "\n--- Résumé ---" << std::endl;
          std::cout << "Liste des deltasMKs : ";
          for (size_t j = 0; j < L2.size(); j++) {
             std::cout << L2[j] << " ";
          }
          std::cout << std::endl;

          std::cout << "Valeur deltaMKs minimale = " << min_x2 << std::endl;
          std::cout << "Indice du deltaMKs  minimum = " << idx_min2 << std::endl;

          std::cout << "----------------\n" << std::endl;


        }

        for (int j = 0; j < 3; j++) {

          pid = c12->getBank(PARTD)->getInt("pid", idx_min2*3 + j);

          if(counter % 100000 == 0){

             cout << "VRAI ATTRIBUTION : le pid avant le if vaut " << pid << " pour le lidentifiant boucle : " << idx_min2*3 + j <<endl;

          }

          if (pid == 310){

              e_vx = c12->getBank(PARTD)->getFloat("ovx", idx_min2*3 + j);
              e_vy = c12->getBank(PARTD)->getFloat("ovy", idx_min2*3 + j);
              e_vz = c12->getBank(PARTD)->getFloat("ovz", idx_min2*3 + j);

              Ks_vx = c12->getBank(PARTD)->getFloat("vx", idx_min2*3 + j);
              Ks_vy = c12->getBank(PARTD)->getFloat("vy", idx_min2*3 + j);
              Ks_vz = c12->getBank(PARTD)->getFloat("vx", idx_min2*3 + j);


              if (counter % 100000 == 0){

                  cout << "dans la boucle Ks "  <<endl;
                  cout << "le pid vaut  : " << pid << " et l'identifiant vaut : " << idx_min2*3 + j << endl;

              }

          }



          if (pid == 211){

               pipD.SetXYZM(c12->getBank(PARTD)->getFloat("px", idx_min2*3 + j), c12->getBank(PARTD)->getFloat("py", idx_min2*3 + j), c12->getBank(PARTD)->getFloat("pz", idx_min2*3 + j), db->GetParticle(211)->Mass());

               pip_vx = c12->getBank(PARTD)->getFloat("vx", idx_min2*3 + j);
               pip_vy = c12->getBank(PARTD)->getFloat("vy", idx_min2*3 + j);
               pip_vz = c12->getBank(PARTD)->getFloat("vx", idx_min2*3 + j);

               index_pip = c12->getBank(PARTD)->getInt("idx", idx_min2*3 + j);

               for(int k=0;k<pips.size();k++){

                  if(pips[k]->getIndex() == index_pip){

                    status_pip = pips[k]->par()->getStatus();

                  }
               }
    

               if(counter % 100000 == 0){

                  cout << "dans la boucle pi+ "  <<endl;
                  cout << "le pid vaut  : " << pid << " et l'identifiant vaut : " << idx_min2*3 + j << endl;

               }

          }

          if (pid == -211){

                pimD.SetXYZM(c12->getBank(PARTD)->getFloat("px", idx_min2*3 + j), c12->getBank(PARTD)->getFloat("py", idx_min2*3 + j), c12->getBank(PARTD)->getFloat("pz", idx_min2*3 + j), db->GetParticle(-211)->Mass());

                pim_vx = c12->getBank(PARTD)->getFloat("vx", idx_min2*3 + j);
                pim_vy = c12->getBank(PARTD)->getFloat("vy", idx_min2*3 + j);
                pim_vz = c12->getBank(PARTD)->getFloat("vx", idx_min2*3 + j);

                index_pim = c12->getBank(PARTD)->getInt("idx", idx_min2*3 + j);

                for(int k=0;k<pims.size();k++){

                  if(pims[k]->getIndex() == index_pim){

                    status_pim = pims[k]->par()->getStatus();

                  }
                }

                if(counter % 100000 == 0){

                    cout << "dans la boucle pi- " <<endl;
                    cout << "le pid vaut  : " << pid << " et l'identifiant vaut : " << idx_min2*3 + j <<endl;

                }
           }
          

        }


     q = beam - el;

     Q2 = -q.M2();
     t = (target - pr).M2();
     W = (target + q).M();

     // e_vz = electrons[0]->par()->getVz();

     // a revoir avec le getindex et l'idx dans DECAYS

    //status_pim = pims[0]->par()->getStatus();
    // status_pip = pips[0]->par()->getStatus();
    // status_pr = protons[0]->par()->getStatus();


     missD = beam+target-el-pr-pipD-pimD;

     Minv_pipD_pimD = (pimD + pipD).M();
     MMD = missD.M();
     Minv_Ks_KlD = (missD + pipD + pimD).M();


     //Cut sur masse inv pi+ pi-
     //Sur missing mass 
     //Sur masse inv Ks + Kl

     if((pimD + pipD).M() < 3.0 && (pimD + pipD).M() > 0.0 && missD.M() < 3.0 && missD.M() > 0.0 && (missD + pipD + pimD).M() < 3.0 && (missD + pipD + pimD).M() > 0.0 ) {

       outT->Fill();

       counter += 1;

     }
     

    }
  }

  cout << "nombre d'evenements retenus : " <<counter<<endl;

   outFile->cd();
   outT->Write("", TObject::kOverwrite);
   outFile->Close();

   gBenchmark->Stop("timer");
   gBenchmark->Print("timer");
  
   auto finish = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed = finish - start;
   std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";

}
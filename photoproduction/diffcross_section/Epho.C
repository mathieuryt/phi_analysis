#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include <fstream>
#include <iostream>
#include <string>

// ce programme permets d'avoir les bornes des bins et les centre des bins de E photon 

void Epho()
{
   gEnv->SetValue("Browser.Name", "TRootBrowser");

   TChain *t1 = new TChain("tree");  // "tree" doit être le nom EXACT du TTree dans les fichiers

   // Ajoute tous tes fichiers ROOT
   t1->Add("/local/home/mr282803/Documents/irfu/diff_cs_new/data/Data_pass2_fall2018_outbending.root");
   t1->Add("/local/home/mr282803/Documents/irfu/diff_cs_new/data/Data_pass2_spring2019_inbending.root");
   t1->Add("/local/home/mr282803/Documents/irfu/diff_cs_new/data/Data_pass2_fall2018_inbending.root");

   TLorentzVector *electron = nullptr;
   TLorentzVector *positron = nullptr;

   Float_t Q2, MMassBeam, counter, Ephovar;

   double plage_a = 2;
   double plage_b = 12;
   double nb_bins = 15000;
   double nb_ev_inter = 7500;
   float taille_bins_ini = (plage_b -plage_a)/nb_bins;

   t1->SetBranchAddress("Electron", &electron);
   t1->SetBranchAddress("Positron", &positron);
   t1->SetBranchAddress("Q2", &Q2);
   t1->SetBranchAddress("MMassBeam", &MMassBeam);
   t1->SetBranchAddress("Epho", &Ephovar);

   counter = 0;
   // Créer un nouvel histogramme
   TH1F *Ephoton = new TH1F("Ephoton", "Distribution de l'energie du photon", nb_bins, plage_a, plage_b);

   Long64_t nentries = t1->GetEntries();

   for (Long64_t i = 0; i < nentries; i++) {
       t1->GetEntry(i);

       if (!electron || !positron) continue;
       double M = (*electron + *positron).M();

       if (M<1.09 && M>0.95 && abs(MMassBeam) < 0.3 && Q2 < 0.3){
           Ephoton->Fill(Ephovar);
       }
   }

   // Affichage de l'histogramme
   Ephoton->SetLineColor(kBlue);
   Ephoton->GetXaxis()->SetTitle("Masse invariante [GeV]");
   Ephoton->GetYaxis()->SetTitle("Nombre d'evenements");
   TCanvas *c1 = new TCanvas("c1", "Masse Invariante", 800, 600);
   Ephoton->Draw();

   std::vector<float> liste;
   liste.push_back(plage_a);

   double counter_pho = 0;
   double moy = 0;
   double new_nb_bin = 0;

   std::vector<float> liste_E_moy;


   int nBins = Ephoton->GetNbinsX();
   for (int i = 1; i <= nBins; ++i) {  // Les bins vont de 1 à nBins
        double binCenter = Ephoton->GetBinCenter(i);
        moy += Ephoton->GetBinCenter(i)*Ephoton->GetBinContent(i);
        double binContent = Ephoton->GetBinContent(i);
        std::cout << "Bin " << i << " (Centre = " << binCenter << ") : " 
                  << binContent << " événements" << std::endl;
        counter_pho += binContent;
        if (counter_pho >= nb_ev_inter) {

            double newbinCenter = moy/nb_ev_inter;
            liste_E_moy.push_back(newbinCenter);

            counter_pho = 0;
            moy = 0;
            liste.push_back(plage_a + taille_bins_ini*i);
            new_nb_bin += 1;
       }
   }



   std::cout << "Liste des bornes bins : ";
   for (float val : liste) std::cout << val << " ";
   std::cout << std::endl;



   std::cout << "le nombre nouveau nombre de bins est = " << new_nb_bin << std::endl;

   TH1F *newEphoton = new TH1F("newEphoton", "Histogramme avec bins variables",
    new_nb_bin, liste.data());

 
    std::cout << "le nombre d'entries avant coupures est = " << nentries << std::endl;
    for (Long64_t i = 0; i < nentries; i++) {
        t1->GetEntry(i);
 
        if (!electron || !positron) continue;
        double M = (*electron + *positron).M();
 
        if (M<1.09 && M>0.95 && abs(MMassBeam) < 0.3 && Q2 < 0.3) {
            counter++;
            newEphoton->Fill(Ephovar);
        }
    }
 
    std::cout << "Nombre de photon interressant total  = " << counter << std::endl;
    counter = 0;
    // Affichage de l'histogramme
    newEphoton->SetLineColor(kBlue);
    newEphoton->GetXaxis()->SetTitle("Masse invariante [GeV]");
    newEphoton->GetYaxis()->SetTitle("Nombre d'evenements");
    newEphoton->GetYaxis()->SetRangeUser(0, 5000);
    TCanvas *c2 = new TCanvas("c2", "new Masse Invariante", 800, 600);
    newEphoton->Draw();


   cout << "Histogramme terminé" << endl;




   std::cout << "Liste des centre de bins en énergie : ";
   for (float val : liste_E_moy) std::cout << val << " ";
   std::cout << std::endl;

   
   std::ofstream fichier("fichier_E_moy.txt");

   if (fichier) {
       for (float valeur : liste_E_moy) {
           fichier << valeur << std::endl;
       }
       fichier.close();
       std::cout << "Fichier des E_moy (centre bins) écrit avec succès !" << std::endl;
   } else {
       std::cerr << "Erreur lors de l'ouverture du fichier !" << std::endl;
   }


   std::ofstream fichier_2("fichier_bornes_bins.txt");

   if (fichier_2) {
       for (float valeur : liste) {
           fichier_2 << valeur << std::endl;
       }
       fichier_2.close();
       std::cout << "Fichier des bornes bins écrit avec succès !" << std::endl;
   } else {
       std::cerr << "Erreur lors de l'ouverture du fichier !" << std::endl;
   }

}


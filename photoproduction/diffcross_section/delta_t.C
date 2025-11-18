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
#include <vector>

// ce programme permets d'avoir les bornes des bins et les centre des bins de E photon 

double Lambda(double x, double y, double z) {
    return (x - y - z)*(x - y - z) - 4 * y*z;
}

//From Byukling Kayanti Formula (5.14) Page 86

double T_min(double E) // arguments are squares of masses of particles in the reaction a+b->1+2, and s is the square of the total c.m. energy i.e. (a+b)^2
{
    const double Mp = 0.938;
    const double M_phi = 1.019;
    double ma_2 = 0;
    double mb_2 = Mp*Mp;

    double m1_2 = 1.019*1.019;
    double m2_2 = Mp*Mp;
    double s = Mp*Mp + 2*Mp*E;
    return abs(ma_2 + m1_2 - (1 / (2 * s))*((s + ma_2 - mb_2)*(s + m1_2 - m2_2) - sqrt(Lambda(s, ma_2, mb_2) * Lambda(s, m1_2, m2_2))));
}


void delta_t()
{  
   //ini data set
   gEnv->SetValue("Browser.Name", "TRootBrowser");

   TChain *t1 = new TChain("tree");  // "tree" doit être le nom EXACT du TTree dans les fichiers

   // Ajoute tous tes fichiers ROOT
   t1->Add("/local/home/mr282803/Documents/irfu/diff_cs_new/data/Data_pass2_fall2018_outbending.root");
   t1->Add("/local/home/mr282803/Documents/irfu/diff_cs_new/data/Data_pass2_spring2019_inbending.root");
   t1->Add("/local/home/mr282803/Documents/irfu/diff_cs_new/data/Data_pass2_fall2018_inbending.root");

   TLorentzVector *electron = nullptr;
   TLorentzVector *positron = nullptr;
   TLorentzVector *proton = nullptr;

   Float_t Q2, MMassBeam, counter, Ephovar;

   t1->SetBranchAddress("Electron", &electron);
   t1->SetBranchAddress("Positron", &positron);
   t1->SetBranchAddress("Proton", &proton);
   t1->SetBranchAddress("Q2", &Q2);
   t1->SetBranchAddress("MMassBeam", &MMassBeam);
   t1->SetBranchAddress("Epho", &Ephovar);

   //physique
   double Mp = 0.938;

   //plage global en t
   double plage_a = 0;
   double plage_b = 2;

   //borne E_gamma
   double elim_a = 0;
   double elim_b = 0;

   //bins sur t initialement
   double nb_bins = 15000;
   float taille_bins_ini = (plage_b -plage_a)/nb_bins;

   //nb interressent par bins variable en t 
   double nb_ev_inter = 1050;

   //on recup la liste des bins en E_gamma

   std::string nomFichier = "/local/home/mr282803/Documents/irfu/diff_cs_new/diff_bins/fichier_bornes_bins.txt";
   std::ifstream fichier(nomFichier);
   std::vector<double> bornes;

   double valeur;
   while (fichier >> valeur) {
        bornes.push_back(valeur);
    }

   fichier.close();

   counter = 0;

   // PARTIE I SEQUANCE AVEC DES TOUT PETIT BIN EN |t| 

   // initialisation des histo initiaux en t

   std::vector<TH1F*> histograms;

   std::cout << "taille de la liste des bornes en E gamma " << bornes.size() << std::endl;

    for (int i = 0; i < bornes.size()-1; i++) {
        std::cout << "Initialisation de l'histogramme numéro : " << i << std::endl;
        histograms.push_back(new TH1F(Form("Histogramme numero %d", i), "", nb_bins, plage_a, plage_b));
    }

    // remplissage des histo ini

    for (int i = 0; i < bornes.size()-1; i++) {
        if (i + 1 >= bornes.size()) {
            std::cerr << "Erreur : accès hors limite à liste[" << i+1 << "]" << std::endl;
            continue;
        }

        elim_a = bornes[i];
        elim_b = bornes[i+1];

        std::cout << "Remplissage histo ini pour Epho compris entre " << elim_a << " et " << elim_b << std::endl;

        Long64_t nentries = t1->GetEntries();

        for (Long64_t j = 0; j < nentries; j++) {
            t1->GetEntry(j);
     
            if (!electron || !positron) continue;

            //masse invariante
            double M = (*electron + *positron).M();


            // calcul du quadrivecteur k = p -p'

            Double_t E = Mp - proton->E();
            Double_t px = 0 - proton->Px();
            Double_t py = 0 - proton->Py();
            Double_t pz = 0 - proton->Pz();

            // calcul de |t| = abs(k²)
            Double_t t = abs(E*E - px*px - py*py - pz*pz);
     
            if (M<1.09 && M>0.95 && abs(MMassBeam) < 0.3 && Q2 < 0.3 && Ephovar > elim_a && Ephovar < elim_b){
                histograms[i]->Fill(t);
            }
        }

        // Affichage de l'histogramme
        histograms[i]->SetLineColor(kBlue);
        histograms[i]->GetXaxis()->SetTitle("|t|[GeV]");
        histograms[i]->GetYaxis()->SetTitle("Nombre d'evenements");
        // Créer un nom dynamique pour le TCanvas avec l'index k
        std::string canvasName = "c" + std::to_string(i);  
        TCanvas *c = new TCanvas(canvasName.c_str(), "Histogram |t| initiaux", 800, 600);
        histograms[i]->Draw();

    }

   // PARTIE 2 ON AJUSTE LA TAILLE DES BINS 

   // liste de liste pour chaque decoupage en t par E_gamma

   std::vector<std::vector<double>> liste_de_listes(bornes.size()-1);
   std::vector<std::vector<double>> liste_E_moy(bornes.size()-1);


   // liste contenant les nb de bins en |t| par bins de E_gamma
   std::vector<float> new_nb_bin(bornes.size()-1);

   //on ajouter la premier borne en |t| pour chaque liste

   for (int i = 0; i < bornes.size()-1; i++) {

    //liste_de_listes[i].push_back(T_min(bornes[i+1]));
    liste_de_listes[i].push_back(0.1);

   }


   for (int i = 0; i < bornes.size()-1; i++) {
    
    double counter_pho = 0; // on re ini bien a chaque fois qu'on change de borne en E gamma
    double moy = 0;


    for (int j = 1; j <= nb_bins; ++j) {  // Les bins vont de 1 à nb_bins

        double binCenter = histograms[i]->GetBinCenter(j);
        moy += histograms[i]->GetBinCenter(j)*histograms[i]->GetBinContent(j);
        double binContent = histograms[i]->GetBinContent(j);

        std::cout << "Bin " << j << " (Centre = " << binCenter << ") : " 
                  << binContent << " événements" << std::endl;

        counter_pho += binContent;
        if (counter_pho >= nb_ev_inter) {

            double newbinCenter = moy/nb_ev_inter;

            liste_E_moy[i].push_back(newbinCenter);
            liste_de_listes[i].push_back(plage_a + taille_bins_ini*j);
            new_nb_bin[i] = new_nb_bin[i] + 1;

            counter_pho = 0;
            moy = 0;
        }
    }

   }

   // PARTIE 3 RESULTATS

   for (int i = 0; i < bornes.size()-1; i++) {

    elim_a = bornes[i];
    elim_b = bornes[i+1];

    std::cout << "Pour E_gamma compris entre " << elim_a << "et " << elim_b <<std::endl;
    std::cout << "La liste des bins en |t| est : " << std::endl;

    for (float val : liste_de_listes[i]) std::cout << val << " ";
    std::cout << std::endl;

    std::cout << "Et la liste des centres de bins moyen (pondere) est : "<< std::endl;

    for (float val : liste_E_moy[i]) std::cout << val << " ";
    std::cout << std::endl;

   }


   std::cout << "Pour chaque bins en E_gamma le nb de bin en |t| est "<< std::endl;

   for (float val : new_nb_bin) std::cout << val << " ";
   std::cout << std::endl;



   // PARTIE 4 VERIFICATION

   std::vector<TH1F*> histograms_verif;

   for (int i = 0; i < bornes.size()-1; i++) {
       std::cout << "Initialisation de l'histogramme de verif numéro : " << i << std::endl;
       histograms_verif.push_back(new TH1F(Form("Histogramme verif numero %d", i), "", new_nb_bin[i], liste_de_listes[i].data()));
   }

     // remplissage des histo verif

    for (int i = 0; i < bornes.size()-1; i++) {
        if (i + 1 >= bornes.size()) {
            std::cerr << "Erreur : accès hors limite à liste[" << i+1 << "]" << std::endl;
            continue;
        }

        elim_a = bornes[i];
        elim_b = bornes[i+1];

        std::cout << "Remplissage histo verif pour Epho compris entre " << elim_a << " et " << elim_b << std::endl;

        Long64_t nentries = t1->GetEntries();

        for (Long64_t j = 0; j < nentries; j++) {
            t1->GetEntry(j);
     
            if (!electron || !positron) continue;

            //masse invariante
            double M = (*electron + *positron).M();


            // calcul du quadrivecteur k = p -p'

            Double_t E = Mp - proton->E();
            Double_t px = 0 - proton->Px();
            Double_t py = 0 - proton->Py();
            Double_t pz = 0 - proton->Pz();

            // calcul de |t| = abs(k²)
            Double_t t = abs(E*E - px*px - py*py - pz*pz);
     
            if (M<1.09 && M>0.95 && abs(MMassBeam) < 0.3 && Q2 < 0.3 && Ephovar > elim_a && Ephovar < elim_b){
                histograms_verif[i]->Fill(t);
            }
        }

        // Affichage de l'histogramme
        histograms_verif[i]->SetLineColor(kBlue);
        histograms_verif[i]->GetXaxis()->SetTitle("|t|[GeV]");
        histograms_verif[i]->GetYaxis()->SetTitle("Nombre d'evenements");
        histograms_verif[i]->GetYaxis()->SetRangeUser(200, 300);
        // Créer un nom dynamique pour le TCanvas avec l'index i
        std::string canvasName = "c" + std::to_string(i);  
        TCanvas *c = new TCanvas(canvasName.c_str(), "Histogram |t| verif", 800, 600);
        histograms_verif[i]->Draw();

    }


    // PARTIE 5 ECRITURE SAUV

    for (size_t i = 0; i < liste_de_listes.size(); ++i) {
        // Créer un nom de fichier unique pour chaque sous-liste
        std::string nom_fichier = "fichier_bornes_bins_" + std::to_string(i) + ".txt";
        std::ofstream fichier(nom_fichier);

        if (fichier) {
            for (float valeur : liste_de_listes[i]) {
                fichier << valeur << std::endl;
            }
            fichier.close();
            std::cout << "Fichier " << nom_fichier << " écrit avec succès !" << std::endl;
        } else {
            std::cerr << "Erreur lors de l'ouverture du fichier : " << nom_fichier << std::endl;
        }
    }

    for (size_t i = 0; i < liste_E_moy.size(); ++i) {
        // Créer un nom de fichier unique pour chaque sous-liste
        std::string nom_fichier = "fichier_E_moy_" + std::to_string(i) + ".txt";
        std::ofstream fichier(nom_fichier);

        if (fichier) {
            for (float valeur : liste_E_moy[i]) {
                fichier << valeur << std::endl;
            }
            fichier.close();
            std::cout << "Fichier " << nom_fichier << " écrit avec succès !" << std::endl;
        } else {
            std::cerr << "Erreur lors de l'ouverture du fichier : " << nom_fichier << std::endl;
        }
    }





}


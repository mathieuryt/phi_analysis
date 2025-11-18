#include "TH1F.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include <iostream>
#include <vector>


void M2Q()  
{
    gEnv->SetValue("Browser.Name", "TRootBrowser");

    float Q2_max = 0.3;
    float MMassBeam_max = 0.3;
    float bins = 100.0;


    TChain *t1 = new TChain("tree");  // "tree" doit être le nom EXACT du TTree dans les fichiers

    // Ajoute tous tes fichiers ROOT
    t1->Add("Data_pass2_fall2018_outbending.root");
    t1->Add("Data_pass2_spring2019_inbending.root");
    t1->Add("Data_pass2_fall2018_inbending.root");

    TLorentzVector *electron = nullptr;
    TLorentzVector *positron = nullptr;
    TLorentzVector *proton = nullptr;

    Float_t Q2, MMassBeam, Epho, elim_a, elim_b;
    double tlim_a = 0;
    double tlim_b = 0;
    double Mp = 0.938;
    Float_t counter = 0; 
    Float_t counterpho = 0;

    t1->SetBranchAddress("Electron", &electron);
    t1->SetBranchAddress("Positron", &positron);
    t1->SetBranchAddress("Proton", &proton);
    t1->SetBranchAddress("Q2", &Q2);
    t1->SetBranchAddress("MMassBeam", &MMassBeam);
    t1->SetBranchAddress("Epho", &Epho);

    std::vector<TH1F*> histograms;

    std::string nomFichier = "/local/home/mr282803/Documents/irfu/diff_cs_new/diff_bins/fichier_bornes_bins.txt";
    std::ifstream fichier(nomFichier);
    std::vector<double> bornes;
 
    double valeur;
    while (fichier >> valeur) {
         bornes.push_back(valeur);
     }
 
    fichier.close();

    // Affichage du contenu du vecteur
    std::cout << "Liste des bins en E_gamma : ";
    for (float val : bornes) std::cout << val << " ";
    std::cout << std::endl;

    int new_nb_bin = bornes.size()-1;

    std::vector<std::vector<double>> liste_de_listes(bornes.size()-1);

    size_t k = bornes.size()-1; // Par exemple, on a 5 fichiers de 0 à 4

    for (size_t i = 0; i < k; ++i) {
        std::string nom_fichier = "/local/home/mr282803/Documents/irfu/diff_cs_new/diff_bins/fichier_bornes_bins_" + std::to_string(i) + ".txt";
        std::ifstream fichier(nom_fichier);

        if (fichier) {
            double valeur;
            while (fichier >> valeur) {
                liste_de_listes[i].push_back(valeur);
            }
            fichier.close();
            std::cout << "Fichier " << nom_fichier << " lu avec succès." << std::endl;
        } else {
            std::cerr << "Erreur lors de l'ouverture du fichier : " << nom_fichier << std::endl;
        }
    }

    // (Optionnel) Vérification : afficher le contenu
    for (size_t i = 0; i < liste_de_listes.size(); ++i) {
        std::cout << "liste des bornes bins en |t| pour le bin en E_gamma numero " << std::endl;
        for (double val : liste_de_listes[i]) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }


    // Initialisation des histogrammes
    for (int i = 0; i < new_nb_bin; i++) {
        for (int j =0; j<liste_de_listes[i].size()-1; j++) {

            std::cout << "Initialisation de l'histogramme pour le bin en E_gamma " << i << " et le bin en |t| " << j << std::endl;
            histograms.push_back(new TH1F(Form("Histogramme_%d_%d", i, j), "", bins, 0.8, 1.2));

        }
    }

    std::cout << "Taille du vecteur d'histogrammes : " << histograms.size() << std::endl;

    double compteur_histo = 0;

    for (int i = 0; i < new_nb_bin; i++) {
  
        elim_a = bornes[i];
        elim_b = bornes[i+1];

        for (int j = 0; j < liste_de_listes[i].size()-1; j++) {

            tlim_a = liste_de_listes[i][j];
            tlim_b = liste_de_listes[i][j+1];

            std::cout << "Remplissage pour E_gamma entre " << elim_a << " et " << elim_b << " et |t| compris entre " << tlim_a << " et " << tlim_b << std::endl;

            Long64_t nentries = t1->GetEntries();
    
            for (Long64_t k = 0; k < nentries; k++) {
                t1->GetEntry(k);
    
                if (!electron || !positron) continue;

                // calcul du quadrivecteur k = p -p'

                Double_t E = Mp - proton->E();
                Double_t px = 0 - proton->Px();
                Double_t py = 0 - proton->Py();
                Double_t pz = 0 - proton->Pz();

                // calcul de |t| = abs(k²)
                Double_t t = abs(E*E - px*px - py*py - pz*pz);


                if (Q2 < Q2_max && abs(MMassBeam) < MMassBeam_max && Epho > elim_a && Epho < elim_b && t > tlim_a && t < tlim_b) {

                    double M = (*electron + *positron).M();

                    histograms[compteur_histo]->Fill(M);
    
                }
            }
            
            compteur_histo += 1;

        }
    }



    for (int k = 0; k < histograms.size(); k++) {
        // Créer un nom dynamique pour le TCanvas avec l'index k
        std::string canvasName = "c" + std::to_string(k);  
        TCanvas *c = new TCanvas(canvasName.c_str(), "Masse Invariante", 800, 600);
        
        histograms[k]->SetLineColor(kBlue);
        histograms[k]->GetXaxis()->SetTitle("Invariant mass [GeV]");
        histograms[k]->GetYaxis()->SetTitle("Number of events");
        histograms[k]->Draw();

    }

    // Création et ouverture d’un fichier ROOT en écriture
    TFile *outputFile = new TFile("Histos_M2Q.root", "RECREATE");

    for (int i = 0; i < histograms.size(); ++i) {
        if (histograms[i]) {
            histograms[i]->Write();  // Sauvegarde de chaque histogramme
        }
    }

    outputFile->Close();  // Ferme proprement le fichier
    std::cout << "Tous les histogrammes ont été sauvegardés dans Histos_M2Q.root" << std::endl;

}


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
    float bins = 200.0;

    int new_nb_bin = 12;


    TChain *t1 = new TChain("tree");  // "tree" doit être le nom EXACT du TTree dans les fichiers

    // Ajoute tous tes fichiers ROOT
    t1->Add("Data_pass2_fall2018_outbending.root");
    t1->Add("Data_pass2_spring2019_inbending.root");
    t1->Add("Data_pass2_fall2018_inbending.root");

    TLorentzVector *electron = nullptr;
    TLorentzVector *positron = nullptr;

    Float_t Q2, MMassBeam, Epho, elim_a, elim_b;
    Float_t counter = 0; 
    Float_t counterpho = 0;

    t1->SetBranchAddress("Electron", &electron);
    t1->SetBranchAddress("Positron", &positron);
    t1->SetBranchAddress("Q2", &Q2);
    t1->SetBranchAddress("MMassBeam", &MMassBeam);
    t1->SetBranchAddress("Epho", &Epho);

    std::vector<TH1F*> histograms;


    std::vector<double> liste = {2, 3.05333, 3.334, 3.566, 3.788, 4.006, 4.23333, 4.49, 4.76933, 5.1, 5.494, 6.076, 7.07933};

    // Affichage du contenu du vecteur
    std::cout << "Liste des éléments : ";
    for (float val : liste) std::cout << val << " ";
    std::cout << std::endl;

    // Initialisation des histogrammes
    for (int i = 0; i < new_nb_bin; i++) {
        std::cout << "Initialisation de l'histogramme numéro : " << i << std::endl;
        histograms.push_back(new TH1F(Form("Histogramme numero %d", i), "", bins, 0.5, 1.5));
    }

    std::cout << "Taille du vecteur d'histogrammes : " << histograms.size() << std::endl;

    for (int i = 0; i < new_nb_bin; i++) {
        if (i + 1 >= liste.size()) {
            std::cerr << "Erreur : accès hors limite à liste[" << i+1 << "]" << std::endl;
            continue;
        }

        elim_a = liste[i];
        elim_b = liste[i+1];

        std::cout << "Pour Epho compris entre " << elim_a << " et " << elim_b << std::endl;

        Long64_t nentries = t1->GetEntries();
        counter = 0;
        counterpho = 0;

        for (Long64_t j = 0; j < nentries; j++) {
            t1->GetEntry(j);

            if (!electron || !positron) continue;

            if (Q2 < Q2_max && abs(MMassBeam) < MMassBeam_max && Epho > elim_a && Epho < elim_b) {
                double M = (*electron + *positron).M();
                if (i < histograms.size()) {
                    counter++;
                    histograms[i]->Fill(M);
                } else {
                    std::cerr << "Erreur : indice i dépasse la taille de histograms" << std::endl;
                }
            }
        }

        for (Long64_t j = 0; j < nentries; j++) {
            t1->GetEntry(j);

            if (!electron || !positron) continue;

            if (Epho > elim_a && Epho < elim_b) {
                if (i < histograms.size()) {
                    counterpho++;
                } else {
                    std::cerr << "Erreur : indice i dépasse la taille de histograms" << std::endl;
                }
            }
        }

        std::cout << "Nombre de photon dans cet intervalle pre coupure = " << counterpho << std::endl;

        std::cout << "Nombre d'entries post-coupures dans cet intervalle = " << counter << std::endl;
        std::cout << "Histogramme terminé" << std::endl;
    }


    for (int k = 0; k < new_nb_bin; k++) {
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

    for (int i = 0; i < new_nb_bin; ++i) {
        if (histograms[i]) {
            histograms[i]->Write();  // Sauvegarde de chaque histogramme
        }
    }

    outputFile->Close();  // Ferme proprement le fichier
    std::cout << "Tous les histogrammes ont été sauvegardés dans Histos_M2Q.root" << std::endl;


}


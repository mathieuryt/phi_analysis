#include "TROOT.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include <fstream>
#include <string>
#include "TFile.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include <array>
#include <iostream>
#include <vector>
#include <tuple>
#include <TString.h>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <iostream>
#include <cmath>
#include "TStyle.h"
#include "TColor.h"


void HistoMinvKpKm() {

    
    gROOT->Reset();
    gSystem->Load("libPhysics.so");
    gSystem->Load("libCore.so");
    gSystem->Load("libTree.so");
    gSystem->Load("libRIO.so");
    gSystem->Load("libHist.so");
    gStyle->SetPalette(kBird); 

    // ROOT -> TREE
    bool isMC = false;
    TString data_adress;


    if(isMC == false){
        data_adress = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_data2/fall2019_outbending_neutronKpKm_bis_angle.root";
    }

    if(isMC == true){
        data_adress = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_data2/fall2019_outbending_neutronKpKm_mc_v2.root";
    }


    TFile *f = TFile::Open(data_adress);
    TTree *tree = (TTree*)f->Get("tree;1");


    TTree *tree_gen = (TTree*)f->Get("tree_gen;1"); // ligne a commenter en finction de si c'est MC ou pas pour prendre les GEN avec les REC-MC


    //VARIABLE DATAs && MONTECARLO
    double Q2, t, W;
    double t_missing_nucleon;
    double MinvKpKm, MM;
    double status_Kp, status_Km, status_n, status_el;
    double e_vx, e_vy, e_vz;
    double Kp_vx, Kp_vy, Kp_vz;
    double Km_vx, Km_vy, Km_vz;
    double n_vx, n_vy, n_vz;
    double angle_neutron_missnucl;
    TLorentzVector *Electron = nullptr;
    TLorentzVector *Neutron = nullptr;
    TLorentzVector *Kp = nullptr;
    TLorentzVector *Km = nullptr;
    TLorentzVector *Missing = nullptr;

    tree->SetBranchAddress("Electron", &Electron);
    tree->SetBranchAddress("Neutron", &Neutron);
    tree->SetBranchAddress("Kp", &Kp);
    tree->SetBranchAddress("Km", &Km);
    tree->SetBranchAddress("Missing", &Missing);
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("t", &t);
    tree->SetBranchAddress("t_missing_nucleon", &t_missing_nucleon);
    tree->SetBranchAddress("W", &W);
    tree->SetBranchAddress("MinvKpKm", &MinvKpKm);
    tree->SetBranchAddress("MissingMass", &MM);
    tree->SetBranchAddress("Status_Kp", &status_Kp);
    tree->SetBranchAddress("Status_Km", &status_Km);
    tree->SetBranchAddress("Status_n", &status_n);
    tree->SetBranchAddress("Status_el", &status_el);
    tree->SetBranchAddress("e_vx", &e_vx);
    tree->SetBranchAddress("e_vy", &e_vy);
    tree->SetBranchAddress("e_vz", &e_vz);
    tree->SetBranchAddress("n_vx", &n_vx);
    tree->SetBranchAddress("n_vy", &n_vy);
    tree->SetBranchAddress("n_vz", &n_vz);
    tree->SetBranchAddress("Km_vx", &Km_vx);
    tree->SetBranchAddress("Km_vy", &Km_vy);
    tree->SetBranchAddress("Km_vz", &Km_vz);
    tree->SetBranchAddress("Kp_vx", &Kp_vx);
    tree->SetBranchAddress("Kp_vy", &Kp_vy);
    tree->SetBranchAddress("Kp_vz", &Kp_vz);
    tree->SetBranchAddress("angle_neutron_missnucl", &angle_neutron_missnucl);


    //VARIABLE SPECIFIQUE MC

    double weight;
    double real_weight;
    double real_weight_correction;

    double weight_gen;
    double real_weight_gen;
    TLorentzVector *Electron_gen = nullptr;
    TLorentzVector *Neutron_gen = nullptr;
    TLorentzVector *Kp_gen = nullptr;
    TLorentzVector *Km_gen = nullptr;
    TLorentzVector *Missing_gen = nullptr;

    double Q2_gen, t_gen, W_gen;
    double MinvKpKm_gen, MM_gen;
    double e_vz_gen, Kp_vz_gen, Km_vz_gen, n_vz_gen;


    if(isMC == true){

        tree->SetBranchAddress("weight", &weight);
        tree->SetBranchAddress("real_weight", &real_weight);
        tree->SetBranchAddress("real_weight_correction", &real_weight_correction);

        tree_gen->SetBranchAddress("weight_gen", &weight_gen);
        tree_gen->SetBranchAddress("real_weight_gen", &real_weight_gen);

        tree_gen->SetBranchAddress("Electron_gen", &Electron_gen);
        tree_gen->SetBranchAddress("Neutron_gen", &Neutron_gen);
        tree_gen->SetBranchAddress("Kp_gen", &Kp_gen);
        tree_gen->SetBranchAddress("Km_gen", &Km_gen);
        tree_gen->SetBranchAddress("Missing_gen", &Missing_gen);

        tree_gen->SetBranchAddress("Q2_gen", &Q2_gen);
        tree_gen->SetBranchAddress("t_gen", &t_gen);
        tree_gen->SetBranchAddress("W_gen", &W_gen);

        tree_gen->SetBranchAddress("e_vz_gen", &e_vz_gen);
        tree_gen->SetBranchAddress("Kp_vz_gen", &Kp_vz_gen);
        tree_gen->SetBranchAddress("Km_vz_gen", &Km_vz_gen);
        tree_gen->SetBranchAddress("n_vz_gen", &n_vz_gen);

        tree_gen->SetBranchAddress("MinvKpKm_gen", &MinvKpKm_gen);
        tree_gen->SetBranchAddress("MissingMass_gen", &MM_gen);

    }


    //TRAITEMENT BORNES EN Q2
    std::vector<TH1F*> histograms;
    std::vector<TH1F*> histograms_Ngen;
    std::vector<double> bornesQ2;
    double valeur; //pr parcourir les fichiers

    std::string name_fichierQ2 = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_analysis/bornes_bins/bornes_Q2.txt";
    std::ifstream fichierQ2(name_fichierQ2);

    while (fichierQ2 >> valeur) {
         bornesQ2.push_back(valeur);
     }
 
    fichierQ2.close();

    std::cout << "Liste des bins en E_gamma : "; // Affichage du contenu du vecteur
    for (float val : bornesQ2) std::cout << val << " ";
    std::cout << std::endl;

    int nb_bin_Q2 = bornesQ2.size()-1;


    //TRAITEMENT BORNES EN t

    std::vector<std::vector<double>> liste_de_listes(bornesQ2.size()-1);

    size_t k = bornesQ2.size()-1; // Par exemple, on a 5 fichiers de 0 à 4

    for (size_t i = 0; i < k; ++i) {
        
        std::string name_fichiert = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_analysis/bornes_bins/bornes_t_" + std::to_string(i) + ".txt";
        std::ifstream fichiert(name_fichiert);

        if (fichiert) {
            double valeur;
            while (fichiert >> valeur) {
                liste_de_listes[i].push_back(valeur);
            }
            fichiert.close();
            std::cout << "Fichier " << name_fichiert << " lu avec succès." << std::endl;
        } else {
            std::cerr << "Erreur lors de l'ouverture du fichier : " << name_fichiert << std::endl;
        }
    }

    for (size_t i = 0; i < liste_de_listes.size(); ++i) { // Vérification : afficher le contenu
        std::cout << "liste des bornes bins en t pour le bin en Q2 numero " << i << std::endl;
        for (double val : liste_de_listes[i]) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }


    // Initialisation des histogrammes

    double bins = 18;
    double bornefill1, bornefill2;
    bornefill1 = 0.992;
    bornefill2 = 1.2;

    if (isMC){

        bornefill1 = 0.992;
        bornefill2 = 1.2;
        bins = 35;
    }

    for (int i = 0; i < nb_bin_Q2; i++) {
        for (int j =0; j<liste_de_listes[i].size()-1; j++) {

            std::cout << "Initialisation de l'histogramme pour le bin en Q2 " << i << " et le bin en t " << j << std::endl;
            histograms.push_back(new TH1F(Form("Histogramme_%d_%d", i, j), "", bins, bornefill1, bornefill2));
            histograms_Ngen.push_back(new TH1F(Form("Histogramme_Ngen_%d_%d", i, j), "", bins, 0.0, 25));

        }
    }

    std::cout << "Taille du vecteur d'histogrammes : " << histograms.size() << std::endl;


    //REMPLISSAGE histograms

    double compteur_histo = 0;

    for (int i = 0; i < nb_bin_Q2; i++) {
  
        double Q2_lim1 = bornesQ2[i];
        double Q2_lim2 = bornesQ2[i+1];

        for (int j = 0; j < liste_de_listes[i].size()-1; j++) {

            double tlim_a = liste_de_listes[i][j];
            double tlim_b = liste_de_listes[i][j+1];

            std::cout << "Remplissage REC pour Q2 entre " << Q2_lim1 << " et " << Q2_lim2 << " et |t| compris entre " << tlim_a << " et " << tlim_b << std::endl;

            Long64_t nentries = tree->GetEntries();
    
            for (Long64_t k = 0; k < nentries; k++) {
                tree->GetEntry(k);
    
                //if (!electron || !positron) continue;

                // calcul du quadrivecteur k = p -p'
                double M = (*Kp + *Km).M();


                if (angle_neutron_missnucl*180./3.14159 < 5 && Q2 > 1 && Missing->M2() < 0.5 && Missing->M2() > -0.5 && Electron->P()>2 && Neutron->Theta()*180./3.141592 > 4 && Neutron->P()>0.25 && status_Kp >= 2000 && status_Kp <= 2999 && status_Km >= 2000 && status_Km <= 2999 && Q2 > Q2_lim1 && Q2 < Q2_lim2 && t_missing_nucleon < tlim_a && t_missing_nucleon > tlim_b && M > 0.987) {

        

                    if(isMC == false){

                        histograms[compteur_histo]->Fill(M);

                    }


                    if(isMC == true){

                        histograms[compteur_histo]->Fill(M, real_weight_correction); // quand c'est MC on fit les rec avec le real_weight 

                    }
    
                }
            }

            //histograms[compteur_histo]->Sumw2();
            
            compteur_histo += 1;

        }
    }


    //REMPLISSAGE histograms GEN

    compteur_histo = 0;
    if(isMC){

        for (int i = 0; i < nb_bin_Q2; i++) {
  
            double Q2_lim1 = bornesQ2[i];
            double Q2_lim2 = bornesQ2[i+1];

            for (int j = 0; j < liste_de_listes[i].size()-1; j++) {

                double tlim_a = liste_de_listes[i][j];
                double tlim_b = liste_de_listes[i][j+1];

                std::cout << "Remplissage GEN pour Q2 entre " << Q2_lim1 << " et " << Q2_lim2 << " et |t| compris entre " << tlim_a << " et " << tlim_b << std::endl;

                histograms_Ngen[compteur_histo]->Sumw2();

                Long64_t nentries = tree_gen->GetEntries();
    
                for (Long64_t k = 0; k < nentries; k++) {

                        tree_gen->GetEntry(k);

                        if(Q2_gen > Q2_lim1 && Q2_gen < Q2_lim2 && t_gen < tlim_a && t_gen > tlim_b)

                            histograms_Ngen[compteur_histo]->Fill(1, real_weight_gen); // on obtient un histo des poids, il faudra juste sommer cet histo pour avoir somme des poid generer.

            
                }

                compteur_histo++;

            }

        }

    }




    for (int k = 0; k < histograms.size(); k++) {
        // Créer un nom dynamique pour le TCanvas avec l'index k
        std::string canvasName = "c" + std::to_string(k);  
        TCanvas *c = new TCanvas(canvasName.c_str(), "Masse Invariante REC", 800, 600);
        
        histograms[k]->SetLineColor(kBlue);
        histograms[k]->GetXaxis()->SetTitle("Invariant mass [GeV]");
        histograms[k]->GetYaxis()->SetTitle("Number of events");
        histograms[k]->Draw();

    }

    if(isMC){

        for (int k = 0; k < histograms_Ngen.size(); k++) {
            // Créer un nom dynamique pour le TCanvas avec l'index k
            std::string canvasName = "c2" + std::to_string(k);  
            TCanvas *c2 = new TCanvas(canvasName.c_str(), "Masse Invariante GEN", 800, 600);
        
            histograms_Ngen[k]->SetLineColor(kBlue);
            histograms_Ngen[k]->GetXaxis()->SetTitle("Invariant mass [GeV]");
            histograms_Ngen[k]->GetYaxis()->SetTitle("Number of events");
            histograms_Ngen[k]->Draw();

        }

    }

    // Création et ouverture d’un fichier ROOT en écriture
    TFile *outputFile = new TFile("Histos_MinvKpKm_data.root", "RECREATE");

    for (int i = 0; i < histograms.size(); ++i) {
        if (histograms[i]) {
            histograms[i]->Write();  // Sauvegarde de chaque histogramme
        }
    }

    if(isMC){

            for (int i = 0; i < histograms_Ngen.size(); ++i) {

                if (histograms_Ngen[i]) {
                    histograms_Ngen[i]->Write();  // Sauvegarde de chaque histogramme
                }

            }

    }

    outputFile->Close();  // Ferme proprement le fichier
    std::cout << "Tous les histogrammes ont été sauvegardés dans Histos_MinvKpKm_data.root" << std::endl;



}

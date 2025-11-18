#include "TH1F.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include <array>
#include <iostream>
#include <vector>
#include <tuple>
#include <TString.h>



void acceptance2() {

    gEnv->SetValue("Browser.Name", "TRootBrowser");

    double elim_a = 0;
    double elim_b = 1;
    double tlim_a = 0;
    double tlim_b = 1;
    double Mp = 0.938;

    double Q_in = 35.667;
    double Q_out = 32.451;
    double Q_s = 45.9;

    double bins_histo = 200;

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

    double nb_bins_Egamma = bornes.size()-1;

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
    for (int i = 0; i < nb_bins_Egamma; i++) {
        for (int j =0; j<liste_de_listes[i].size()-1; j++) {

            std::cout << "Initialisation de l'histogramme pour le bin en E_gamma " << i << " et le bin en |t| " << j << std::endl;
            histograms.push_back(new TH1F(Form("Histogramme_%d_%d", i, j), "", bins_histo, 0.8, 1.2));

        }
    }

    std::cout << "Taille du vecteur d'histogrammes : " << histograms.size() << std::endl;

    double compteur_histo = 0;

    gEnv->SetValue("Browser.Name", "TRootBrowser");
    TFile *file_in = TFile::Open("outputTCS_Phi_inbending.root");
    TTree *t1_in = (TTree*)file_in->Get("tree;1");
    TTree *t2_in = (TTree*)file_in->Get("tree_Gen;3");

    TLorentzVector *electron = nullptr;
    TLorentzVector *positron = nullptr;
    TLorentzVector *proton = nullptr;
 
    Float_t Q2, MMassBeam, Ephovar_rec, weightvar_rec;

    t1_in->SetBranchAddress("Electron", &electron);
    t1_in->SetBranchAddress("Positron", &positron);
    t1_in->SetBranchAddress("Proton", &proton);
    t1_in->SetBranchAddress("Q2", &Q2);
    t1_in->SetBranchAddress("MMassBeam", &MMassBeam);
    t1_in->SetBranchAddress("Epho", &Ephovar_rec);
    t1_in->SetBranchAddress("weight", &weightvar_rec);

    Float_t Ephovar_gen, weightvar_gen, t_gen;

    t2_in->SetBranchAddress("t_Gen", &t_gen);
    t2_in->SetBranchAddress("Epho_Gen", &Ephovar_gen);
    t2_in->SetBranchAddress("weight", &weightvar_gen);

    TFile *file_out = TFile::Open("outputTCS_Phi_outbending.root");
    TTree *t1_out = (TTree*)file_out->Get("tree;1");
    TTree *t2_out = (TTree*)file_out->Get("tree_Gen;3");

    TLorentzVector *electron_out = nullptr;
    TLorentzVector *positron_out = nullptr;
    TLorentzVector *proton_out = nullptr;
 
    Float_t Q2_out, MMassBeam_out, Ephovar_rec_out, weightvar_rec_out;

    t1_out->SetBranchAddress("Electron", &electron_out);
    t1_out->SetBranchAddress("Positron", &positron_out);
    t1_out->SetBranchAddress("Proton", &proton_out);
    t1_out->SetBranchAddress("Q2", &Q2_out);
    t1_out->SetBranchAddress("MMassBeam", &MMassBeam_out);
    t1_out->SetBranchAddress("Epho", &Ephovar_rec_out);
    t1_out->SetBranchAddress("weight", &weightvar_rec_out);

    Float_t Ephovar_gen_out, weightvar_gen_out, t_gen_out;

    t2_out->SetBranchAddress("t_Gen", &t_gen_out);
    t2_out->SetBranchAddress("Epho_Gen", &Ephovar_gen_out);
    t2_out->SetBranchAddress("weight", &weightvar_gen_out);

    for (int i = 0; i < nb_bins_Egamma; i++) {
  
        elim_a = bornes[i];
        elim_b = bornes[i+1];

        for (int j = 0; j < liste_de_listes[i].size()-1; j++) {

            tlim_a = liste_de_listes[i][j];
            tlim_b = liste_de_listes[i][j+1];

            std::cout << "Remplissage REC pour E_gamma entre " << elim_a << " et " << elim_b << " et |t| compris entre " << tlim_a << " et " << tlim_b << std::endl;

            Long64_t nentries = t1_in->GetEntries();
    
            for (Long64_t k = 0; k < nentries; k++) {
                t1_in->GetEntry(k);
    
                if (!electron || !positron) continue;

                // calcul du quadrivecteur k = p -p'

                Double_t E = Mp - proton->E();
                Double_t px = 0 - proton->Px();
                Double_t py = 0 - proton->Py();
                Double_t pz = 0 - proton->Pz();

                // calcul de |t| = abs(k²)
                Double_t t = abs(E*E - px*px - py*py - pz*pz);

                if (Ephovar_rec > elim_a && Ephovar_rec < elim_b && t < tlim_b && t > tlim_a) {

                    double M = (*electron + *positron).M();
                    histograms[compteur_histo]->Fill(M, weightvar_rec*(Q_in+Q_s));
    
                }

            }

            Long64_t nentries_bis = t1_out->GetEntries();


            for (Long64_t k = 0; k < nentries_bis; k++) {

                t1_out->GetEntry(k);
    
                if (!electron_out || !positron_out) continue;

                // calcul du quadrivecteur k = p -p'

                Double_t E = Mp - proton_out->E();
                Double_t px = 0 - proton_out->Px();
                Double_t py = 0 - proton_out->Py();
                Double_t pz = 0 - proton_out->Pz();

                // calcul de |t| = abs(k²)
                Double_t t = abs(E*E - px*px - py*py - pz*pz);

                if (Ephovar_rec_out > elim_a && Ephovar_rec_out < elim_b && t < tlim_b && t > tlim_a) {

                    double M = (*electron_out + *positron_out).M();
                    histograms[compteur_histo]->Fill(M, weightvar_rec_out*Q_out);
    
                }

            }

            
            compteur_histo += 1;

        }
    }

    double compteur_liste = 0;

    std::vector<double> Nphi_GEN(histograms.size());

    for (int i = 0; i < nb_bins_Egamma; i++) {
  
        elim_a = bornes[i];
        elim_b = bornes[i+1];

        for (int j = 0; j < liste_de_listes[i].size()-1; j++) {

            tlim_a = liste_de_listes[i][j];
            tlim_b = liste_de_listes[i][j+1];

            std::cout << "Remplissage GEN pour E_gamma entre " << elim_a << " et " << elim_b << " et |t| compris entre " << tlim_a << " et " << tlim_b << std::endl;

            Long64_t nentries = t2_in->GetEntries();
    
            for (Long64_t k = 0; k < nentries; k++) {

                t2_in->GetEntry(k);


                if (Ephovar_gen > elim_a && Ephovar_gen < elim_b && abs(t_gen) < tlim_b && abs(t_gen) > tlim_a) {
    
                    Nphi_GEN[compteur_liste] = Nphi_GEN[compteur_liste] + weightvar_gen*(Q_in+Q_s);

                }
        
            }

            Long64_t nentries_bis = t2_out->GetEntries();


            for (Long64_t k = 0; k < nentries_bis; k++) {

                t2_out->GetEntry(k);


                if (Ephovar_gen_out > elim_a && Ephovar_gen_out < elim_b && abs(t_gen_out) < tlim_b && abs(t_gen_out) > tlim_a) {
    
                    Nphi_GEN[compteur_liste] = Nphi_GEN[compteur_liste] + weightvar_gen_out*Q_out;

                }

            }
            
            compteur_liste += 1;

        }
    }

    std::cout << "Nombre de phi GEN: ";
    for (float val : Nphi_GEN) std::cout << val << " ";
    std::cout << std::endl;

    std::ofstream fichier2("nb_phi_GEN.txt");

    if (fichier2) {
        for (float valeur : Nphi_GEN) {
            fichier2 << valeur << std::endl;
        }
        fichier2.close();
        std::cout << "Fichier nb de phi GEN écrit avec succès !" << std::endl;
    } else {
        std::cerr << "Erreur lors de l'ouverture du fichier !" << std::endl;
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

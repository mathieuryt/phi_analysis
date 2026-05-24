#ifndef FILLBINNING_H
#define FILLBINNING_H

#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TH2D.h"

#include <vector>

bool frontier(double Pkp, double Pkm){

        const double low1  = (23.0/36.0)*Pkp - (5.0/36.0);
        const double low2  = -2.0*Pkp + 2.0;
        const double high1 = (68.0/47.0)*Pkp + (125.0/235.0);
        const double high2 = (-8.0/9.0)*Pkp + (66.0/9.0);

        return (Pkm < high1 && Pkm < high2 && Pkm > low1 && Pkm > low2);

}


void fillbinning(const char* chemin, bool isMC)
{

    if(isMC){

        std::cout << "Fillbininig MC...... " << std::endl;

    }

    if(isMC == false){

        std::cout << "Fillbininig DATA...... " << std::endl;

    }
    

 

    // Ouvre le fichier ROOT d'entrée
    TFile* fin = TFile::Open(chemin, "READ");

    // ===== A ADAPTER : nom du TTree =====
    TTree* tree = (TTree*)fin->Get("tree");

    TLorentzVector *Kp = nullptr;
    TLorentzVector *Km = nullptr;
    TLorentzVector Kp_smear;
    TLorentzVector Km_smear;
    TLorentzVector *Missing = nullptr;
    TLorentzVector *Electron = nullptr;
    TLorentzVector *Neutron = nullptr;

    TVector3 Kp_smear_vect;
    TVector3 Km_smear_vect;

    double angle_neutron_missnucl, status_Km, status_Kp, Q2, real_weight_effcorrection;
    double Pkp_smear,Pkm_smear, M_smear;
    double bins = 28;
    double bornefill1, bornefill2;
    bornefill1 = 0.992;
    bornefill2 = 1.1;

    std::vector<double> bornesKp;
    std::vector<TH1F*> histograms;
    std::vector<TH2F*> histograms2D;
    std::vector<TTree*> trees;


    TLorentzVector Kp_out, Km_out;

    

    tree->SetBranchAddress("Kp", &Kp);
    tree->SetBranchAddress("Km", &Km);
    tree->SetBranchAddress("Missing", &Missing);
    tree->SetBranchAddress("Electron", &Electron);
    tree->SetBranchAddress("Neutron", &Neutron);
    tree->SetBranchAddress("Status_Kp", &status_Kp);
    tree->SetBranchAddress("Status_Km", &status_Km);
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("angle_neutron_missnucl", &angle_neutron_missnucl);

    if(isMC){

        tree->SetBranchAddress("real_weight_effcorrection", &real_weight_effcorrection);

    }


    std::string name_fichierQ2 = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_analysis/acceptance/binningPKpKm/bornesKp.txt";
    std::ifstream fichierQ2(name_fichierQ2);

    double valeur;
    while (fichierQ2 >> valeur) {
         bornesKp.push_back(valeur);
     }
 
    fichierQ2.close();

    std::vector<std::vector<double>> liste_de_listes_Km(bornesKp.size()-1);

    size_t k = bornesKp.size()-1; // 5 nombre de le fichier donc 4 bins

    for (size_t i = 0; i < k; ++i) {
        
        std::string name_fichiert = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_analysis/acceptance/binningPKpKm/bornesKm_" + std::to_string(i) + ".txt";
        std::ifstream fichiert(name_fichiert);

        if (fichiert) {
            double valeur;
            while (fichiert >> valeur) {
                liste_de_listes_Km[i].push_back(valeur);
            }
            fichiert.close();
            //std::cout << "Fichier " << name_fichiert << " lu avec succès." << std::endl;
        } else {
            std::cerr << "Erreur lors de l'ouverture du fichier : " << name_fichiert << std::endl;
        }
    }

    //std::cout << "Liste des bornes en Kp (phase de FILL): "; // Affichage du contenu du vecteur
    //for (float val : bornesKp) std::cout << val << " ";
    //std::cout << std::endl;

    //for (size_t i = 0; i < liste_de_listes_Km.size(); ++i) { // Vérification : afficher le contenu
        //std::cout << "liste des bornes bins en Km pour le bin en Kp numero (phase de FILL) " << i << std::endl;
        //for (double val : liste_de_listes_Km[i]) {
            //std::cout << val << " ";
        //}
        //std::cout << std::endl;
    //}

    //Inititalisation des histogrames

    TFile* fout = nullptr;

    if (isMC) {
        fout = TFile::Open("TreeKpKm_MC.root", "RECREATE");
    } else {
        fout = TFile::Open("HistoMinvbinP_data.root", "RECREATE");
    }

    fout->cd();

    for (int i = 0; i < k; i++) {

        for (int j = 0; j < liste_de_listes_Km[i].size()-1; j++) {

            TTree* t = new TTree(Form("Tree_%d_%d", i, j), "Kp Km per bin");

            t->Branch("Kp", &Kp_out);
            t->Branch("Km", &Km_out);

            trees.push_back(t);
        }
    }

    
    for (int i = 0; i < k; i++) { //k is bornesKp.size()-1; 
        for (int j =0; j<liste_de_listes_Km[i].size()-1; j++) {

            //std::cout << "Initialisation de l'histogramme pour le bin en Kp " << i << " et le bin en Km " << j << std::endl;
            histograms.push_back(new TH1F(Form("Histogramme_%d_%d", i, j), "", bins, bornefill1, bornefill2));

        }
    }

    for (int i = 0; i < k; i++) { //k is bornesKp.size()-1; 
        for (int j =0; j<liste_de_listes_Km[i].size()-1; j++) {

            //std::cout << "Initialisation de l'histogramme 2D pour le bin en Kp " << i << " et le bin en Km " << j << std::endl;
            histograms2D.push_back(new TH2F(Form("Histogramme2D_%d_%d", i, j), "",
                              150, 0, 8,
                              150, 0, 8));

        }
    }
    

    double compteur_histo = 0;
    double Kp_lim1, Kp_lim2;
    double Km_lim1, Km_lim2;
    double sigmaSmearKp;
    double sigmaSmearKm;

    //sigmaSmearKp = 0.05;

    for (int i = 0; i < k; i++) {

        double Kp_lim1 = bornesKp[i];
        double Kp_lim2 = bornesKp[i+1];

        for (int j = 0; j < liste_de_listes_Km[i].size()-1; j++) {

            double Km_lim1 = liste_de_listes_Km[i][j];
            double Km_lim2 = liste_de_listes_Km[i][j+1];

            //std::cout << "Remplissage REC pour Kp entre " << Kp_lim1 << " et " << Kp_lim2 << " et Km compris entre " << Km_lim1 << " et " << Km_lim2 << std::endl;

            Long64_t nentries = tree->GetEntries();

            if(isMC){
                nentries = nentries/1.0;
            }
    
            for (Long64_t ent = 0; ent < nentries; ent++) {

                tree->GetEntry(ent);
    
                //if (!electron || !positron) continue;

                // calcul du quadrivecteur k = p -p'
                double M = (*Kp + *Km).M();


                if (angle_neutron_missnucl*180./3.14159 < 5 && Q2 > 1 && Missing->M2() < 0.5 && Missing->M2() > -0.5 && Electron->P()>2 && Neutron->Theta()*180./3.141592 > 4 && Neutron->P()>0.25 && status_Kp >= 2000 && status_Kp <= 2999 && status_Km >= 2000 && status_Km <= 2999 && Kp->P() > Kp_lim1 && Kp->P() < Kp_lim2 && Km->P() > Km_lim1 && Km->P() < Km_lim2 && M > 0.987) {

    
                    if(isMC == false){

                        histograms[compteur_histo]->Fill(M);
                        histograms2D[compteur_histo]->Fill(Kp->P(), Km->P());

                    }


                    if(isMC == true){

                        //sigmaSmearKp = Sigma(Kp->P(), 0, 1.2, 1.6, 2.2, 6.0, A1, B1, C1, D1);
                        //sigmaSmearKm = Sigma(Km->P(), 0, 1.2, 1.6, 2.2, 6.0, A2, B2, C2, D2);

                        //Pkp_smear = Kp->P() * rng.Gaus(1.0, sigmaSmearKp);
                        //Pkm_smear = Km->P() * rng.Gaus(1.0, sigmaSmearKm);
                        
                        //Kp_smear_vect.SetMagThetaPhi(Pkp_smear, Kp->Theta(), Kp->Phi());
                        //Km_smear_vect.SetMagThetaPhi(Pkm_smear, Km->Theta(), Km->Phi());

                        //Kp_smear.SetVectM(Kp_smear_vect, m_kaons);
                        //Km_smear.SetVectM(Km_smear_vect, m_kaons);        

                        //M_smear = (Kp_smear + Km_smear).M();
                   
                        //histograms[compteur_histo]->Fill(M_smear, real_weight_effcorrection); // quand c'est MC on fit les rec avec le real_weight 

                        Kp_out = *Kp;
                        Km_out = *Km;

                        trees[compteur_histo]->Fill();

                        histograms2D[compteur_histo]->Fill(Kp->P(), Km->P(), real_weight_effcorrection);

                    }
    
                }
            }

            compteur_histo += 1;

        }

    }


    for (auto t : trees) {
        t->Write();
    }

    if(isMC == false){ // on inscrit que les data en Minv car pour les MC on garde les tree. 

        for (int k = 0; k < histograms.size(); k++) {

            histograms[k]->GetXaxis()->SetTitle("Invariant mass [GeV]");
            histograms[k]->GetYaxis()->SetTitle("Number of events");
            histograms[k]->Write();

        }

    }

    for (int k = 0; k < histograms2D.size(); k++) {

        histograms2D[k]->GetXaxis()->SetTitle("K+ momentum");
        histograms2D[k]->GetYaxis()->SetTitle("K- momentum");
        histograms2D[k]->Write();

    }

    fout->Close();
    fin->Close();

}

#endif
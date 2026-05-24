#ifndef CREATEROOT_H
#define CREATEROOT_H

#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TH2D.h"

void createroot(const char* chemin, double s, bool isMC)
{
    // s est inutilisé pour l'instant
    (void)s;

    // Ouvre le fichier ROOT d'entrée
    TFile* fin = TFile::Open(chemin, "READ");

    // ===== A ADAPTER : nom du TTree =====
    TTree* tree = (TTree*)fin->Get("tree");


    // ===== A ADAPTER : noms des branches =====
    // Exemple : composantes impulsion du K+ et du K-
    TLorentzVector *Kp = nullptr;
    TLorentzVector *Km = nullptr;
    TLorentzVector *Missing = nullptr;
    TLorentzVector *Electron = nullptr;
    TLorentzVector *Neutron = nullptr;

    double angle_neutron_missnucl, status_Km, status_Kp, Q2, real_weight_effcorrection;

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


    // Histogramme de masse invariante
    TH1D* hMinv = new TH1D("hMinv", "Masse invariante K^{+}K^{-};M_{K^{+}K^{-}} (GeV);Entries", 200, 0.9, 1.3);
    TH1D* hPkp = new TH1D("hPkp", "Momemtum K^{+} ;Momemtum K^{+} (GeV);Entries", 200, 0.0, 8);
    TH1D* hPkm = new TH1D("hPkm", "Momemtum K^{-} ;Momemtum K^{-} (GeV);Entries", 200, 0.0, 8);
    TH2D* hPkp_vs_Pkm = new TH2D("hPkp_vs_Pkm", "P(K^{+}) vs P(K^{-}); P(K^{+}) (GeV);P(K^{-}) (GeV)", 200, 0.0, 8.0, 200, 0.0, 8.0);
    TH2D* hPkp_vs_Pkm_v2 = new TH2D("hPkp_vs_Pkm_v2", "P(K^{+}) vs P(K^{-}) with cut on Minv; P(K^{+}) (GeV);P(K^{-}) (GeV)", 200, 0.0, 8.0, 200, 0.0, 8.0);

    Long64_t nentries = tree->GetEntries();

    if(isMC == false){

        for (Long64_t i = 0; i < nentries; i++) {

            tree->GetEntry(i);

            double Minv = (*Kp + *Km).M();

            if(angle_neutron_missnucl*180./3.14159 < 5 && Q2 > 1 && Missing->M2() < 0.5 && Missing->M2() > -0.5 && Electron->P()>2 && Neutron->Theta()*180./3.141592 > 4 && Neutron->P()>0.25 && status_Kp >= 2000 && status_Kp <= 2999 && status_Km >= 2000 && status_Km <= 2999){

                hMinv->Fill(Minv);
                hPkp->Fill(Kp->P());
                hPkm->Fill(Km->P());
                hPkp_vs_Pkm->Fill(Kp->P(), Km->P());

            }

            if(angle_neutron_missnucl*180./3.14159 < 5 && Q2 > 1 && Missing->M2() < 0.5 && Missing->M2() > -0.5 && Electron->P()>2 && Neutron->Theta()*180./3.141592 > 4 && Neutron->P()>0.25 && status_Kp >= 2000 && status_Kp <= 2999 && status_Km >= 2000 && status_Km <= 2999 && Minv < 1.04 && Minv > 1.00){
                hPkp_vs_Pkm_v2->Fill(Kp->P(), Km->P());
            }

        }
    }

    if(isMC == true){

        for (Long64_t i = 0; i < nentries; i++) {

            tree->GetEntry(i);

            double Minv = (*Kp + *Km).M();

            if(angle_neutron_missnucl*180./3.14159 < 5 && Q2 > 1 && Missing->M2() < 0.5 && Missing->M2() > -0.5 && Electron->P()>2 && Neutron->Theta()*180./3.141592 > 4 && Neutron->P()>0.25 && status_Kp >= 2000 && status_Kp <= 2999 && status_Km >= 2000 && status_Km <= 2999){

                hMinv->Fill(Minv, real_weight_effcorrection);
                hPkp->Fill(Kp->P(), real_weight_effcorrection);
                hPkm->Fill(Km->P(), real_weight_effcorrection);
                hPkp_vs_Pkm->Fill(Kp->P(), Km->P(), real_weight_effcorrection);

            }

        }
    }


    // Sauvegarde dans un nouveau fichier ROOT
    
    if(isMC){

        TFile* fout = TFile::Open("HistoMinv_mc.root", "RECREATE");
        hMinv->Write();
        hPkp->Write();
        hPkm->Write();
        hPkp_vs_Pkm->Write();

        fout->Close();
        fin->Close();

    }else{

        TFile* fout = TFile::Open("HistoMinv_data.root", "RECREATE");
        hMinv->Write();
        hPkp->Write();
        hPkm->Write();
        hPkp_vs_Pkm->Write();
        hPkp_vs_Pkm_v2->Write();

        fout->Close();
        fin->Close();
    }

    if(isMC){

        std::cout << "Histogramme sauvegarde dans HistoMinv_mc.root" << std::endl;

    }else{

        std::cout << "Histogramme sauvegarde dans HistoMinv_data.root" << std::endl;

    }

}

#endif
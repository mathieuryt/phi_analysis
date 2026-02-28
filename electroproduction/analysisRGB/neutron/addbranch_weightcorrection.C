#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TCut.h"
#include "TStyle.h"
#include "TLine.h"
#include "TString.h"
#include <vector>
#include <iostream>
#include "TLatex.h"
#include "TLorentzVector.h"

void addbranch_weightcorrection() {

    // ------------------------------
    // Ouvrir le fichier TTree MC
    // ------------------------------
    TString filename_tree = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_data2/fall2019_outbending_neutronKpKm_mc_v2.root";
    TFile *f_tree = TFile::Open(filename_tree, "UPDATE");
    TTree *t = (TTree*)f_tree->Get("tree");

    // ------------------------------
    // Ouvrir le fichier avec h_ratio_p
    // ------------------------------
    TString filename_ratio = "ratio_efficiency.root"; // ton fichier ROOT sauvegardé
    TFile *f_ratio = TFile::Open(filename_ratio, "READ");
    TH1F *h_ratio_p = (TH1F*)f_ratio->Get("h_ratio_p");

    if (!h_ratio_p) {
        std::cout << "Erreur: h_ratio_p non trouvé dans " << filename_ratio << std::endl;
        return;
    }

    h_ratio_p->SetDirectory(0);
    f_ratio->Close();

    // ------------------------------
    // Variables existantes
    // ------------------------------
    TLorentzVector *Electron = nullptr;
    TLorentzVector *Kp = nullptr;
    TLorentzVector *Km = nullptr;
    TLorentzVector *Missing_nucleon = nullptr;

    double weight;
    double real_weight;

    t->SetBranchAddress("Electron", &Electron);
    t->SetBranchAddress("Kp", &Kp);
    t->SetBranchAddress("Km", &Km);
    t->SetBranchAddress("Missing_nucleon", &Missing_nucleon);
    t->SetBranchAddress("weight", &weight);
    t->SetBranchAddress("real_weight", &real_weight);


    // ------------------------------
    // Déclarer la nouvelle branch
    // ------------------------------

    double real_weight_correction;

    TBranch *breal_weight_correction = t->Branch("real_weight_correction", &real_weight_correction, "real_weight_correction/D");

    // ------------------------------
    // Boucle sur le TTree


    Long64_t nentries = t->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);

        double P_missing_nucleon = Missing_nucleon->P();

        // valeur de correction dans h_ratio_p
        int bin = h_ratio_p->FindBin(P_missing_nucleon);
        double corr = 1.0;
        if (bin >= 1 && bin <= h_ratio_p->GetNbinsX())
            corr = h_ratio_p->GetBinContent(bin);

        // poids corrigé
        real_weight_correction = real_weight * corr;
        if(i % 500 == 0){
            std::cout<< "P miss = " << P_missing_nucleon << " and correction = " << corr << std::endl;
        
        }

        // remplir la branch
        breal_weight_correction->Fill();

    }

    // ------------------------------
    //  Écrire et fermer
    // ------------------------------
    f_tree->cd(); 
    t->Write("", TObject::kOverwrite);
    f_tree->Close();
    f_ratio->Close();

    std::cout << "Branch weight_mc_correction ajoutée avec succès !" << std::endl;


}

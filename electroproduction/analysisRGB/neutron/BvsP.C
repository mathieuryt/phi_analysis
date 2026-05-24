#include "TROOT.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include "TSystem.h"

#include <fstream>
#include <string>
#include <array>
#include <iostream>
#include <vector>
#include <tuple>
#include <iomanip>
#include <algorithm>
#include <cmath>

void BvsP() {

    gROOT->Reset();
    gStyle->SetOptStat(0);
    gSystem->Load("libPhysics.so");
    gSystem->Load("libCore.so");
    gSystem->Load("libTree.so");
    gSystem->Load("libRIO.so");
    gSystem->Load("libHist.so");
    gStyle->SetPalette(kBird);

    TString data_adress;
    data_adress = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_data2/fall2019_outbending_neutronKpKm_data_Bvsp.root";

    TFile *f = TFile::Open(data_adress);

    TTree *tree_electron = (TTree*)f->Get("tree_electron;1");
    TTree *tree_Kp = (TTree*)f->Get("tree_Kp;1");
    TTree *tree_Km = (TTree*)f->Get("tree_Km;1");
    TTree *tree_pip = (TTree*)f->Get("tree_pip;1");
    TTree *tree_pim = (TTree*)f->Get("tree_pim;1");
    TTree *tree_proton = (TTree*)f->Get("tree_proton;1");

    // VARIABLES DATA
    double Status_Kp, Status_Km, Status_pip, Status_pim, Status_proton;

    TLorentzVector *Electron = nullptr;
    TLorentzVector *Kp = nullptr;
    TLorentzVector *Km = nullptr;
    TLorentzVector *pip = nullptr;
    TLorentzVector *pim = nullptr;
    TLorentzVector *proton = nullptr;

    tree_electron->SetBranchAddress("Electron", &Electron);
    tree_Kp->SetBranchAddress("Kp", &Kp);
    tree_Km->SetBranchAddress("Km", &Km);
    tree_pip->SetBranchAddress("pip", &pip);
    tree_pim->SetBranchAddress("pim", &pim);
    tree_proton->SetBranchAddress("proton", &proton);

    tree_Kp->SetBranchAddress("Status_Kp", &Status_Kp);
    tree_Km->SetBranchAddress("Status_Km", &Status_Km);
    tree_pip->SetBranchAddress("Status_pip", &Status_pip);
    tree_pim->SetBranchAddress("Status_pim", &Status_pim);
    tree_proton->SetBranchAddress("Status_proton", &Status_proton);

    // Histogrammes 2D : beta vs p
    TH2F *h_BvsP_pos = new TH2F(
        "h_BvsP_pos",
        "#beta vs p for positive particles (FD);p (GeV/c);#beta",
        300, 0.0, 5.0,
        300, 0.0, 1.2
    );

    TH2F *h_BvsP_neg = new TH2F(
        "h_BvsP_neg",
        "#beta vs p for negative particles (FD);p (GeV/c);#beta",
        300, 0.0, 5.0,
        300, 0.0, 1.2
    );

    //proton

    Long64_t nentries = tree_proton->GetEntries();

    for (Long64_t k = 0; k < nentries; k++) {
        tree_proton->GetEntry(k);

        double B_proton = proton->Beta();


        double p_proton = proton->P();


        // Tous les hadrons dans le Forward Detector (FD)
        if (Status_proton >= 2000 && Status_proton <= 2999) {

            // Particules positives : proton, K+, pi+
            h_BvsP_pos->Fill(p_proton, B_proton);

        }
    }


    //Kp

    nentries = tree_Kp->GetEntries();

    for (Long64_t k = 0; k < nentries; k++) {
        tree_Kp->GetEntry(k);

        double B_Kp = Kp->Beta();
        double p_Kp = Kp->P();


        // Tous les hadrons dans le Forward Detector (FD)
        if (Status_Kp >= 2000 && Status_Kp <= 2999) {

            // Particules positives : proton, K+, pi+
            h_BvsP_pos->Fill(p_Kp, B_Kp);

        }
    }

    //pip

    nentries = tree_pip->GetEntries();

    for (Long64_t k = 0; k < nentries; k++) {
        tree_pip->GetEntry(k);

        double B_pip = pip->Beta();
        double p_pip = pip->P();


        // Tous les hadrons dans le Forward Detector (FD)
        if (Status_pip >= 2000 && Status_pip <= 2999) {

            // Particules positives : proton, K+, pi+
            h_BvsP_pos->Fill(p_pip, B_pip);

        }
    }

    //Km

    nentries = tree_Km->GetEntries();

    for (Long64_t k = 0; k < nentries; k++) {
        tree_Km->GetEntry(k);

        double B_Km = Km->Beta();
        double p_Km = Km->P();


        // Tous les hadrons dans le Forward Detector (FD)
        if (Status_Km >= 2000 && Status_Km <= 2999) {

            // Particules positives : proton, K+, pi+
            h_BvsP_neg->Fill(p_Km, B_Km);

        }
    }

    //pim

    nentries = tree_pim->GetEntries();

    for (Long64_t k = 0; k < nentries; k++) {
        tree_pim->GetEntry(k);

        double B_pim = pim->Beta();
        double p_pim = pim->P();


        // Tous les hadrons dans le Forward Detector (FD)
        if (Status_pim >= 2000 && Status_pim <= 2999) {

            // Particules positives : proton, K+, pi+
            h_BvsP_neg->Fill(p_pim, B_pim);

        }
    }


    // Canvas et sauvegarde PDF
    TCanvas *c1 = new TCanvas("c1", "Beta vs Momentum", 1000, 800);
    c1->Print("BvsP.pdf[");

    h_BvsP_pos->Draw("COLZ");
    c1->Print("BvsP.pdf");

    h_BvsP_neg->Draw("COLZ");
    c1->Print("BvsP.pdf");

    c1->Print("BvsP.pdf]");

    // Sauvegarde optionnelle dans un fichier ROOT
    TFile *outFile = new TFile("BvsP_histos.root", "RECREATE");
    h_BvsP_pos->Write();
    h_BvsP_neg->Write();
    outFile->Close();

    f->Close();

    std::cout << "Termine. Les plots ont ete sauvegardes dans BvsP.pdf" << std::endl;
}
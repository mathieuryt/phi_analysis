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


void momentum_corr_kaons() {

    
    gROOT->Reset();
    gStyle->SetOptStat(0);
    gSystem->Load("libPhysics.so");
    gSystem->Load("libCore.so");
    gSystem->Load("libTree.so");
    gSystem->Load("libRIO.so");
    gSystem->Load("libHist.so");
    gStyle->SetPalette(kBird); 

    // ROOT -> TREE
    bool isMC = true;
    TString data_adress;


    if(isMC == false){
        data_adress = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_data2/fall2019_outbending_neutronKpKm_bis_angle.root";
    }

    if(isMC == true){
        data_adress = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_data2/fall2019_outbending_neutronKpKm_mc_v5.root";
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
    TLorentzVector *el_gen_associated = nullptr;
    TLorentzVector *Neutron = nullptr;
    TLorentzVector *Kp = nullptr;
    TLorentzVector *Km = nullptr;
    TLorentzVector *Kp_gen_associated = nullptr;
    TLorentzVector *Km_gen_associated = nullptr;
    TLorentzVector *Missing = nullptr;

    tree->SetBranchAddress("Electron", &Electron);
    tree->SetBranchAddress("el_gen_associated", &el_gen_associated);
    tree->SetBranchAddress("Neutron", &Neutron);
    tree->SetBranchAddress("Kp", &Kp);
    tree->SetBranchAddress("Km", &Km);
    tree->SetBranchAddress("Kp_gen_associated", &Kp_gen_associated);
    tree->SetBranchAddress("Km_gen_associated", &Km_gen_associated);
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
    double real_weight_effcorrection;
    double real_weight_SDME_effcorrection;

    double weight_gen;
    double real_weight_gen_SDME;
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
        tree->SetBranchAddress("real_weight_SDME_effcorrection", &real_weight_SDME_effcorrection); //le seul qui sert ici 
        //tree->SetBranchAddress("real_weight_effcorrection", &real_weight_effcorrection); //ou lui

        tree_gen->SetBranchAddress("weight_gen", &weight_gen);
        tree_gen->SetBranchAddress("real_weight_gen_SDME", &real_weight_gen_SDME); //same 
        tree_gen->SetBranchAddress("real_weight_gen", &real_weight_gen); //same 

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

    //creer une liste p_rec
    // creer un histo 1d avec pour bin en P_rec la liste precedante

    TProfile *prof = new TProfile(
    "prof",
    ";p_{rec} (GeV);#LT100*(p_{rec} - p_{gen}) / p_{gen}#GT (\%)",
    10, 0.5, 4.5
    );

    TProfile *prof_Km = new TProfile(
    "prof_Km",
    ";p_{rec} (GeV);#LT100*(p_{rec} - p_{gen}) / p_{gen}#GT (\%)",
    10, 0.5, 4.5
    );

    TProfile *prof_e = new TProfile(
    "prof_e",
    ";p_{rec} (GeV);#LT100*(p_{rec} - p_{gen}) / p_{gen}#GT (\%)",
    10, 2, 8
    );

    TH2D *h2_ratio_vs_prec = new TH2D(
    "h2_ratio_vs_prec",
    "ratio vs p_{rec};p_{rec} (GeV);ratio",
    50, 0.5, 4.5,     // bins en p_rec
    50, -4.5, 2.0    // bins en ratio, à adapter
    );



    Long64_t nentries = tree->GetEntries();
    
    for (Long64_t k = 0; k < nentries; k++) {
        
        tree->GetEntry(k);

        double p_rec = Kp->P();
        double p_gen = Kp_gen_associated->P();
        double ratio = (p_rec - p_gen)/p_gen;

        double p_rec_Km = Km->P();
        double p_gen_Km = Km_gen_associated->P();
        double ratio_Km = (p_rec_Km - p_gen_Km)/p_gen_Km;

        double p_rec_e = Electron->P();
        double p_gen_e = el_gen_associated->P();
        double ratio_e = (p_rec_e - p_gen_e)/p_gen_e;


        if (angle_neutron_missnucl*180./3.14159 < 5 && Q2 > 1 && Missing->M2() < 0.5 && Missing->M2() > -0.5 && Electron->P()>2 && Neutron->Theta()*180./3.141592 > 4 && Neutron->P()>0.25 && status_Kp >= 2000 && status_Kp <= 2999 && status_Km >= 2000 && status_Km <= 2999) {

                prof->Fill(p_rec, 100*ratio, real_weight_SDME_effcorrection); //tprofil stock des valeur et calcul automatiquement la moyenne (ponderé avec les poid) de ratio par bin de p_rec
                prof_Km->Fill(p_rec_Km, 100*ratio_Km, real_weight_SDME_effcorrection);
                prof_e->Fill(p_rec_e, 100*ratio_e, real_weight_SDME_effcorrection);
                h2_ratio_vs_prec->Fill(p_rec, 100*ratio, real_weight_SDME_effcorrection);


        }    

    }

    int bin1 = h2_ratio_vs_prec->GetXaxis()->FindBin(2.0);
    int bin2 = h2_ratio_vs_prec->GetXaxis()->FindBin(2.5);

    TH1D *h_ratio_2_25 = h2_ratio_vs_prec->ProjectionY(
     "h_ratio_2_25",
     bin1,
     bin2
    );

    h_ratio_2_25->SetTitle("Distribution de ratio pour 2 < p_{rec} < 2.5;ratio;Counts");
  


    // Canvas
    TCanvas *c1 = new TCanvas("c1", "Momentum correction", 800, 600);
    // Ouvrir le PDF (début)
c1->SaveAs("momentum_correction_kaons.pdf(");

// ===== PAGE 1 : Kaons =====
prof->SetLineColor(kRed);
prof->SetMarkerColor(kRed);
prof->SetMarkerStyle(24);
prof->SetMinimum(-3);

prof_Km->SetLineColor(kGreen);
prof_Km->SetMarkerColor(kGreen);
prof_Km->SetMarkerStyle(24);

prof->Draw("E1");
prof_Km->Draw("E1 SAME");

// Légende
TLegend *leg1 = new TLegend(0.6, 0.75, 0.88, 0.88);
leg1->AddEntry(prof, "K+", "lep");
leg1->AddEntry(prof_Km, "K-", "lep");
leg1->Draw();

// Sauvegarde page 1
c1->SaveAs("momentum_correction_kaons.pdf");

// ===== PAGE 2 : Electron =====
c1->Clear();  // IMPORTANT : nouvelle page propre

prof_e->SetLineColor(kBlue);
prof_e->SetMarkerColor(kBlue);
prof_e->SetMarkerStyle(24);

prof_e->Draw("E1");

// Légende électron
TLegend *leg2 = new TLegend(0.6, 0.75, 0.88, 0.88);
leg2->AddEntry(prof_e, "e^{-}", "lep");
leg2->Draw();
c1->SaveAs("momentum_correction_kaons.pdf");


// Fermer le PDF

c1->Clear(); h2_ratio_vs_prec->Draw("HIST"); c1->SaveAs("momentum_correction_kaons.pdf");

c1->Clear(); h_ratio_2_25->Draw("HIST");



c1->SaveAs("momentum_correction_kaons.pdf)");




}

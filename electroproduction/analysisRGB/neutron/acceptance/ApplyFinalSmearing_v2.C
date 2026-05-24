
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
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "TRandom3.h"
#include "TDatabasePDG.h"
#include "TVector3.h"

using namespace RooFit;

double Sigma(double x, double p1, double p2, double p3, double p4, double p5,
         double A, double B, double C, double D)
{
    //if (x < p1 || x > p5){
        //cout << x << endl;
    //}

    if (x < p2) return A;
    if (x < p3) return B;
    if (x < p4) return C;
    return D;
}


void ApplyFinalSmearing_v2()
{

    TString filename_tree = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_data2/fall2019_outbending_neutronKpKm_mc_v7.root";
    TFile *f_tree = TFile::Open(filename_tree, "UPDATE");
    TTree *t = (TTree*)f_tree->Get("tree");


    TRandom3 rng(0);

    TLorentzVector *Kp = nullptr;
    TLorentzVector *Km = nullptr;
    TLorentzVector Kp_smear;
    TLorentzVector Km_smear;
    TLorentzVector Missing_nucleon_smeared;
    TLorentzVector *Electron = nullptr;


    auto db=TDatabasePDG::Instance();
    double Eb = 10.40;
    double m_target = db->GetParticle(2112)->Mass(); 

    TLorentzVector beam(0,0, Eb, Eb);
    TLorentzVector target(0,0,0,m_target);

    TVector3 Kp_smear_vect;
    TVector3 Km_smear_vect;

    double t_missing_smeared;

    double Pkp_smear, Pkm_smear, sigmaSmearKp, sigmaSmearKm;

    double m_kaons = 0.4937;
    
    double A1 = 0.033;
    double B1 = 0.102;
    double C1 = 0.086;
    double D1 = 0.126;
    
    double A2 = 0.031;
    double B2 = 0.083;
    double C2 = 0.033;
    double D2 = 0.085;


    t->SetBranchAddress("Kp", &Kp);
    t->SetBranchAddress("Km", &Km);
    t->SetBranchAddress("Electron", &Electron);

    TBranch *bKp_smear = t->Branch("Kp_smear", "TLorentzVector", &Kp_smear);
    TBranch *bKm_smear = t->Branch("Km_smear", "TLorentzVector", &Km_smear);
    TBranch *bMissing_nucleon_smeared = t->Branch("Missing_nucleon_smeared", "TLorentzVector", &Missing_nucleon_smeared);
    TBranch *bt_missing_smeared = t->Branch("t_missing_smeared", &t_missing_smeared, "t_missing_smeared/D");


    Long64_t nentries = t->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {

        t->GetEntry(i);

        sigmaSmearKp = Sigma(Kp->P(), 0, 1.2, 1.6, 2.2, 6.0, A1, B1, C1, D1);
        sigmaSmearKm = Sigma(Km->P(), 0, 1.2, 1.6, 2.2, 6.0, A2, B2, C2, D2);

        Pkp_smear = Kp->P() * rng.Gaus(1.0, sigmaSmearKp);
        Pkm_smear = Km->P() * rng.Gaus(1.0, sigmaSmearKm);
                        
        Kp_smear_vect.SetMagThetaPhi(Pkp_smear, Kp->Theta(), Kp->Phi());
        Km_smear_vect.SetMagThetaPhi(Pkm_smear, Km->Theta(), Km->Phi());

        Kp_smear.SetVectM(Kp_smear_vect, m_kaons);
        Km_smear.SetVectM(Km_smear_vect, m_kaons); 

        Missing_nucleon_smeared =  beam + target - *Electron - Kp_smear - Km_smear;

        t_missing_smeared = (target - Missing_nucleon_smeared).M2();

        bKp_smear->Fill();
        bKm_smear->Fill();
        bMissing_nucleon_smeared->Fill();
        bt_missing_smeared->Fill();

    }

    t->Write("", TObject::kOverwrite);
    f_tree->Close();


    std::cout << "Branch smearing Kp and Km ajoutée avec succès !" << std::endl;


}
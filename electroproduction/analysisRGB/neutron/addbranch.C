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

void addbranch()
{
    TString filename = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_data2/fall2019_outbending_neutronKpKm_data_v2.root";
    TFile *f = TFile::Open(filename, "UPDATE");
    TTree *t = (TTree*)f->Get("tree");

    double Eb = 10.410; //rbg fall2019 outbending
    double m_target = 0.9395654;

    double cos_theta_H;
 

    TLorentzVector beam(0,0, Eb, Eb);
    TLorentzVector target(0,0,0,m_target); //cible neutron

    TLorentzVector *Neutron  = nullptr;
    TLorentzVector *Electron = nullptr;
    TLorentzVector *Kp = nullptr;
    TLorentzVector *Km = nullptr;
    TLorentzVector *Missing = nullptr;
    TLorentzVector *Missing_nucleon = nullptr;

    t->SetBranchAddress("Neutron",  &Neutron);
    t->SetBranchAddress("Electron", &Electron);
    t->SetBranchAddress("Kp", &Kp);
    t->SetBranchAddress("Km", &Km);
    t->SetBranchAddress("Missing", &Missing);
    t->SetBranchAddress("Missing_nucleon", &Missing_nucleon);

    TLorentzVector Phi;
    TBranch *bcos_theta_H = t->Branch("cos_theta_H", &cos_theta_H);

    Long64_t n = t->GetEntries();
    for (Long64_t i = 0; i < n; i++) {
        t->GetEntry(i);

        Phi = *Kp + *Km;
        TVector3 boost_phi = -Phi.BoostVector();
         
        TLorentzVector neutron_calcul = *Neutron;
        TLorentzVector Kp_calcul = *Kp;

        neutron_calcul.Boost(boost_phi);
        Kp_calcul.Boost(boost_phi);

        TVector3 k_dir = Kp_calcul.Vect().Unit(); 

        TVector3 z_axis = -neutron_calcul.Vect().Unit();

        cos_theta_H = k_dir.Dot(z_axis);
        bcos_theta_H->Fill();
    }

    t->Write("", TObject::kOverwrite);
    f->Close();
}

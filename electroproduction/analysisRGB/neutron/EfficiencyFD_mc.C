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


void EfficiencyFD_mc() {

    
    gROOT->Reset();
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
        data_adress = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_data2/fall2019_outbending_neutronKpKm_mc_efficiency.root";
    }


    TFile *f = TFile::Open(data_adress);
    TTree *tree = (TTree*)f->Get("tree;1");
    TTree *tree_gen = (TTree*)f->Get("tree_gen;1"); // ligne a commenter en finction de si c'est MC ou pas pour prendre les GEN avec les REC-MC


    //VARIABLE DATAs && MONTECARLO
    double Q2, t, W;
    double MinvKpKm, MM;
    double status_Kp, status_Km, status_n, status_el;
    double e_vx, e_vy, e_vz;
    double Kp_vx, Kp_vy, Kp_vz;
    double Km_vx, Km_vy, Km_vz;
    double n_vx, n_vy, n_vz;
    double angle_neutron_missnucl;
    double angle_neutron_neutrongen;
    double Lv_PCAL_n, Lu_PCAL_n, Lw_PCAL_n, x_PCAL_n, Lv_ECIN_n, Lu_ECIN_n, Lw_ECIN_n, x_ECIN_n, Lv_ECOUT_n, Lw_ECOUT_n, Lu_ECOUT_n, x_ECOUT_n;


    TLorentzVector *Electron = nullptr;
    TLorentzVector *Neutron = nullptr;
    TLorentzVector *Neutron_gen_associated = nullptr;
    TLorentzVector *Kp = nullptr;
    TLorentzVector *Km = nullptr;
    TLorentzVector *Missing = nullptr;

    tree->SetBranchAddress("Electron", &Electron);
    tree->SetBranchAddress("Neutron", &Neutron);
    tree->SetBranchAddress("Kp", &Kp);
    tree->SetBranchAddress("Km", &Km);
    tree->SetBranchAddress("Missing", &Missing);
    tree->SetBranchAddress("Neutron_gen_associated", &Neutron_gen_associated);
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("t", &t);
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
    tree->SetBranchAddress("angle_neutron_neutrongen", &angle_neutron_neutrongen);

    tree->SetBranchAddress("Lv_PCAL_n", &Lv_PCAL_n);
    tree->SetBranchAddress("Lu_PCAL_n", &Lu_PCAL_n);
    tree->SetBranchAddress("Lw_PCAL_n", &Lw_PCAL_n);
    tree->SetBranchAddress("x_PCAL_n", &x_PCAL_n);

    tree->SetBranchAddress("Lv_ECIN_n", &Lv_ECIN_n);
    tree->SetBranchAddress("Lu_ECIN_n", &Lu_ECIN_n);
    tree->SetBranchAddress("Lw_ECIN_n", &Lw_ECIN_n);
    tree->SetBranchAddress("x_ECIN_n", &x_ECIN_n);

    tree->SetBranchAddress("Lv_ECOUT_n", &Lv_ECOUT_n);
    tree->SetBranchAddress("Lu_ECOUT_n", &Lu_ECOUT_n);
    tree->SetBranchAddress("Lw_ECOUT_n", &Lw_ECOUT_n);
    tree->SetBranchAddress("x_ECOUT_n", &x_ECOUT_n);




    //VARIABLE SPECIFIQUE MC

    double weight;
    double real_weight;

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

    double is_e_rec;


    if(isMC == true){

        tree->SetBranchAddress("weight", &weight);
        tree->SetBranchAddress("real_weight", &real_weight);

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

        tree_gen->SetBranchAddress("is_e_rec", &is_e_rec);

    }

    
    TH2F *thetavsphi_neutron_rec = new TH2F("thetavsphi_neutron_rec", "theta vs phi neutron REC", 400, -200.0, 250, 400, -10.0, 45);
    TH2F *thetavsphi_neutron_gen = new TH2F("thetavsphi_neutron_gen", "theta vs phi neutron GEN", 400, -200.0, 250, 400, -10.0, 45);

    auto setTitles = [](TH1 *h, const char *xtitle, const char *ytitle) {
               h->GetXaxis()->SetTitle(xtitle);
               h->GetYaxis()->SetTitle(ytitle);
               h->GetXaxis()->CenterTitle();
               h->GetYaxis()->CenterTitle();
   };

   setTitles(thetavsphi_neutron_rec, "#phi[degree]", "#theta [degree]");
   setTitles(thetavsphi_neutron_gen, "#phi[degree]", "#theta [degree]");


   const int nPbins = 30;
   double pmin = 0.3;
   double pmax = 3.5;


   TH1F *h_rec_p = new TH1F("h_rec_p", "p neutron REC", nPbins, pmin, pmax);
   TH1F *h_gen_p = new TH1F("h_gen_p", "p neutron GEN", nPbins, pmin, pmax);
   TH1F *h_eff_p = new TH1F("h_eff_p", "Efficiency vs p", nPbins, pmin, pmax);

   TH1F *h_rec_p_2 = new TH1F("h_rec_p_2", "p neutron REC", nPbins, pmin, pmax);
   TH1F *h_gen_p_2 = new TH1F("h_gen_p_2", "p neutron GEN", nPbins, pmin, pmax);
   TH1F *h_eff_p_2 = new TH1F("h_eff_p_2", "Efficiency vs p", nPbins, pmin, pmax);

   TH1F *h_rec_p_tot = new TH1F("h_rec_p_tot", "p neutron REC", nPbins, pmin, pmax);
   TH1F *h_gen_p_tot = new TH1F("h_gen_p_tot", "p neutron GEN", nPbins, pmin, pmax);
   TH1F *h_eff_p_tot = new TH1F("h_eff_p_tot", "Efficiency vs p", nPbins, pmin, pmax);

   double pt = 2.146;  // seuil de transition (à ajuster)

    TF1 *f_eff_data = new TF1("f_eff_data",  [pt](double *x, double *par) { double P = x[0];

        if (P < pt) {
            return par[0]
                 + par[1]*P
                 + par[2]*P*P
                 + par[3]*P*P*P;
        } else {
            return par[4] * (1.0 - 1.0/(1.0 + exp((P - par[5])/par[6])));
        }

    },
    0.25, 3.5, 7);

    f_eff_data->SetParameters(
    -0.1817,   // a0
    0.6187,   // a1
   -0.0605,  // a2
    -0.0179,  // a3
    0.7884,   // a4
    0.0086,   // a5
    1.0796    // a6
    );

   h_rec_p->Sumw2();
   h_gen_p->Sumw2();
   h_eff_p->Sumw2();

   h_rec_p_2->Sumw2();
   h_gen_p_2->Sumw2();
   h_eff_p_2->Sumw2();

   h_rec_p_tot->Sumw2();
   h_gen_p_tot->Sumw2();
   h_eff_p_tot->Sumw2();

   h_rec_p->GetXaxis()->SetTitle("p_{neutron_gen} [GeV]");
   h_rec_p->GetYaxis()->SetTitle("Counts");

   h_rec_p_2->GetXaxis()->SetTitle("p_{neutron_gen} [GeV]");
   h_rec_p_2->GetYaxis()->SetTitle("Counts");

   h_rec_p_tot->GetXaxis()->SetTitle("p_{neutron_gen} [GeV]");
   h_rec_p_tot->GetYaxis()->SetTitle("Counts");

   h_gen_p->GetXaxis()->SetTitle("p_{neutron_gen} [GeV]");
   h_gen_p->GetYaxis()->SetTitle("Counts");

   h_gen_p_2->GetXaxis()->SetTitle("p_{neutron_gen} [GeV]");
   h_gen_p_2->GetYaxis()->SetTitle("Counts");

   h_gen_p_tot->GetXaxis()->SetTitle("p_{neutron_gen} [GeV]");
   h_gen_p_tot->GetYaxis()->SetTitle("Counts");

   h_eff_p->GetXaxis()->SetTitle("p_{neutron_gen} [GeV]");
   h_eff_p->GetYaxis()->SetTitle("Efficiency");

   h_eff_p_2->GetXaxis()->SetTitle("p_{neutron_gen} [GeV]");
   h_eff_p_2->GetYaxis()->SetTitle("Efficiency");

   h_eff_p_tot->GetXaxis()->SetTitle("p_{neutron_gen} [GeV]");
   h_eff_p_tot->GetYaxis()->SetTitle("Efficiency");



   Long64_t nentries_bis = tree->GetEntries();
    for (Long64_t k = 0; k < nentries_bis; k++) {

        tree->GetEntry(k);

        //acceptance cut for REC

        if(status_n >= 2000 && status_n <= 2999 && ((Lv_PCAL_n > 9.0 && Lw_PCAL_n > 9.0) || x_PCAL_n == 0) && ((Lv_ECIN_n > 9.0 && Lw_ECIN_n > 9.0) || x_ECIN_n == 0) && ((Lv_ECOUT_n > 9.0 && Lw_ECOUT_n > 9.0) || x_ECOUT_n == 0) && (x_PCAL_n != 0 || x_ECIN_n != 0 || x_ECOUT_n != 0) && Neutron->Theta()*180./3.141592 > 4 &&  angle_neutron_neutrongen*180./3.14159 < 11.48 ){

            thetavsphi_neutron_rec->Fill(Neutron->Phi()*180./3.141592, Neutron->Theta()*180./3.141592, weight);

            if(Neutron_gen_associated->P()< pmax && Neutron_gen_associated->P() > pmin){

                h_rec_p_tot->Fill(Neutron_gen_associated->P(), weight);

                if(Neutron_gen_associated->Theta()*180./3.141592 >= 4 && Neutron_gen_associated->Theta()*180./3.141592 <= 20 ){

                    h_rec_p->Fill(Neutron_gen_associated->P(), weight);

                }

                if(Neutron_gen_associated->Theta()*180./3.141592 > 20 && Neutron_gen_associated->Theta()*180./3.141592 <= 35 ){

                    h_rec_p_2->Fill(Neutron_gen_associated->P(), weight);

                }

            }

        }

    }


    nentries_bis = tree_gen->GetEntries();
    for (Long64_t k = 0; k < nentries_bis; k++) {

        tree_gen->GetEntry(k);

        //acceptance cut for GEN with histo rec
        int bin = thetavsphi_neutron_rec->FindBin(Neutron_gen->Phi()*180./3.141592, Neutron_gen->Theta()*180./3.141592);

        if (thetavsphi_neutron_rec->GetBinContent(bin) > 0 && Neutron_gen->Theta()*180./3.141592 < 40 && Neutron_gen->P()>0.25 && is_e_rec == 1) {

            thetavsphi_neutron_gen->Fill(Neutron_gen->Phi()*180./3.141592,  Neutron_gen->Theta()*180./3.141592, weight_gen);

            if(Neutron_gen->P() < pmax && Neutron_gen->P() > pmin){

                h_gen_p_tot->Fill(Neutron_gen->P(), weight_gen);

                if(Neutron_gen->Theta()*180./3.141592 >= 4 && Neutron_gen->Theta()*180./3.141592 <= 20 ){

                     h_gen_p->Fill(Neutron_gen->P(), weight_gen);

                }

                if(Neutron_gen->Theta()*180./3.141592 >= 20 && Neutron_gen->Theta()*180./3.141592 <= 35 ){

                     h_gen_p_2->Fill(Neutron_gen->P(), weight_gen);

                }

            }

        }

    }

   h_eff_p->Divide(h_rec_p, h_gen_p);
   h_eff_p_2->Divide(h_rec_p_2, h_gen_p_2);
   h_eff_p_tot->Divide(h_rec_p_tot, h_gen_p_tot);

   h_eff_p->SetMinimum(0.0);
   h_eff_p->SetMaximum(1.1);

   h_eff_p_2->SetMinimum(0.0);
   h_eff_p_2->SetMaximum(1.1);

   h_eff_p_tot->SetMinimum(0.0);
   h_eff_p_tot->SetMaximum(1.1);

   h_rec_p->SetMaximum(4*h_rec_p->GetMaximum());
   h_gen_p->SetMaximum(4*h_gen_p->GetMaximum());


   TH1F *h_ratio_p = (TH1F*)h_eff_p_tot->Clone("h_ratio_p");
   h_ratio_p->Reset();
   h_ratio_p->SetTitle("Data/MC Efficiency Ratio");
   h_ratio_p->GetYaxis()->SetTitle("Data / MC");

   for (int i = 1; i <= h_eff_p_tot->GetNbinsX(); i++) {

    double p_center = h_eff_p_tot->GetBinCenter(i);

    double eff_mc = h_eff_p_tot->GetBinContent(i);
    double err_mc = h_eff_p_tot->GetBinError(i);

    double eff_data = f_eff_data->Eval(p_center);

    if (eff_data > -1) {

        double ratio = eff_data / eff_mc;

        // propagation d'erreur (data supposée sans erreur ici)
        double err_ratio = 0.01;

        h_ratio_p->SetBinContent(i, ratio);
        h_ratio_p->SetBinError(i, err_ratio);

    } 
    
    }

    h_ratio_p->SetMinimum(0.0);
    h_ratio_p->SetMaximum(1.5);

    h_ratio_p->SetMarkerStyle(20);
    h_ratio_p->SetMarkerColor(kGreen);
    h_ratio_p->SetLineColor(kGreen);


    //pdf


   TString pdfFile = "EfficiencyFD_mc_plots.pdf";
   TCanvas *c = new TCanvas("c", "Plots", 800, 600);

   thetavsphi_neutron_rec->Draw("COLZ"); c->Print(pdfFile + "("); 
   thetavsphi_neutron_gen->Draw("COLZ"); c->Print(pdfFile); 

   h_rec_p->SetMarkerStyle(20);
   h_rec_p->SetMarkerColor(kBlue);
   h_rec_p->SetLineColor(kBlue);

   h_rec_p_2->SetMarkerStyle(24);
   h_rec_p_2->SetMarkerColor(kRed);
   h_rec_p_2->SetLineColor(kRed);

   h_rec_p_tot->SetMarkerStyle(21);
   h_rec_p_tot->SetMarkerColor(kBlack);
   h_rec_p_tot->SetLineColor(kBlack);

   h_rec_p->Draw("E1");  
   h_rec_p_2->Draw("E1 SAME");
   h_rec_p_tot->Draw("E1 SAME");

   TLegend *leg1 = new TLegend(0.6,0.7,0.85,0.85);
   leg1->AddEntry(h_rec_p,"4 < #theta < 20","lep");
   leg1->AddEntry(h_rec_p_2,"20 < #theta < 35","lep");
   leg1->AddEntry(h_rec_p_tot,"all #theta","lep");
   leg1->Draw();
   
   
   c->Print(pdfFile); 

   h_gen_p->SetMarkerStyle(20);
   h_gen_p->SetMarkerColor(kBlue);
   h_gen_p->SetLineColor(kBlue);

   h_gen_p_2->SetMarkerStyle(24);
   h_gen_p_2->SetMarkerColor(kRed);
   h_gen_p_2->SetLineColor(kRed);

   h_gen_p_tot->SetMarkerStyle(21);
   h_gen_p_tot->SetMarkerColor(kBlack);
   h_gen_p_tot->SetLineColor(kBlack);

   h_gen_p->Draw("E1");
   h_gen_p_2->Draw("E1 SAME");
   h_gen_p_tot->Draw("E1 SAME");

   TLegend *leg2 = new TLegend(0.6,0.7,0.85,0.85);
   leg2->AddEntry(h_gen_p,"4 < #theta < 20","lep");
   leg2->AddEntry(h_gen_p_2,"20 < #theta < 35","lep");
   leg2->AddEntry(h_gen_p_tot,"all #theta","lep");
   leg2->Draw();
   
   
   c->Print(pdfFile); 

   h_eff_p->SetMarkerStyle(20);
   h_eff_p->SetMarkerColor(kBlue);
   h_eff_p->SetLineColor(kBlue);

   h_eff_p_2->SetMarkerStyle(24);
   h_eff_p_2->SetMarkerColor(kRed);
   h_eff_p_2->SetLineColor(kRed);

   h_eff_p_tot->SetMarkerStyle(21);
   h_eff_p_tot->SetMarkerColor(kBlack);
   h_eff_p_tot->SetLineColor(kBlack);

   h_eff_p->Draw("E1");
   h_eff_p_2->Draw("E1 SAME");
   h_eff_p_tot->Draw("E1 SAME");

   TLegend *leg3 = new TLegend(0.6,0.7,0.85,0.85);
   leg3->AddEntry(h_eff_p,"4 < #theta < 20","lep");
   leg3->AddEntry(h_eff_p_2,"20 < #theta < 35","lep");
   leg3->AddEntry(h_eff_p_tot,"all #theta","lep");
   leg3->Draw();

   c->Print(pdfFile); 

   f_eff_data->SetRange(0.01, 7.5);
   f_eff_data->Draw();
   c->Print(pdfFile);

   h_ratio_p->Draw("E1");

   c->Print(pdfFile + ")");

   delete c;

   TFile *fout = new TFile("ratio_efficiency.root", "RECREATE");
   h_ratio_p->Write();
   fout->Close();
   delete fout;



}

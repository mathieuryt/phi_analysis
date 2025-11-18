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


void postrootMC_DisplacedVertex() {
    gROOT->Reset();
    gSystem->Load("libPhysics.so");
    gSystem->Load("libCore.so");
    gSystem->Load("libTree.so");
    gSystem->Load("libRIO.so");

    gSystem->Load("libHist.so");

    gStyle->SetPalette(kBird); 
    // ouverture fichier root
    TFile *f = TFile::Open("/Users/mr282803/Documents/generateur/inbendingMC_DisplacedVertex.root");
    // recuperer t tree
    TTree *tree = (TTree*)f->Get("tree;1");

    // Variables pour lire les branches
    double W, Q2, t;
    double weight;
    double Sweight = 0.0;
    double MM, Minv_pip_pim, Minv_Ks_Kl;
    double status_pim, status_el, status_pip, status_pr;

    double e_vx, e_vy, e_vz;
    double pip_vx, pip_vy, pip_vz;
    double pim_vx, pim_vy, pim_vz;
    double Ks_vx, Ks_vy, Ks_vz;


    TLorentzVector *el = nullptr;
    TLorentzVector *pr = nullptr;
    TLorentzVector *pip = nullptr;
    TLorentzVector *pim = nullptr;

    double bin = 150;
    double bin2 = 50;
    double bin3 = 80;
    double bin4 = 200;
    double bin5 = 70;
    double bin6 = 180;
    double echelle_plot_couleur = 3000;

    // Lier les branches
    tree->SetBranchAddress("e_vx", &e_vx);
    tree->SetBranchAddress("e_vy", &e_vy);
    tree->SetBranchAddress("e_vz", &e_vz);

    tree->SetBranchAddress("Ks_vx", &Ks_vx);
    tree->SetBranchAddress("Ks_vy", &Ks_vy);
    tree->SetBranchAddress("Ks_vz", &Ks_vz);
    
    tree->SetBranchAddress("pip_vx", &pip_vx);
    tree->SetBranchAddress("pip_vy", &pip_vy);
    tree->SetBranchAddress("pip_vz", &pip_vz);

    tree->SetBranchAddress("pim_vx", &pim_vx);
    tree->SetBranchAddress("pim_vy", &pim_vy);
    tree->SetBranchAddress("pim_vz", &pim_vz);

    tree->SetBranchAddress("W", &W);
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("t", &t);

    tree->SetBranchAddress("Electron", &el);
    tree->SetBranchAddress("Proton", &pr);

    tree->SetBranchAddress("PiPlus", &pip);
    tree->SetBranchAddress("PiMinus", &pim);

    tree->SetBranchAddress("MinvKsKl", &Minv_Ks_Kl);
    tree->SetBranchAddress("MissingMass", &MM);
    tree->SetBranchAddress("MinvPiPlusMoins", &Minv_pip_pim);

    tree->SetBranchAddress("Status_pim", &status_pim);
    tree->SetBranchAddress("Status_pip", &status_pip);
    tree->SetBranchAddress("Status_pr", &status_pr);
    tree->SetBranchAddress("Status_el", &status_el);

    tree->SetBranchAddress("weight", &weight);

    cout << "Variables initialisees" <<endl;

    // Créer l’histogramme 
    TH1F *vz_diff1 = new TH1F("vz_diff1", "Vz_{e-} - Vz_{#pi+}", bin, -15, 15);
    TH1F *vz_diff2 = new TH1F("vz_diff2", "Vz_{#pi-} - Vz_{#pi+}", bin, -15, 15);
    TH1F *vz_diff3 = new TH1F("vz_diff3", "Vz_{e-} - Vz_{#pi-}", bin, -15, 15);
    TH1F *vz_diff4 = new TH1F("vz_diff4", "Vz_{e-} - Vz_{Ks}", bin, -15, 15);

    TH1F *Dist1 = new TH1F("Dist1", "Distance vertex e- and Ks", bin, -1, 20);
    TH1F *Dist2 = new TH1F("Dist2", "Distance vertex #pi+ and #pi-", bin, -1, 20);

    TH1F *t_hist = new TH1F("t_hist", "t", bin2, -10, 0);
    TH1F *W_hist = new TH1F("W_hist", "W", bin2, 0, 5);
    TH1F *Q2_hist = new TH1F("Q2_hist", "Q^{2}", bin2, 0, 10);

    TH1F *weight_hist = new TH1F("weight_hist", "weight", 70, 0, 0.001);

    TH1F *MinvPiPlusMoins_hist = new TH1F("MinvPiPlusMoins_hist", "Invariant Mass #pi^{-} #pi^{+}", 100, 0, 1.5);
    TH1F *MissingMass_hist = new TH1F("MissingMass_hist", "Missing Mass", 100, 0, 1.5);
    TH1F *MinvKsKl_hist = new TH1F("MinvKsKl_hist", "Invariant Mass K_{s} K_{L}", 100, 0.7, 1.5);

    TH2F *MinvKsKl_W = new TH2F("MinvKsKl_W", "Invariant Mass Ks Kl vs W", bin3, 0, 7, bin3, 0, 5.5); 
    TH2F *MinvKsKl_Q2 = new TH2F("MinvKsKl_Q2", "Invariant Mass Ks Kl vs Q^{2}", bin3, 0, 12, bin3, 0, 5.5); 
    TH2F *MinvKsKl_Thetae = new TH2F("MinvKsKl_Thetae", "Invariant Mass Ks Kl vs Theta electron", bin3, 0, 50, bin3, 0, 3.5); 

    TH1F *MinvPiPlusMoins_hist2 = new TH1F("MinvPiPlusMoins_hist2", "Invariant Mass #pi^{-} #pi^{+} with cut on MM", 100, 0, 1.5);
    TH1F *MissingMass_hist2 = new TH1F("MissingMass_hist2", "Missing Mass with cut on M_{inv} of #pi^{+}#pi^{-}", 100, 0, 1.5);

    TH2F *theta_p_el = new TH2F("theta_p_el", "Theta vs p for electron", bin5, 0, 100, bin5, 0, 13);
    TH2F *theta_p_pr = new TH2F("theta_p_pr", "Theta vs p for proton", bin5, 0, 60, bin5, 0, 9); 
    TH2F *theta_p_pip = new TH2F("theta_p_pip", "Theta vs p for #pi^{+}", bin5, 0, 140, bin5, 0, 7); 
    TH2F *theta_p_pim = new TH2F("theta_p_pim", "Theta vs p for #pi^{-}", bin5, 0, 140, bin5, 0, 7);

    TH2F *theta_p_el2 = new TH2F("theta_p_el2", "Theta vs p for electron with all cuts", bin5, 0, 100, bin5, 0, 13);
    TH2F *theta_p_pr2 = new TH2F("theta_p_pr2", "Theta vs p for proton with all cuts", bin5, 0, 60, bin5, 0, 9); 
    TH2F *theta_p_pip2 = new TH2F("theta_p_pip2", "Theta vs p for #pi^{+} with all cuts", bin5, 0, 140, bin5, 0, 7); 
    TH2F *theta_p_pim2 = new TH2F("theta_p_pim2", "Theta vs p for #pi^{-} with all cuts", bin5, 0, 140, bin5, 0, 7); 

    //TH2F *MKs_R_Ks_e = new TH2F("MKs_R_Ks_e", "M_{Ks} vs Dist(Ks-electron)", bin6, 0, 10, bin6, 0, 2); 
    //TH2F *MKl_R_Ks_e = new TH2F("MKl_R_Ks_e", "M_{Kl} vs Dist(Ks-electron)", bin6, 0, 10, bin6, 0, 2.5); 
    //TH2F *MPhi_R_Ks_e = new TH2F("MPhi_R_Ks_e", "M_{#phi} vs Dist(Ks-electron)", bin6, 0, 10, bin6, 0, 3.5);

    //TH2F *MKs_R_pip_pim = new TH2F("MKs_R_pip_pim", "M_{Ks} vs Dist(#pi^{+} #pi^{-})", bin6, 0, 6, bin6, 0, 2); 
    //TH2F *MKl_R_pip_pim = new TH2F("MKl_R_pip_pim", "M_{Kl} vs Dist(#pi^{+} #pi^{-})", bin6, 0, 6, bin6, 0, 2.5); 
    //TH2F *MPhi_R_pip_pim = new TH2F("MPhi_R_pip_pim", "M_{#phi} vs Dist(#pi^{+} #pi^{-})", bin6, 0, 6, bin6, 0, 3.5);

    //TH2F *MKs_R_Ks_e2 = new TH2F("MKs_R_Ks_e2", "M_{Ks} vs Dist(Ks-electron) with cut on MM", bin6, 0, 6, bin6, 0, 2); 
    //TH2F *MKl_R_Ks_e2 = new TH2F("MKl_R_Ks_e2", "M_{Kl} vs Dist(Ks-electron) with cut on M_{inv} #pi^{+} #pi^{-}", bin6, 0, 6, bin6, 0, 2.5); 
    //TH2F *MPhi_R_Ks_e2 = new TH2F("MPhi_R_Ks_e2", "M_{#phi} vs Dist(Ks-electron) with cut on MM && M_{inv} #pi^{+} #pi^{-}", bin6, 0, 6, bin6, 0, 3.5);

    //TH2F *MKs_R_pip_pim2 = new TH2F("MKs_R_pip_pim2", "M_{Ks} vs Dist(#pi^{+} #pi^{-}) with cut on MM", bin6, 0, 3, bin6, 0, 2); 
    //TH2F *MKl_R_pip_pim2 = new TH2F("MKl_R_pip_pim2", "M_{Kl} vs Dist(#pi^{+} #pi^{-}) with cut on M_{inv} #pi^{+} #pi^{-}", bin6, 0, 3, bin6, 0, 2.5); 
    //TH2F *MPhi_R_pip_pim2 = new TH2F("MPhi_R_pip_pim2", "M_{#phi} vs Dist(#pi^{+} #pi^{-}) with cut on MM && M_{inv} #pi^{+} #pi^{-}", bin6, 0, 3, bin6, 0, 3.5); 

    

    auto setTitles = [](TH1 *h, const char *xtitle, const char *ytitle) {
        h->GetXaxis()->SetTitle(xtitle);
        h->GetYaxis()->SetTitle(ytitle);
        h->GetXaxis()->CenterTitle();
        h->GetYaxis()->CenterTitle();
    };

    setTitles(vz_diff1, "#Delta v_{z} (e - #pi^{+}) [cm]", "");
    setTitles(vz_diff2, "#Delta v_{z} (#pi^{-} - #pi^{+}) [cm]", "");
    setTitles(vz_diff3, "#Delta v_{z} (e - #pi^{-}) [cm]", "");
    setTitles(vz_diff4, "#Delta v_{z} (e - Ks) [cm]", "");

    setTitles(Dist1, "#sqrt{(v_{Ksx} - v_{ex})^{2} + (v_{Ksy} - v_{ey})^{2} + (v_{Ksz} - v_{ez})^{2}} [cm]", "");
    setTitles(Dist2, "#sqrt{(v_{#pi^{+}_{x}} - v_{#pi^{-}_{x}})^{2} + (v_{#pi^{+}_{y}} - v_{#pi^{-}_{y}})^{2} + (v_{#pi^{+}_{z}} - v_{#pi^{-}_{z}})^{2}} [cm]", "");

    setTitles(MinvKsKl_W, "W [GeV]", "M_{inv}(K_{s}K_{L}) [GeV]");
    setTitles(MinvKsKl_Q2, "Q^{2} [GeV^{2}]", "M_{inv}(K_{s}K_{L}) [GeV]");
    setTitles(MinvKsKl_Thetae, "#theta_{e} [degree]", "M_{inv}(K_{s}K_{L}) [GeV]");

    setTitles(t_hist, "t [GeV^{2}]", "");
    setTitles(Q2_hist, "Q^{2} [GeV^{2}]", "");
    setTitles(W_hist, "W [GeV]", "");

    setTitles(weight_hist, "weight of events", "");

    setTitles(MinvPiPlusMoins_hist, "M_{inv}(#pi^{-} #pi^{+})", "");
    setTitles(MissingMass_hist, "MissingMass [GeV]", "");
    setTitles(MinvKsKl_hist, "M_{inv}(K_{s}K_{L}) [GeV]", "");

    setTitles(MinvPiPlusMoins_hist2, "M_{inv}(#pi^{-} #pi^{+})", "");
    setTitles(MissingMass_hist2, "MissingMass [GeV]", "");

    setTitles(theta_p_el, "#theta [degree]", "p [GeV]");
    setTitles(theta_p_pr, "#theta [degree]", "p [GeV]");
    setTitles(theta_p_pip, "#theta [degree]", "p [GeV]");
    setTitles(theta_p_pim, "#theta [degree]", "p [GeV]");

    setTitles(theta_p_el2, "#theta [degree]", "p [GeV]");
    setTitles(theta_p_pr2, "#theta [degree]", "p [GeV]");
    setTitles(theta_p_pip2, "#theta [degree]", "p [GeV]");
    setTitles(theta_p_pim2, "#theta [degree]", "p [GeV]");

    //setTitles(MKs_R_Ks_e,"Dist(Ks-electron) [cm]" , "M_{Ks} [GeV]");
    //setTitles(MKl_R_Ks_e,"Dist(Ks-electron) [cm]", "M_{Kl} [GeV]");
    //setTitles(MPhi_R_Ks_e,"Dist(Ks-electron) [cm]", "M_{#phi} [GeV]");

    //setTitles(MKs_R_pip_pim,"Dist(#pi^{+} #pi^{-}) [cm]", "M_{Ks} [GeV]");
    //setTitles(MKl_R_pip_pim,"Dist(#pi^{+} #pi^{-}) [cm]", "M_{Kl} [GeV]");
    //setTitles(MPhi_R_pip_pim,"Dist(#pi^{+} #pi^{-}) [cm]", "M_{Phi} [GeV]");

    //setTitles(MKs_R_Ks_e2,"Dist(Ks-electron) [cm]" , "M_{Ks} [GeV]");
    //setTitles(MKl_R_Ks_e2,"Dist(Ks-electron) [cm]", "M_{Kl} [GeV]");
    //setTitles(MPhi_R_Ks_e2,"Dist(Ks-electron) [cm]", "M_{#phi} [GeV]");

    //setTitles(MKs_R_pip_pim2,"Dist(#pi^{+} #pi^{-}) [cm]", "M_{Ks} [GeV]");
    //setTitles(MKl_R_pip_pim2,"Dist(#pi^{+} #pi^{-}) [cm]", "M_{Kl} [GeV]");
    //setTitles(MPhi_R_pip_pim2,"Dist(#pi^{+} #pi^{-}) [cm]", "M_{Phi} [GeV]");

    cout << "Histogrammes initialises" <<endl;



    // Remplir l’histogramme
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {

        tree->GetEntry(i);

        if(Q2 > 1){

        double delta_vz1 = e_vz - pip_vz;
        double delta_vz2 = pim_vz - pip_vz;
        double delta_vz3 = e_vz - pim_vz;
        double delta_vz4 = e_vz - Ks_vz;

        vz_diff1->Fill(delta_vz1);
        vz_diff2->Fill(delta_vz2);
        vz_diff3->Fill(delta_vz3);
        vz_diff4->Fill(delta_vz4);

        double R_Ks_e = sqrt((Ks_vx - e_vx)*(Ks_vx - e_vx) + (Ks_vy - e_vy)*(Ks_vy - e_vy) + (Ks_vz - e_vz)*(Ks_vz - e_vz));
        double R_pip_pim = sqrt((pip_vx - pim_vx)*(pip_vx - pim_vx) + (pip_vy - pim_vy)*(pip_vy - pim_vy) + (pip_vz - pim_vz)*(pip_vz - pim_vz));

        Dist1->Fill(R_Ks_e);
        Dist2->Fill(R_pip_pim);

        //status_pim >= 2000 && status_pim <= 2999 && status_pip >= 2000 && status_pip <= 2999 && status_pr >= 2000 && status_pr <= 2999

        t_hist->Fill(t, weight);
        Q2_hist->Fill(Q2, weight);
        W_hist->Fill(W, weight);

        weight_hist->Fill(weight);

        MinvKsKl_W->Fill(W, Minv_Ks_Kl, weight);
        MinvKsKl_Q2->Fill(Q2, Minv_Ks_Kl, weight);
        MinvKsKl_Thetae->Fill(el->Theta()*(180/3.14), Minv_Ks_Kl, weight);

        theta_p_el->Fill(el->Theta()*(180/3.14), el->P(), weight);
        theta_p_pr->Fill(pr->Theta()*(180/3.14), pr->P(), weight);
        theta_p_pim->Fill(pim->Theta()*(180/3.14), pim->P(), weight);
        theta_p_pip->Fill(pip->Theta()*(180/3.14), pip->P(), weight);

        //MKs_R_Ks_e->Fill(R_Ks_e, Minv_pipD_pimD);
        //MKl_R_Ks_e->Fill(R_Ks_e, MMD);
        //MPhi_R_Ks_e->Fill(R_Ks_e, Minv_Ks_KlD);

        //MKs_R_pip_pim->Fill(R_pip_pim, Minv_pipD_pimD);
        //MKl_R_pip_pim->Fill(R_pip_pim, MMD);
        //MPhi_R_pip_pim->Fill(R_pip_pim, Minv_Ks_KlD);

        MinvPiPlusMoins_hist->Fill(Minv_pip_pim, weight);
        MissingMass_hist->Fill(MM, weight);


        if(MM < 0.7 && MM > 0.4){ // cut FD + cut missing mass

                    MinvPiPlusMoins_hist2->Fill(Minv_pip_pim, weight);
                    //MKs_R_Ks_e2->Fill(R_Ks_e, Minv_pipD_pimD);
                    //MKs_R_pip_pim2->Fill(R_pip_pim, Minv_pipD_pimD);

        }

        if(Minv_pip_pim < 0.6 && Minv_pip_pim > 0.3){ // cut FD + cut masse invariante pi+ pi-

                    MissingMass_hist2->Fill(MM, weight);
                    //MKl_R_Ks_e2->Fill(R_Ks_e, MMD);
                    //MKl_R_pip_pim2->Fill(R_pip_pim, MMD);
                    

        }

        if(MM < 0.7 && MM > 0.4 && Minv_pip_pim > 0.3 && Minv_pip_pim < 0.6){ // cut FD + cut masse invariante pi+ pi- + cut missing mass

                    MinvKsKl_hist->Fill(Minv_Ks_Kl, weight);
                    //MPhi_R_Ks_e2->Fill(R_Ks_e, Minv_Ks_KlD);
                    //MPhi_R_pip_pim2->Fill(R_pip_pim, Minv_Ks_KlD);

        }

        if(Minv_Ks_Kl < 1.1 && Minv_Ks_Kl > 0.95 && MM < 0.7 && MM > 0.4 && Minv_pip_pim > 0.3 && Minv_pip_pim < 0.6){
                        
                    theta_p_el2->Fill(el->Theta()*(180/3.14), el->P(), weight);
                    theta_p_pr2->Fill(pr->Theta()*(180/3.14), pr->P(), weight);
                    theta_p_pim2->Fill(pim->Theta()*(180/3.14), pim->P(), weight);
                    theta_p_pip2->Fill(pip->Theta()*(180/3.14), pip->P(), weight);

                    Sweight += weight;
        }


        //cout << delta_vz1 << endl;
        //cout << delta_vt1 << endl;

        }

        
    }

    cout << "Histogrammes remplis" <<endl;

    TString pdfFile = "analysisMC_inbending_DisplacedVertex.pdf";

    // Premier canvas (ouvre le PDF)
    TCanvas *c = new TCanvas("c", "Plots", 800, 600);

    //Q2_hist->Draw();

    c->Clear();
    c->SetFillColor(0);
    c->cd();

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextAlign(13);

    latex.DrawLatexNDC(0.1, 0.9, "File : MC inbending with vertex correction");
    latex.DrawLatexNDC(0.1, 0.8, "Number of generated events: 2 000 000");

    latex.DrawLatexNDC(0.1, 0.6, "Summary of cuts for the next plots:");
    latex.DrawLatexNDC(0.1, 0.5, "- Only 1 proton, 1 #pi^{+}, 1 #pi^{-}, and 1 e^{-}");
    latex.DrawLatexNDC(0.1, 0.4, "- Q^{2} > 1 -> electron FD");
    latex.DrawLatexNDC(0.1, 0.3, "- Very large cut on Missing mass, Invariant mass #pi^{+} #pi^{-}");
    latex.DrawLatexNDC(0.1, 0.2, "  and Invariant mass Ks Kl (cut between 0 and 3 GeV)");

    c->Print(pdfFile + "("); // ouvre le fichier PDF


    Q2_hist->Draw("HIST");
    c->Print(pdfFile);


    // Ajouter les autres plots
    t_hist->Draw("HIST"); c->Print(pdfFile);
    W_hist->Draw("HIST"); c->Print(pdfFile);

    vz_diff1->Draw(); c->Print(pdfFile);

    vz_diff3->Draw();
    c->Print(pdfFile);

    vz_diff4->Draw(); 
    c->Print(pdfFile);

    Dist1->Draw(); 
    c->Print(pdfFile);

    vz_diff2->Draw(); 
    TLine *line6 = new TLine(-5, 0, -5, vz_diff2->GetMaximum());
    line6->SetLineColor(kBlack);     // couleur rouge
    line6->SetLineStyle(2);        // 2 = pointillé
    line6->SetLineWidth(2);
    line6->Draw("same");
    TLine *line7 = new TLine(5, 0, 5, vz_diff2->GetMaximum());
    line7->SetLineColor(kBlack);     // couleur rouge
    line7->SetLineStyle(2);        // 2 = pointillé
    line7->SetLineWidth(2);
    line7->Draw("same");
    c->Print(pdfFile);

    Dist2->Draw(); 
    c->Print(pdfFile);

    weight_hist->Draw("HIST"); c->Print(pdfFile);

    MinvPiPlusMoins_hist->Draw("HIST"); 
    TLine *line = new TLine(0.496, 0, 0.496, MinvPiPlusMoins_hist->GetMaximum());
    line->SetLineColor(kRed);     // couleur rouge
    line->SetLineStyle(2);        // 2 = pointillé
    line->SetLineWidth(2);
    line->Draw("same");
    c->Print(pdfFile);


    MissingMass_hist->Draw("HIST"); 
    TLine *line2 = new TLine(0.496, 0, 0.496, MissingMass_hist->GetMaximum());
    line2->SetLineColor(kRed);     // couleur rouge
    line2->SetLineStyle(2);        // 2 = pointillé
    line2->SetLineWidth(2);
    line2->Draw("same");
    c->Print(pdfFile);

    theta_p_el->Draw("COLZ"); 
    theta_p_el->SetMinimum(0);
    c->Print(pdfFile);

    theta_p_pr->Draw("COLZ"); c->Print(pdfFile);
    theta_p_pip->Draw("COLZ"); c->Print(pdfFile);
    theta_p_pim->Draw("COLZ"); c->Print(pdfFile);

    //MKs_R_Ks_e->Draw("COLZ"); c->Print(pdfFile);
    //MKl_R_Ks_e->Draw("COLZ"); c->Print(pdfFile);
   // MPhi_R_Ks_e->Draw("COLZ"); c->Print(pdfFile);

    //MKs_R_pip_pim->Draw("COLZ"); c->Print(pdfFile);
   // MKl_R_pip_pim->Draw("COLZ"); c->Print(pdfFile);
    //MPhi_R_pip_pim->Draw("COLZ"); c->Print(pdfFile);

    MinvKsKl_W->Draw("COLZ"); 
    MinvKsKl_W->SetMinimum(0);
    c->Print(pdfFile);

    MinvKsKl_Q2->Draw("COLZ"); 
    MinvKsKl_Q2->SetMinimum(0);
    c->Print(pdfFile);

    MinvKsKl_Thetae->Draw("COLZ"); c->Print(pdfFile);

    // --- Page avec texte ---

    c->Clear();
    c->SetFillColor(0);
    c->cd();

    latex.SetTextSize(0.04);
    latex.SetTextAlign(13);

    latex.DrawLatexNDC(0.1, 0.9, "Summary of cuts for the next plots:");
    latex.DrawLatexNDC(0.1, 0.8, "- Add a cut on missing mass : 0.4 < MM < 0.7 GeV");
    //latex.DrawLatexNDC(0.1, 0.7, "- Add a cut on distance (on x y z) of vertex e- and Ks : 1.5 < R_{1} < 5.0 cm");
    //latex.DrawLatexNDC(0.1, 0.6, "- Add a cut on distance (on x y z) of vertex pi+ pi- : 0 < R_{2} < 2 cm");
    //latex.DrawLatexNDC(0.1, 0.5, "- pi+ and pi- need to be in FD");

    c->Print(pdfFile);



    MinvPiPlusMoins_hist2->Draw("HIST");
    TLine *line3 = new TLine(0.496, 0, 0.496, MinvPiPlusMoins_hist2->GetMaximum());
    line3->SetLineColor(kRed);     // couleur rouge
    line3->SetLineStyle(2);        // 2 = pointillé
    line3->SetLineWidth(2);
    line3->Draw("same");
    c->Print(pdfFile);

    //MKs_R_Ks_e2->Draw("COLZ"); c->Print(pdfFile);
    //MKs_R_pip_pim2->Draw("COLZ"); c->Print(pdfFile);

    c->Clear();
    c->SetFillColor(0);
    c->cd();

    latex.SetTextSize(0.04);
    latex.SetTextAlign(13);

    latex.DrawLatexNDC(0.1, 0.9, "Summary of cuts for the next plots:");
    latex.DrawLatexNDC(0.1, 0.8, "- Replace the cut on missing mass by the cut on invariant mass pi+ pi- : ");
    latex.DrawLatexNDC(0.1, 0.7, "  0.3 < M_{inv} < 0.6 GeV");
    

    c->Print(pdfFile);


    MissingMass_hist2->Draw("HIST");
    TLine *line4 = new TLine(0.496, 0, 0.496, MissingMass_hist2->GetMaximum());
    line4->SetLineColor(kRed);     // couleur rouge
    line4->SetLineStyle(2);        // 2 = pointillé
    line4->SetLineWidth(2);
    line4->Draw("same");
    c->Print(pdfFile);

    //MKl_R_Ks_e2->Draw("COLZ"); c->Print(pdfFile);
    //MKl_R_pip_pim2->Draw("COLZ"); c->Print(pdfFile);


    
    c->Clear();
    c->SetFillColor(0);
    c->cd();

    latex.SetTextSize(0.04);
    latex.SetTextAlign(13);

    latex.DrawLatexNDC(0.1, 0.9, "Summary of cuts for the next plots:");
    latex.DrawLatexNDC(0.1, 0.8, "- both cut are present (in invariant mass pi+ pi- and missing mass) ");
    

    c->Print(pdfFile);



    MinvKsKl_hist->Draw("HIST");
    TLine *line5 = new TLine(1.019, 0, 1.019, MinvKsKl_hist->GetMaximum());
    line5->SetLineColor(kGreen);     // couleur rouge
    line5->SetLineStyle(2);        // 2 = pointillé
    line5->SetLineWidth(2);
    line5->Draw("same");
    c->Print(pdfFile);

    //MPhi_R_Ks_e2->Draw("COLZ"); c->Print(pdfFile);
    //MPhi_R_pip_pim2->Draw("COLZ"); c->Print(pdfFile);


    c->Clear();
    c->SetFillColor(0);
    c->cd();

    latex.SetTextSize(0.04);
    latex.SetTextAlign(13);

    latex.DrawLatexNDC(0.1, 0.9, "Summary of cuts for the next plots:");
    latex.DrawLatexNDC(0.1, 0.7, "- All cuts (like invariant mass pi+ pi- and missing mass) + 0.8 < M_{#phi} < 1.2 GeV");

    latex.DrawLatexNDC(0.1, 0.6, "objective : see where the particles reconstructed as a phi went");
    

    c->Print(pdfFile);

    theta_p_el2->Draw("COLZ"); theta_p_el2->SetMinimum(0); c->Print(pdfFile);
    theta_p_pr2->Draw("COLZ"); c->Print(pdfFile);
    theta_p_pip2->Draw("COLZ"); c->Print(pdfFile);
    theta_p_pim2->Draw("COLZ");
    
    c->Print(pdfFile + ")"); // ferme le PDF

    cout << "PDF MC analysis inbending ecris" <<endl;

    cout << "La somme des poids des evenements correspondants a des phi correctement reconstruit est : " << Sweight <<endl;

    // Nettoyage
    f->Close();
    delete f;
    delete c;


}

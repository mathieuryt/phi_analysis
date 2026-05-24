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
#include "TLine.h"
#include "CreateRoot.h"
#include "FillBinning.h"
#include "FitforSigma.h"

double best_chi2 = 1e99;
std::vector<std::vector<double>> bestSigmaDATA;
std::vector<std::vector<double>> bestSigmaMC;

double bestA1, bestB1, bestC1, bestD1;
double bestA2, bestB2, bestC2, bestD2;

bool changebest = false;
int nBestPlot = 0;

double ComputeChi2(double A1, double B1, double C1, double D1,
                   double A2, double B2, double C2, double D2)
{

    auto result = FitforSigma("/Users/mr282803/Documents/analysisRGB/neutron/neutron_analysis/acceptance/HistoMinvbinP_data.root", "/Users/mr282803/Documents/analysisRGB/neutron/neutron_analysis/acceptance/TreeKpKm_MC.root", A1, B1, C1, D1, A2, B2, C2, D2, 20);
    auto sigmaDATA = result.first;
    auto sigmaMC   = result.second;

    for (size_t i = 0; i < sigmaDATA.size(); ++i) { // Vérification : afficher le contenu
        std::cout << "liste des sigma DATA en Km pour le bin en Kp numero " << i << std::endl;
        for (double val : sigmaDATA[i]) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }

    for (size_t i = 0; i < sigmaMC.size(); ++i) { // Vérification : afficher le contenu
        std::cout << "liste des sigma MC en Km pour le bin en Kp numero " << i << std::endl;
        for (double val : sigmaMC[i]) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }

    //CALCUL DE CHI2
    double chi2_tot = 0.0;

    for (size_t i = 0; i < sigmaDATA.size(); ++i) { // Vérification : afficher le contenu

        for (int j = 0; j < sigmaDATA[i].size(); j++) {

            chi2_tot += (sigmaDATA[i][j] - sigmaMC[i][j])*(sigmaDATA[i][j] - sigmaMC[i][j]);
            
        }
        std::cout << std::endl;
    }

    if (chi2_tot < best_chi2) {

        best_chi2 = chi2_tot;

        bestSigmaDATA = sigmaDATA;
        bestSigmaMC   = sigmaMC;

        bestA1 = A1; bestB1 = B1; bestC1 = C1; bestD1 = D1;
        bestA2 = A2; bestB2 = B2; bestC2 = C2; bestD2 = D2;


        std::cout << "New best: "
                  << "chi2_tot = " << best_chi2
                  << " A1 = " << A1
                  << " B1 = " << B1
                  << " C1 = " << C1
                  << std::endl;

        changebest = true;

    }


    return chi2_tot;
}


void FindSmearing()
{
    double s = 0.0;
    gStyle->SetOptStat(0);

    double A1, B1, C1, D1;
    double A2, B2, C2, D2;

    gROOT->Reset();
    gSystem->Load("libPhysics.so");
    gSystem->Load("libCore.so");
    gSystem->Load("libTree.so");
    gSystem->Load("libRIO.so");
    gSystem->Load("libHist.so");
    gStyle->SetPalette(kBird); 

    //createroot("/Users/mr282803/Documents/analysisRGB/neutron/neutron_data2/fall2019_outbending_neutronKpKm_data_v2.root", s, false);
    //createroot("/Users/mr282803/Documents/analysisRGB/neutron/neutron_data2/fall2019_outbending_neutronKpKm_mc_v6.root", s, true);

    TCanvas* c = new TCanvas("c","c",1000,800);

    // ouvrir le PDF
    c->Print("plots.pdf[");

    TFile* f = TFile::Open("HistoMinv_data.root", "READ");
    TH1D* hMinv = (TH1D*)f->Get("hMinv");
    TH2D* h2 = (TH2D*)f->Get("hPkp_vs_Pkm");
    TH2D* h3 = (TH2D*)f->Get("hPkp_vs_Pkm_v2");

    h2->Draw("COLZ");
    c->Print("plots.pdf");

    h3->Draw("COLZ");
    c->Print("plots.pdf");

    hMinv->Draw();
    c->Print("plots.pdf");

    TFile* f_mc = TFile::Open("HistoMinv_mc.root", "READ");
    TH1D* hMinv_mc = (TH1D*)f_mc->Get("hMinv");
    TH2D* h2_mc = (TH2D*)f_mc->Get("hPkp_vs_Pkm");

    TF1* f1 = new TF1("f1", "(23./36.)*x - (5./36.)", 0, 5.2);
    TF1* f2 = new TF1("f2", "(68./47.)*x + (125./235.)", 0, 3.1);
    TF1* f3 = new TF1("f3", "-2.*x + 2.", 0.15, 5.6);
    TF1* f4 = new TF1("f4", "(-8./9.)*x + (66./9.)", 2.5, 5.2);

    h2_mc->Draw("COLZ");

    f1->SetLineColor(kRed);
    f1->SetLineWidth(2);
    f1->SetLineStyle(2);
    f1->Draw("SAME");

    f2->SetLineColor(kRed);
    f2->SetLineWidth(2);
    f2->SetLineStyle(2);
    f2->Draw("SAME");

    f3->SetLineColor(kRed);
    f3->SetLineWidth(2);
    f3->SetLineStyle(2);
    f3->Draw("SAME");

    f4->SetLineColor(kRed);
    f4->SetLineWidth(2);
    f4->SetLineStyle(2);
    f4->Draw("SAME");

    c->Print("plots.pdf");

    hMinv_mc->Draw();
    c->Print("plots.pdf");


    //la fonction fillbinning rempli des histo root par bin de momemtum pour les data (false) et les MC (true).
    // le but est de trouver les paraemttre A,B... (indice 1 pour Kp et 2 pour Km) de la fonction de smearing qui minimise le chi2.

    A1 = B1 = C1 = D1 = 0.05;
    A2 = B2 = C2 = D2 = 0.05;


    fillbinning("/Users/mr282803/Documents/analysisRGB/neutron/neutron_data2/fall2019_outbending_neutronKpKm_data_v2.root", false);
    fillbinning("/Users/mr282803/Documents/analysisRGB/neutron/neutron_data2/fall2019_outbending_neutronKpKm_mc_v6.root", true);

    //double chi2_tot = ComputeChi2(A1, B1, C1, D1, A2, B2, C2, D2);
    //std::cout << "chi2 = " << chi2_tot << std::endl;


    for (double A1 = 0.029; A1 < 0.030; A1 += 0.001) {
        for (double B1 = 0.071; B1 < 0.074; B1 += 0.0005) {
            for (double C1 = 0.055; C1 < 0.058; C1 += 0.0005) {
                for (double D1 = 0.095; D1 < 0.100; D1 += 0.0005) {

                    A2 = A1;
                    B2 = B1;
                    C2 = C1;
                    D2 = D1;

                    std::cout << "\n=====================\n"
                    << Form("-- NEW CALCUL --\n")
                    << "=====================\n\n";


                    std::cout << "Calcul pour pour A1 =  " << A1 << ", B1 = "  << B1 << ", C1 = "<< C1 << ", D1 = " << D1 << "......" << std::endl;
                    double chi2_tot = ComputeChi2(A1, B1, C1, D1,
                                        A2, B2, C2, D2);

                    std::cout << "Chi2 is equal to 1 : " << chi2_tot << std::endl;

                    if(changebest){

                        int nBins = 0;
                        for (size_t i = 0; i < bestSigmaDATA.size(); ++i) {
                            nBins += bestSigmaDATA[i].size();
                        }

                        nBestPlot++;

                        TH1D* hData = new TH1D(Form("hData_best_%d", nBestPlot),
                            "#sigma data et MC par bin",
                            nBins, 0, nBins);

                        TH1D* hMC = new TH1D(Form("hMC_best_%d", nBestPlot),
                            "#sigma data et MC par bin",
                            nBins, 0, nBins);

                        hData->SetDirectory(0);
                        hMC->SetDirectory(0);

                        int binIndex = 1;

                        for (size_t i = 0; i < bestSigmaDATA.size(); ++i) {
                            for (size_t j = 0; j < bestSigmaDATA[i].size(); ++j) {

                            hData->SetBinContent(binIndex, bestSigmaDATA[i][j]);
                            hMC->SetBinContent(binIndex, bestSigmaMC[i][j]);

                             binIndex++;
                            }
                        }

                        TCanvas* cBest = new TCanvas(Form("cBest_%d", nBestPlot),
                             Form("cBest_%d", nBestPlot),
                             1000, 800);

                        cBest->cd();

                        hData->SetMarkerStyle(20);
                        hData->SetMarkerSize(1.2);
                        hData->SetMarkerColor(kBlack);
                        hData->SetLineColor(kBlack);

                        hMC->SetMarkerStyle(21);
                        hMC->SetMarkerSize(1.2);
                        hMC->SetMarkerColor(kRed);
                        hMC->SetLineColor(kRed);

                        hData->GetXaxis()->SetTitle("Bin");
                        hData->GetYaxis()->SetTitle("#sigma");

                        double maxY = std::max(hData->GetMaximum(), hMC->GetMaximum());
                        hData->SetMaximum(2.0 * maxY);
                        hData->SetMinimum(0.0);

                        hData->Draw("P");
                        hMC->Draw("P SAME");

                        TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
                        leg->AddEntry(hData, "Data", "p");
                        leg->AddEntry(hMC, "MC smeared", "p");
                        leg->Draw();

                        TLatex latex;
                        latex.SetNDC();
                        latex.SetTextSize(0.03);
                        latex.DrawLatex(0.15, 0.85, Form("best chi2 = %.6g", best_chi2));
                        latex.DrawLatex(0.15, 0.80, Form("A1=%.4f B1=%.4f C1=%.4f D1=%.4f",
                                  bestA1, bestB1, bestC1, bestD1));

                        cBest->Modified();
                        cBest->Update();
                        cBest->Print("plots.pdf");

                        delete cBest;
                        changebest = false;

                    }

                }


            }
        }
    }

    // fermer le PDF

    f->Close();
    f_mc->Close();

    delete f;
    delete f_mc;
    c->Print("plots.pdf]");


}
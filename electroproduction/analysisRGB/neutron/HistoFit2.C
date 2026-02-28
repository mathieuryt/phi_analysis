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


#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"

using namespace RooFit;

void HistoFit2()
{


    std::vector<double> bornesQ2;
    double valeur; //pr parcourir les fichiers

    std::string name_fichierQ2 = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_analysis/bornes_bins/bornes_Q2.txt";
    std::ifstream fichierQ2(name_fichierQ2);

    while (fichierQ2 >> valeur) {
         bornesQ2.push_back(valeur);
     }
 
    fichierQ2.close();

    int nb_bin_Q2 = bornesQ2.size()-1;


    std::vector<std::vector<double>> liste_de_listes(bornesQ2.size()-1);
    size_t k = bornesQ2.size()-1; // Par exemple, on a 5 fichiers de 0 à 4

    for (size_t i = 0; i < k; ++i) {
        
        std::string name_fichiert = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_analysis/bornes_bins/bornes_t_" + std::to_string(i) + ".txt";
        std::ifstream fichiert(name_fichiert);

        if (fichiert) {
            double valeur;
            while (fichiert >> valeur) {
                liste_de_listes[i].push_back(valeur);
            }
            fichiert.close();
            std::cout << "Fichier " << name_fichiert << " lu avec succès." << std::endl;
        } else {
            std::cerr << "Erreur lors de l'ouverture du fichier : " << name_fichiert << std::endl;
        }
    }

    //liste des integral fit : 

    std::vector<std::vector<double>> liste_de_listes_integral(bornesQ2.size()-1);
    std::vector<std::vector<double>> liste_de_listes_error(bornesQ2.size()-1);

    std::vector<std::vector<double>> liste_de_listes_Nintegral_mc(bornesQ2.size()-1);
    std::vector<std::vector<double>> liste_de_listes_Nerror_mc(bornesQ2.size()-1);

    std::vector<std::vector<double>> liste_de_listes_Nintegral_gen(bornesQ2.size()-1);

    

    // -----------------------
    // Ouvrir le fichier ROOT
    // -----------------------
    TFile *f = TFile::Open("Histos_MinvKpKm_data.root");
    if (!f || f->IsZombie()) {
        cout << "Erreur ouverture fichier" << endl;
        return;
    }

    TFile *fMC = TFile::Open("Histos_MinvKpKm_mc.root");


    //ouvrir le PDF de sortie 

    TCanvas *c = new TCanvas("c", "Fit Minv", 800, 600);
        c->Print("Fits_Minv.pdf[");


    for (int i = 0; i < nb_bin_Q2; i++) {
  
        double Q2_lim1 = bornesQ2[i];
        double Q2_lim2 = bornesQ2[i+1];

            for (int j = 0; j < liste_de_listes[i].size()-1; j++) {

                double tlim_a = liste_de_listes[i][j];
                double tlim_b = liste_de_listes[i][j+1];

                // Nom de l'histogramme
                TString hname = Form("Histogramme_%d_%d", i, j);
                TH1F *h = (TH1F*)f->Get(hname);

                if (!h) {
                    cout << "Histogramme non trouvé : " << hname << endl;
                    return;
                }

                // -----------------------
                // Variable observable
                // -----------------------
                RooRealVar mass("mass", "M(K^{+}K^{-}) [GeV]", 0.992, 1.2);


                // -----------------------
                // Données (binnées)
                // -----------------------
                RooDataHist data("data", "data", mass, Import(*h));

                // -----------------------
                // Signal (Gaussian)
                // -----------------------
                RooRealVar mean("mean", "mean", 1.02, 1.01, 1.03);
                RooRealVar sigma("sigma", "sigma", 0.006, 0.001, 0.007);

                RooGaussian signal("signal", "signal pdf", mass, mean, sigma);

                RooRealVar mK("mK","mK",0.493677);
                mK.setConstant(true);

                // paramètres libres
                RooRealVar a0("a0","a0",0.0,-20,20);
                RooRealVar a1("a1","a1",0.0,-20,20);

                // background physique
                RooGenericPdf background("background",  "sqrt(mass*mass - 4*mK*mK) * exp(a0*mass + a1*mass*mass)", RooArgList(mass,mK,a0,a1));
    
                // -----------------------
                // Rendements
                // -----------------------
                RooRealVar nsig("nsig", "signal yield", h->Integral()*0.5, 0, h->Integral());
                RooRealVar nbkg("nbkg", "background yield", h->Integral()*0.5, 0, h->Integral());

                // -----------------------
                // Modèle total
                // -----------------------
                RooAddPdf model("model", "signal + background",
                    RooArgList(signal, background),
                    RooArgList(nsig, nbkg));

                // -----------------------
                // Fit
                // -----------------------
                RooFitResult* res = model.fitTo(data, Extended(true), Save(), PrintLevel(-1));
                //model.fitTo(data, Extended(true), PrintLevel(-1));

                // -----------------------
                // Plot
                // -----------------------
                RooPlot *frame = mass.frame(Title(" "));
                data.plotOn(frame, Name("myData"));
                model.plotOn(frame, Name("myModel"));
                model.plotOn(frame, Components(background), LineStyle(kDashed), LineColor(kRed));
                model.plotOn(frame, Components(signal), LineStyle(kDashed), LineColor(kBlue));


                
                frame->Draw();
                TLatex latex;
                latex.SetNDC();
                latex.SetTextSize(0.035);
                latex.DrawLatex(0.65, 0.85, Form("Q^{2} #in [%.2f, %.2f] GeV^{2}", Q2_lim1, Q2_lim2));
                latex.DrawLatex(0.65, 0.80,Form("t #in [%.2f, %.2f] GeV^{2}", tlim_a, tlim_b));

                TLatex latex2;
                latex2.SetNDC();              // coordonnées normalisées [0,1]
                latex2.SetTextSize(0.035);    // taille du texte
                latex2.SetTextColor(kRed);    // couleur rouge

                // Position à droite du plot : x ~ 0.65-0.85, y ~ 0.8-0.9
                latex2.DrawLatex(0.65, 0.75, Form("N_{#phi} = %.0f #pm %.0f", nsig.getVal(), nsig.getError()));

                TLatex latex3;
                latex3.SetNDC();
                latex3.SetTextSize(0.035);
                //latex3.DrawLatex(0.75, 0.85, Form("#chi^{2}/ndf = %.2f", frame->chiSquare()));

                int nPar = res->floatParsFinal().getSize();

                double chi2_reduced = frame->chiSquare("myModel", "myData", nPar);

                latex3.DrawLatex(0.65, 0.70, Form("#chi^{2}/ndf = %.2f", chi2_reduced));



                c->Print("Fits_Minv.pdf");

                cout << "Signal yield = " << nsig.getVal()
                << " +/- " << nsig.getError() << endl;

                cout << "Background yield = " << nbkg.getVal()
                << " +/- " << nbkg.getError() << endl;


                liste_de_listes_integral[i].push_back(nsig.getVal());
                liste_de_listes_error[i].push_back(nsig.getError());

                delete frame;


                
                // MC// MC// MC //FIT




                TH1F *hMC = (TH1F*)fMC->Get(hname);
                TH1F *hMC_withBkg = (TH1F*)hMC->Clone(Form("%s_withBkg", hname.Data()));
                int nBkgToGenerate = (int) nbkg.getVal();

                TF1 fbkg("fbkg", "sqrt(x*x - 4*0.493677*0.493677) * exp([0]*x + [1]*x*x)", 0.992, 1.2);

                fbkg.SetParameters(a0.getVal(), a1.getVal());
                hMC_withBkg->FillRandom("fbkg", nBkgToGenerate);

                RooRealVar mass_mc("mass_mc", "M(K^{+}K^{-}) [GeV]", 0.992, 1.2);

                RooDataHist data_mc("data_mc", "data_mc", mass_mc, Import(*hMC_withBkg));

                // -----------------------
                // Signal (Gaussian)
                // -----------------------
                RooRealVar mean_mc("mean_mc", "mean_mc", 1.02, 1.01, 1.03);
                RooRealVar sigma_mc("sigma_mc", "sigma_mc", 0.005, 0.003, 0.01);

                RooGaussian signal_mc("signal_mc", "signal_mc pdf", mass_mc, mean_mc, sigma_mc);

                RooRealVar mK_mc("mK_mc","mK_mc",0.493677);
                mK_mc.setConstant(true);

                // paramètres libres
                RooRealVar a0_mc("a0_mc","a0_mc",0.0,-20,20);
                RooRealVar a1_mc("a1_mc","a1_mc",0.0,-20,20);

                // background physique
                RooGenericPdf background_mc("background_mc",  "sqrt(mass_mc*mass_mc - 4*mK_mc*mK_mc) * exp(a0_mc*mass_mc + a1_mc*mass_mc*mass_mc)", RooArgList(mass_mc,mK_mc,a0_mc,a1_mc));
    
                // -----------------------
                // Rendements
                // -----------------------
                RooRealVar nsig_mc("nsig_mc", "signal_mc yield", hMC_withBkg->Integral()*0.5, 0, hMC_withBkg->Integral());
                RooRealVar nbkg_mc("nbkg_mc", "background_mc yield", hMC_withBkg->Integral()*0.5, 0, hMC_withBkg->Integral());

                // -----------------------
                // Modèle total
                // -----------------------
                RooAddPdf model_mc("model_mc", "signal_mc + background_mc",
                    RooArgList(signal_mc, background_mc),
                    RooArgList(nsig_mc, nbkg_mc));

                // -----------------------
                // Fit
                // -----------------------
                RooFitResult* res_mc = model_mc.fitTo(data_mc, Extended(true), Save(), PrintLevel(-1));
                //model.fitTo(data, Extended(true), PrintLevel(-1));

                // -----------------------
                // Plot
                // -----------------------
                RooPlot *frame_mc = mass_mc.frame(Title(" "));
                data_mc.plotOn(frame_mc, Name("myData_mc"));
                model_mc.plotOn(frame_mc, Name("myModel_mc"));
                model_mc.plotOn(frame_mc, Components(background_mc), LineStyle(kDashed), LineColor(kRed));
                model_mc.plotOn(frame_mc, Components(signal_mc), LineStyle(kDashed), LineColor(kBlue));


                
                frame_mc->Draw();
                TLatex latex_mc;
                latex_mc.SetNDC();
                latex_mc.SetTextSize(0.035);
                latex_mc.DrawLatex(0.65, 0.85, Form("Q^{2} #in [%.2f, %.2f] GeV^{2}", Q2_lim1, Q2_lim2));
                latex_mc.DrawLatex(0.65, 0.80,Form("t #in [%.2f, %.2f] GeV^{2}", tlim_a, tlim_b));

                TLatex latex2_mc;
                latex2_mc.SetNDC();              // coordonnées normalisées [0,1]
                latex2_mc.SetTextSize(0.035);    // taille du texte
                latex2_mc.SetTextColor(kRed);    // couleur rouge

                // Position à droite du plot : x ~ 0.65-0.85, y ~ 0.8-0.9
                latex2_mc.DrawLatex(0.65, 0.75, Form("N_{#phi} MC = %.0f #pm %.0f", nsig_mc.getVal(), nsig_mc.getError()));

                TLatex latex3_mc;
                latex3_mc.SetNDC();
                latex3_mc.SetTextSize(0.035);
                //latex3.DrawLatex(0.75, 0.85, Form("#chi^{2}/ndf = %.2f", frame->chiSquare()));

                int nPar_mc = res_mc->floatParsFinal().getSize();

                double chi2_reduced_mc = frame_mc->chiSquare("myModel_mc", "myData_mc", nPar_mc);

                latex3_mc.DrawLatex(0.65, 0.70, Form("#chi^{2}/ndf MC = %.2f", chi2_reduced_mc));

                liste_de_listes_Nintegral_mc[i].push_back(nsig_mc.getVal());
                liste_de_listes_Nerror_mc[i].push_back(nsig_mc.getError());



                c->Print("Fits_Minv.pdf");

                delete frame_mc;

                TString hname_mc = Form("Histogramme_Ngen_%d_%d", i, j);
                TH1F *hMC_gen = (TH1F*)fMC->Get(hname_mc);
                double Ngen = hMC_gen->GetSumOfWeights();
                std::cout << "Ngen =  " << Ngen << std::endl;
                
                liste_de_listes_Nintegral_gen[i].push_back(Ngen);






            }
    }

    c->Print("Fits_Minv.pdf]");


    //SAUVEGARDE DANS DES FICHIERS N DATA .txt : 

     std::string chemin = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_crosssection/N/";


    for (size_t i = 0; i < liste_de_listes_integral.size(); ++i) {
    // Nom du fichier : integral_t_0.txt, integral_t_1.txt, ...
    std::string filename = chemin + Form("integralN_t_%zu.txt", i);
    std::ofstream outfile(filename);

    if (!outfile) {
        std::cerr << "Erreur ouverture fichier " << filename << std::endl;
        continue;
    }
    // Écrire tous les nsig pour ce Q² en une seule ligne
    for (size_t j = 0; j < liste_de_listes_integral[i].size(); ++j) {
        outfile << liste_de_listes_integral[i][j] << std::endl;
    }
    outfile << std::endl;

    outfile.close();
    std::cout << "Fichier integral N data signal ecris : " << filename << std::endl;
    }



    for (size_t i = 0; i < liste_de_listes_error.size(); ++i) {
    // Nom du fichier : integral_t_0.txt, integral_t_1.txt, ...
    std::string filename2 = chemin + Form("errorN_t_%zu.txt", i);
    std::ofstream outfile2(filename2);

    if (!outfile2) {
        std::cerr << "Erreur ouverture fichier " << filename2 << std::endl;
        continue;
    }
    // Écrire tous les nsig pour ce Q² en une seule ligne
    for (size_t j = 0; j < liste_de_listes_error[i].size(); ++j) {
        outfile2 << liste_de_listes_error[i][j] << std::endl;
    }
    outfile2 << std::endl;

    outfile2.close();
    std::cout << "Fichier erreur sur N data écrit : " << filename2 << std::endl;
    }




    //SAUVGARDE DANS FICHIER Nrec et Ngen MC

    std::string chemin_mc = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_crosssection/A/";


    for (size_t i = 0; i < liste_de_listes_Nintegral_mc.size(); ++i) {
    // Nom du fichier : integral_t_0.txt, integral_t_1.txt, ...
    std::string filename_mc = chemin_mc + Form("integralNrec_mc_t_%zu.txt", i);
    std::ofstream outfile_mc(filename_mc);

    if (!outfile_mc) {
        std::cerr << "Erreur ouverture fichier " << filename_mc << std::endl;
        continue;
    }
    // Écrire tous les nsig pour ce Q² en une seule ligne
    for (size_t j = 0; j < liste_de_listes_Nintegral_mc[i].size(); ++j) {
        outfile_mc << liste_de_listes_Nintegral_mc[i][j] << std::endl;
    }
    outfile_mc << std::endl;

    outfile_mc.close();
    std::cout << "Fichier integral Nrec MC signal ecris : " << filename_mc << std::endl;
    }



    for (size_t i = 0; i < liste_de_listes_Nerror_mc.size(); ++i) {
    // Nom du fichier : integral_t_0.txt, integral_t_1.txt, ...
    std::string filename2_mc = chemin_mc + Form("errorNrec_mc_t_%zu.txt", i);
    std::ofstream outfile2_mc(filename2_mc);

    if (!outfile2_mc) {
        std::cerr << "Erreur ouverture fichier " << filename2_mc << std::endl;
        continue;
    }
    // Écrire tous les nsig pour ce Q² en une seule ligne
    for (size_t j = 0; j < liste_de_listes_Nerror_mc[i].size(); ++j) {
        outfile2_mc << liste_de_listes_Nerror_mc[i][j] << std::endl;
    }
    outfile2_mc << std::endl;

    outfile2_mc.close();
    std::cout << "Fichier erreur sur integral Nrec MC écrit : " << filename2_mc << std::endl;
    }




    for (size_t i = 0; i < liste_de_listes_Nintegral_gen.size(); ++i) {
    // Nom du fichier : integral_t_0.txt, integral_t_1.txt, ...
    std::string filename3_mc = chemin_mc + Form("integralNgen_mc_t_%zu.txt", i);
    std::ofstream outfile3_mc(filename3_mc);

    if (!outfile3_mc) {
        std::cerr << "Erreur ouverture fichier " << filename3_mc << std::endl;
        continue;
    }
    // Écrire tous les nsig pour ce Q² en une seule ligne
    for (size_t j = 0; j < liste_de_listes_Nintegral_gen[i].size(); ++j) {
        outfile3_mc << liste_de_listes_Nintegral_gen[i][j] << std::endl;
    }
    outfile3_mc << std::endl;

    outfile3_mc.close();
    std::cout << "Fichier integral Ngen MC signal ecris : " << filename3_mc << std::endl;
    }

    

}


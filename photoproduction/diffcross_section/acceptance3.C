#include "TH1F.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include <array>
#include <iostream>
#include <vector>
#include <tuple>
#include <TString.h>
#include "Fit/FitResult.h"
#include "Math/MinimizerOptions.h"


void acceptance3() {

    int nb_histo = 5;
    int nb_bin_E_gamma = 5;
    int nb_bin_t = 5;

    double bin = 200;
    double plage = 1;

    double borne_fit_a = 0.95;
    double borne_fit_b = 1.075;

    double borne_plot_a = 0.50;
    double borne_plot_b = 1.5;
    std::vector<float> liste_nb_phi;

    // Charger le fichier ROOT et l'histogramme
    gEnv->SetValue("Browser.Name", "TRootBrowser");
    gStyle->SetOptStat(0); // supr tt les stats attention
    TFile *file = TFile::Open("Histos_M2Q.root");


    for (int u = 0; u < nb_bin_E_gamma; u++) {

        for (int k = 0; k < nb_bin_t; k++) {

        TH1F *hist = (TH1F*)file->Get(Form("Histogramme_%d_%d", u, k));

        // sous gauss
        TF1 *fit_model = new TF1("fit_model","[0]*ROOT::Math::crystalball_pdf(x, [1], [2], [3], [4])", borne_fit_a, borne_fit_b);
        
        // Paramètres initiaux N, alpha, n, sigma, mu
        fit_model->SetParameters(0.002, -2, 2.12, 0.005, 1.017);

        fit_model->SetParLimits(4, 1, 1.05); // Moyenne
        fit_model->SetParLimits(3, 0.001, 0.008); // Sigma ecart type
        fit_model->SetParLimits(2, 1.5, 5); // n decroissance queue 
        fit_model->SetParLimits(1, -5, 5); // alpha position queue
        fit_model->SetParLimits(0, 0.001, 6.8); // N devant crystall

    
        hist->GetXaxis()->SetRangeUser(borne_plot_a, borne_plot_b);
        hist->GetYaxis()->SetRangeUser(0, hist->GetMaximum() * 1.1); 
        TCanvas *c1 = new TCanvas(Form("c1_%d_%d", u, k), Form("Masse Invariante de l'histogramme %d %d", u, k), 800, 600);
        hist->Draw("pe");  // Tracer l'histogramme
    
        // faire le fit
        hist->Fit(fit_model, "R");  // "R" pour ajuster uniquement dans la région de 0.9 à 1.1 GeV
    
        // Afficher le résultat sur un canvas
        fit_model->Draw("same");
        fit_model->SetLineColor(kRed);        // Couleur rouge pour le fit total
    
    
        double integral = fit_model->Integral(borne_fit_a, borne_fit_b);
 
    
        cout << "Intégrale de la gaussienne entre 0.95 et 1.05 GeV = " << integral << endl;

    
        double bin_width = plage/bin; // Intervalle total divisé par le nombre de bins
        double integral_norm = integral / bin_width;

        cout << "Intégrale normalisée = " << integral_norm << endl;


    
    
        // Ajouter les infos du fit
        TPaveText *pt2 = new TPaveText(0.6, 0.7, 0.75, 0.8, "NDC"); 
        pt2->SetTextColor(kBlue);
        pt2->AddText(Form("param0 = %.2f", fit_model->GetParameter(0)));
        pt2->AddText(Form("param1 = %.3f GeV", fit_model->GetParameter(1)));
        pt2->AddText(Form("param2 = %.3f GeV", fit_model->GetParameter(2)));
        pt2->AddText(Form("param3 = %.2f", fit_model->GetParameter(3)));
        pt2->AddText(Form("param4 = %.3f GeV", fit_model->GetParameter(4)));

        cout << "param0 = " << fit_model->GetParameter(0) << endl;
    
        pt2->Draw("same");
    
        // Création du pavé de texte
        TPaveText *pt3 = new TPaveText(0.6, 0.7, 0.75, 0.8, "NDC");
        pt3->SetFillColor(0); // Fond transparent
        pt3->SetTextSize(0.04);
        pt3->SetBorderSize(1);
        pt3->SetTextAlign(12);
    

        double chi2 = fit_model->GetChisquare();
        double ndf = fit_model->GetNDF();
    
    
        // Ajouter les infos du fit
        pt3->SetTextColor(kBlue);
        pt3->AddText(Form("N_{#phi} = %.2f +- %.2f", integral_norm, 0.0));
        pt3->AddText(Form("#chi^{2} reduced = %.2f", chi2/ndf));
        liste_nb_phi.push_back(integral_norm);
      
    
        // Dessiner le pavé
        pt3->Draw("same");


        }

    }


    // On va chercher les GEN
    std::cout << "Nombre de phi REC: ";
    for (float val : liste_nb_phi) std::cout << val << " ";
    std::cout << std::endl;

    // On va chercher les GEN

    std::string nomFichier = "/local/home/mr282803/Documents/irfu/diff_cs_new/AE/nb_phi_GEN.txt";
    std::ifstream fichier(nomFichier);
    std::vector<double> NphiGEN;
 
    double valeur;
    while (fichier >> valeur) {
         NphiGEN.push_back(valeur);
     }
 
    fichier.close();

    std::cout << "Nombre de phi GEN: ";
    for (float val : NphiGEN) std::cout << val << " ";
    std::cout << std::endl;

    // calcul AE avec rec/gen

    std::vector<std::vector<double>> AE(nb_bin_E_gamma);

    double compteur = 0;

    for (int u = 0; u < nb_bin_E_gamma; u++) {

        for (int k = 0; k < nb_bin_t; k++) {

            AE[u].push_back(liste_nb_phi[compteur]/NphiGEN[compteur]);
            compteur += 1;


        }

    }

    // on crrer des fichier par bin Egamma

    for (size_t i = 0; i < AE.size(); ++i) {
        // Les noms de fichiers commencent à 1
        std::string nom_fichier = "AE_bin_E_" + std::to_string(i + 1) + ".txt";
        std::ofstream fichier(nom_fichier);

        if (fichier) {
            for (double valeur : AE[i]) {
                fichier << valeur << std::endl;
            }
            fichier.close();
            std::cout << "Fichier " << nom_fichier << " écrit avec succès." << std::endl;
        } else {
            std::cerr << "Erreur lors de l'ouverture du fichier : " << nom_fichier << std::endl;
        }
    }


    
}
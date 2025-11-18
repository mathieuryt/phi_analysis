#include "TH1F.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
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

void fit2() {
    int nb_histo = 5;
    int nb_bin_E_gamma = 5;
    int nb_bin_t = 5;

    double bin = 100;
    double plage = 0.4;

    double borne_fit_a = 0.90;
    double borne_fit_b = 1.18;

    double borne_plot_a = 0.8;
    double borne_plot_b = 1.2;

    std::vector<std::vector<double>> liste_nb_phi(nb_bin_E_gamma);
    std::vector<std::vector<double>> liste_error_phi(nb_bin_E_gamma);

    // Charger le fichier ROOT et l'histogramme
    gEnv->SetValue("Browser.Name", "TRootBrowser");
    gStyle->SetOptStat(0); // supr tt les stats attention
    TFile *file = TFile::Open("Histos_M2Q.root");

    for (int u = 0; u < nb_bin_E_gamma; u++) {

        for (int k = 0; k < nb_bin_t; k++) {

            TH1F *hist = (TH1F*)file->Get(Form("Histogramme_%d_%d", u, k));
        
            // Fonction combinée : fond quadratique + Crystal Ball
            TF1 *fit_model = new TF1("fit_model",
                "[7]*x*x + [0]*x + [1] + [2]*ROOT::Math::crystalball_pdf(x, [3], [4], [5], [6])",
                borne_fit_a, borne_fit_b);
        
            fit_model->SetParameters(1, 200, 2.05, -1, 2, 0.017, 1.014, 1);
        
            fit_model->SetParLimits(6, 1, 1.05);   // mu
            fit_model->SetParLimits(5, 0.01, 0.03); // sigma
            fit_model->SetParLimits(4, 1.5, 2.5);   // n
            fit_model->SetParLimits(3, -5, -0.5);   // alpha
            fit_model->SetParLimits(2, 1.0, 3.0);   // amplitude CB
        
            // === SETUP Minuit2 ===
            ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
            ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);
            ROOT::Math::MinimizerOptions::SetDefaultErrorDef(1.0);
        
            hist->GetXaxis()->SetRangeUser(borne_plot_a, borne_plot_b);
            hist->GetYaxis()->SetRangeUser(0, hist->GetMaximum() * 1.1);
        
            TCanvas *c1 = new TCanvas(Form("c1_%d_%d", u, k),
                                      Form("Masse Invariante de l'histogramme %d %d", u, k), 800, 600);
            hist->Draw("pe");
        
            // === FIT avec stockage du résultat ===
            TFitResultPtr fit_result = hist->Fit(fit_model, "RS");
        
            // === Forcer Minos sur les paramètres 2 à 6 ===
            TVirtualFitter *fitter = TVirtualFitter::GetFitter();
            if (fitter) {
                for (int i = 2; i <= 6; ++i) {
                    fitter->GetMinuit()->mnmnos(); // force Minos globalement (optionnel)
                    double errLow = 0, errHigh = 0, eparab, gcc;
                    ((TMinuit*)fitter)->mnerrs(i, errHigh, errLow, eparab, gcc);
                    std::cout << "Param " << i << ": erreur Minos - = " << errLow << ", + = " << errHigh << std::endl;
                }
            }
        
            // === Reconstruire les fonctions séparées ===
            TF1 *background = new TF1("background", "[2]*x*x + [0]*x + [1]", borne_fit_a, borne_fit_b);
            background->SetParameters(fit_model->GetParameter(0),
                                      fit_model->GetParameter(1),
                                      fit_model->GetParameter(7));
        
            TF1 *crystal = new TF1("crystal", "[0]*ROOT::Math::crystalball_pdf(x, [1], [2], [3], [4])", borne_fit_a, borne_fit_b);

            crystal->SetParameters(fit_model->GetParameter(2),
                                   fit_model->GetParameter(3),
                                   fit_model->GetParameter(4),
                                   fit_model->GetParameter(5),
                                   fit_model->GetParameter(6));
        
            // === Affichage ===
            fit_model->Draw("same");
            crystal->SetLineColor(kBlue);
            background->SetLineColor(kGreen);
            crystal->Draw("same");
            background->Draw("same");
        
            fit_model->SetLineColor(kRed);
        
            TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
            legend->AddEntry(hist, "Histogramme", "l");
            legend->AddEntry(fit_model, "Fit total", "l");
            legend->AddEntry(crystal, "Signal (Crystal Ball)", "l");
            legend->AddEntry(background, "Background", "l");
            legend->SetTextSize(0.035);
            legend->Draw("same");



            // Paramètres d'origine du signal
            std::vector<double> params_cb(5);
            for (int i = 0; i < 5; i++) {
                params_cb[i] = fit_result->Parameter(i + 2);
            }

            // Intégrale centrale (valeur nominale)
            double integral_central = crystal->Integral(borne_fit_a, borne_fit_b);

            // Vecteurs pour stocker les intégrales ±
            std::vector<double> integrals_up(5), integrals_down(5);

            // On évalue l’impact de chaque paramètre ± erreur Minos
            for (int i = 0; i < 5; i++) {
                double err_low = fit_result->LowerError(i + 2);
                double err_up  = fit_result->UpperError(i + 2);

                // + erreur Minos
                std::vector<double> params_plus = params_cb;
                params_plus[i] += err_up;
                crystal->SetParameters(params_plus.data());
                integrals_up[i] = crystal->Integral(borne_fit_a, borne_fit_b);

                // - erreur Minos
                std::vector<double> params_minus = params_cb;
                params_minus[i] += err_low;
                crystal->SetParameters(params_minus.data());
                integrals_down[i] = crystal->Integral(borne_fit_a, borne_fit_b);
            }

            // Écart max en + et en - (approximation conservative)
            double err_up_integral = 0;
            double err_dn_integral = 0;
            for (int i = 0; i < 5; i++) {
                err_up_integral += std::pow(integrals_up[i] - integral_central, 2);
                err_dn_integral += std::pow(integrals_down[i] - integral_central, 2);
            }
            err_up_integral = std::sqrt(err_up_integral);
            err_dn_integral = std::sqrt(err_dn_integral);

            // Affichage final
            cout << "Intégrale du signal = " << integral_central << endl;
            cout << "Erreur Minos - = " << err_dn_integral << ", Erreur Minos + = " << err_up_integral << endl;

            // Option : normalisée
            double bin_width = plage / bin;
            cout << "Intégrale normalisée = " << integral_central / bin_width
                << "  +" << err_up_integral / bin_width
                << "  -" << err_dn_integral / bin_width << endl;

        
            double bin_width = plage / bin;
            double integral_norm = integral / bin_width;
            double error_integral_norm = error_integral / bin_width;
        
            std::cout << "Intégrale de la Crystal Ball = " << integral
                      << " ± " << error_integral << std::endl;
            std::cout << "Intégrale normalisée = " << integral_norm
                      << " ± " << error_integral_norm << std::endl;
    




        
            // Création du pavé de texte
            TPaveText *pt = new TPaveText(0.6, 0.7, 0.75, 0.8, "NDC"); 
            pt->SetFillColor(0); // Fond transparent
            pt->SetTextSize(0.02);
            pt->SetBorderSize(1);
            pt->SetTextAlign(12);
        
            // Ajouter les infos du fit
            pt->AddText("Parametres background"); // Titre du pavé
            pt->SetTextColor(kGreen);
            pt->AddText(Form("a = %.2f", fit_model->GetParameter(0))); 
            pt->AddText(Form("b = %.2f", fit_model->GetParameter(1)));
        
        
            // Dessiner le pavé
            //pt->Draw("same");
        
            // Création du pavé de texte
            TPaveText *pt2 = new TPaveText(0.6, 0.7, 0.75, 0.8, "NDC");
            pt2->SetFillColor(0); // Fond transparent
            pt2->SetTextSize(0.02);
            pt2->SetBorderSize(1);
            pt2->SetTextAlign(12);
            pt2->AddText("Parametres signal"); // Titre du pavé
        
        
            // Ajouter les infos du fit
            pt2->SetTextColor(kBlue);
            pt2->AddText(Form("param0 = %.2f", fit_model->GetParameter(0)));
            pt2->AddText(Form("param1 = %.2f", fit_model->GetParameter(1)));
            pt2->AddText(Form("param2 = %.2f", fit_model->GetParameter(2)));
            pt2->AddText(Form("pram3 = %.3f GeV", fit_model->GetParameter(3)));
            pt2->AddText(Form("param4 = %.3f GeV", fit_model->GetParameter(4)));
            pt2->AddText(Form("param5 = %.3f GeV", fit_model->GetParameter(5)));
            pt2->AddText(Form("param6 = %.3f GeV", fit_model->GetParameter(6)));
            pt2->AddText(Form("param7 = %.2f", fit_model->GetParameter(7)));
        
            // Dessiner le pavé
            pt2->Draw("same");
    
            double chi2 = fit_model->GetChisquare();
            double ndf = fit_model->GetNDF();
        
            // Création du pavé de texte
            TPaveText *pt3 = new TPaveText(0.6, 0.7, 0.75, 0.8, "NDC");
            pt3->SetFillColor(0); // Fond transparent
            pt3->SetTextSize(0.035);
            pt3->SetBorderSize(1);
            pt3->SetTextAlign(12);
            //pt3->AddText("Nb de phi"); // Titre du pavé
        
        
            // Ajouter les infos du fit
            pt3->SetTextColor(kBlue);
            //pt3->AddText(Form("Integral simple = %.2f +- %.2f", integral, 0.0));
            pt3->AddText(Form("N_{#phi} = %.2f +- %.2f", integral_norm, error_integral_norm));
            pt3->AddText(Form("#chi^{2} reduced = %.2f", chi2/ndf));
           
            // Dessiner le pavé
            pt3->Draw("same");

            liste_nb_phi[u].push_back(integral_norm);
            liste_error_phi[u].push_back(error_integral_norm);

        }
    
    }



    for (size_t i = 0; i < liste_nb_phi.size(); ++i) {
        // Les noms de fichiers commencent à 1
        std::string nom_fichier = "nb_phi_bin_E_" + std::to_string(i + 1) + ".txt";
        std::ofstream fichier(nom_fichier);

        if (fichier) {
            for (double valeur : liste_nb_phi[i]) {
                fichier << valeur << std::endl;
            }
            fichier.close();
            std::cout << "Fichier " << nom_fichier << " écrit avec succès." << std::endl;
        } else {
            std::cerr << "Erreur lors de l'ouverture du fichier : " << nom_fichier << std::endl;
        }
    }

    for (size_t i = 0; i < liste_error_phi.size(); ++i) {
        // Les noms de fichiers commencent à 1
        std::string nom_fichier2 = "error_phi_bin_E_" + std::to_string(i + 1) + ".txt";
        std::ofstream fichier2(nom_fichier2);

        if (fichier2) {
            for (double valeur : liste_error_phi[i]) {
                fichier2 << valeur << std::endl;
            }
            fichier2.close();
            std::cout << "Fichier " << nom_fichier2 << " écrit avec succès." << std::endl;
        } else {
            std::cerr << "Erreur lors de l'ouverture du fichier : " << nom_fichier2 << std::endl;
        }
    }


    
}
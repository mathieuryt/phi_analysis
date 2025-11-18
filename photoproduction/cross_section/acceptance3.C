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
    int nb_histo = 12;

    double bin = 200;
    double plage = 1;

    double borne_fit_a = 0.90;
    double borne_fit_b = 1.18;

    double borne_plot_a = 0.50;
    double borne_plot_b = 1.5;
    std::vector<float> liste_nb_phi;

    // Charger le fichier ROOT et l'histogramme
    gEnv->SetValue("Browser.Name", "TRootBrowser");
    gStyle->SetOptStat(0); // supr tt les stats attention
    TFile *file = TFile::Open("Histos_M2_REC.root");

    for (int u = 0; u < nb_histo; u++) {

        TH1F *hist = (TH1F*)file->Get(Form("Histogramme numero %d", u));

    
        // sous gauss
        TF1 *fit_model = new TF1("fit_model","[0]*ROOT::Math::crystalball_pdf(x, [1], [2], [3], [4])", borne_fit_a, borne_fit_b);
        
        // Paramètres initiaux N, alpha, n, sigma, mu
        fit_model->SetParameters(0.0006, 1.1, 1.7, 0.01, 1.019);
    
    
        hist->GetXaxis()->SetRangeUser(borne_plot_a, borne_plot_b);
        hist->GetYaxis()->SetRangeUser(0, hist->GetMaximum() * 1.1); 

        TCanvas *c1 = new TCanvas(Form("c1_%d", u), Form("Masse Invariante de l'histogramme %d", u), 800, 600);
        hist->Draw("hist");  // Tracer l'histogramme
    
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
        //TPaveText *pt2 = new TPaveText(0.6, 0.7, 0.75, 0.8, "NDC"); 
        //pt2->SetTextColor(kBlue);
        //pt2->AddText(Form("param0 = %.2f", fit_model->GetParameter(0)));
        //pt2->AddText(Form("param1 = %.3f GeV", fit_model->GetParameter(1)));
        //pt2->AddText(Form("param2 = %.3f GeV", fit_model->GetParameter(2)));
        //pt2->AddText(Form("param3 = %.2f", fit_model->GetParameter(3)));
        //pt2->AddText(Form("param4 = %.3f GeV", fit_model->GetParameter(4)));

        cout << "param0 = " << fit_model->GetParameter(0) << endl;
    
        //pt2->Draw("same");
    
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
        pt3->AddText(Form("N_{#phi} = %.2f", integral_norm));
        pt3->AddText(Form("#chi^{2} reduced = %.2f", chi2/ndf));
        liste_nb_phi.push_back(integral_norm);
      
    
        // Dessiner le pavé
        pt3->Draw("same");
    
    }

    // Affichage du contenu du vecteur
    std::cout << "Nombre de phi : ";
    for (float val : liste_nb_phi) std::cout << val << " ";
    std::cout << std::endl;


    std::ofstream fichier("nb_phi.txt");

    if (fichier) {
        for (float valeur : liste_nb_phi) {
            fichier << valeur << std::endl;
        }
        fichier.close();
        std::cout << "Fichier nb de phi écrit avec succès !" << std::endl;
    } else {
        std::cerr << "Erreur lors de l'ouverture du fichier !" << std::endl;
    }


    
}
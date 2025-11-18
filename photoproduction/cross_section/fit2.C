#include "TH1F.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"

void fit2() {
    int nb_histo = 12;

    double bin = 200;
    double plage = 1;

    double borne_fit_a = 0.90;
    double borne_fit_b = 1.18;

    double borne_plot_a = 0.50;
    double borne_plot_b = 1.5;
    std::vector<float> liste_nb_phi;
    std::vector<float> liste_phi_error;

    // Charger le fichier ROOT et l'histogramme
    gEnv->SetValue("Browser.Name", "TRootBrowser");
    gStyle->SetOptStat(0); // supr tt les stats attention
    TFile *file = TFile::Open("Histos_M2Q.root");

    for (int u = 0; u < nb_histo; u++) {

        TH1F *hist = (TH1F*)file->Get(Form("Histogramme numero %d", u));

        //fonction combinée : fond + gaussienne
        TF1 *fit_model = new TF1("fit_model", "[7]*x*x + [0]*x + [1] + [2]*ROOT::Math::crystalball_pdf(x, [3], [4], [5], [6])", borne_fit_a, borne_fit_b);
    
        fit_model->SetParameters(0, 0, 3.37, -2.618, 3.305, 0.019, 1.017, 1);  // Valeurs initiales pour ax+b et la Gaussienne

        if (u==0) {

            fit_model->SetParameters(-1300, 972, 5.14, -4.498, 3, 0.005, 1.019, 461);

        }

        if (u==8) {

            fit_model->SetParameters(1200, -540, 4.37, -4.618, 4.305, 0.015, 1.017, -650);

        }

        if (u==7) {

            fit_model->SetParameters(0, 0, 3.37, -2.618, 3.305, 0.019, 1.017, 1);

        }

        fit_model->SetParLimits(6, 1, 1.05); // Moyenne
        fit_model->SetParLimits(5, 0.001, 0.03); // Sigma ecart type
        fit_model->SetParLimits(4, 1.5, 5); // n decroissance queue 
        fit_model->SetParLimits(3, -5, -1); // alpha position queue
        fit_model->SetParLimits(2, 1.0, 6.8); // N devant crystall

        if (u==11) {

            fit_model->SetParLimits(3, -5, 5); // alpha position queue
            fit_model->SetParameters(-1500, 700, 3.37, 2.618, 3.305, 0.019, 1.017, 900); 

        }
        

        ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);


        hist->GetXaxis()->SetRangeUser(borne_plot_a, borne_plot_b);
        hist->GetYaxis()->SetRangeUser(0, hist->GetMaximum() * 1.1); 
        TCanvas *c1 = new TCanvas(Form("c1_%d", u), Form("Masse Invariante de l'histogramme %d", u), 800, 600);
        hist->Draw("hist");  // Tracer l'histogramme
    
        // faire le fit
        TFitResultPtr fit_result = hist->Fit(fit_model, "RS");  // "R" pour ajuster uniquement dans la région de 0.9 à 1.1 GeV
    
        // fonction background (ax+b)
        TF1 *background = new TF1("background", " [2]*x*x + [0]*x + [1]", borne_fit_a, borne_fit_b);  // ax+b entre 0.9 et 1.1 GeV
        background->SetParameters(fit_model->GetParameter(0), fit_model->GetParameter(1), fit_model->GetParameter(7));  // Valeurs initiales pour a et b
    
        // fonction gaussienne (signal)
        TF1 *crystal = new TF1("crystal", "[0]*ROOT::Math::crystalball_pdf(x, [1], [2], [3], [4])", borne_fit_a, borne_fit_b);
        crystal->SetParameters(fit_model->GetParameter(2), fit_model->GetParameter(3), fit_model->GetParameter(4), fit_model->GetParameter(5), fit_model->GetParameter(6));  // Valeurs initiales pour A, mu, sigma
    
        // Afficher le résultat sur un canvas
        fit_model->Draw("same");
        crystal->Draw("same");
        background->Draw("same"); // Ajouter le fit
    
        // Afficher la légende
        TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->AddEntry(hist, "Histogramme", "l");
        legend->AddEntry(fit_model, "Fit total", "l");
        legend->AddEntry(crystal, "Signal", "l");
        legend->AddEntry(background, "Background", "l");
        legend->SetTextSize(0.035);
    
        legend->Draw("same");
    
        fit_model->SetLineColor(kRed);        // Couleur rouge pour le fit total
        crystal->SetLineColor(kBlue);      // Couleur bleue pour la gaussienne
        background->SetLineColor(kGreen);
    
    
    
        double integral = crystal->Integral(borne_fit_a, borne_fit_b);


        TMatrixD cov_cb = fit_result->GetCovarianceMatrix().GetSub(2,6,2,6);
        Double_t error_integral = crystal->IntegralError(borne_fit_a, borne_fit_b, fit_result->GetParams() + 2, cov_cb.GetMatrixArray());


    
        cout << "Intégrale de la gaussienne entre 0.95 et 1.05 GeV = " << integral << endl;

    
        double bin_width = plage/bin; // Intervalle total divisé par le nombre de bins
        double integral_norm = integral / bin_width;
        double error_integral_norm = error_integral / bin_width;


        cout << "Intégrale normalisée = " << integral_norm << endl;

    
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
        pt2->AddText(Form("param2 = %.2f", fit_model->GetParameter(2)));
        pt2->AddText(Form("pram3 = %.3f GeV", fit_model->GetParameter(3)));
        pt2->AddText(Form("param4 = %.3f GeV", fit_model->GetParameter(4)));
        pt2->AddText(Form("param5 = %.3f GeV", fit_model->GetParameter(5)));
        pt2->AddText(Form("param6 = %.3f GeV", fit_model->GetParameter(6)));
    
        // Dessiner le pavé
        //pt2->Draw("same");

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
        liste_nb_phi.push_back(integral_norm);
        liste_phi_error.push_back(error_integral_norm);

    
        // Dessiner le pavé
        pt3->Draw("same");

    
    }

    // Affichage du contenu du vecteur
    std::cout << "Nombre de phi : ";
    for (float val : liste_nb_phi) std::cout << val << " ";
    std::cout << std::endl;

    std::cout << "Erreur nombre de phi : ";
    for (float val : liste_phi_error) std::cout << val << " ";
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



    std::ofstream fichier2("err_nb_phi.txt");

    if (fichier2) {
        for (float valeur : liste_phi_error) {
            fichier2 << valeur << std::endl;
        }
        fichier2.close();
        std::cout << "Fichier ERR nb de phi écrit avec succès !" << std::endl;
    } else {
        std::cerr << "Erreur lors de l'ouverture du fichier !" << std::endl;
    }

    
}
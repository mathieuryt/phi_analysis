#ifndef FITFORSIGMA_H
#define FITFORSIGMA_H

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

using namespace RooFit;

double Sigma(double x, double p1, double p2, double p3, double p4, double p5,
         double A, double B, double C, double D)
{
    if (x < p1 || x > p5)
        throw std::out_of_range("x hors des bornes");

    if (x < p2) return A;
    if (x < p3) return B;
    if (x < p4) return C;
    return D;
}

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
FitforSigma(const char* chemin, const char* chemin2, double A1, double B1, double C1, double D1, double A2, double B2, double C2, double D2, double nb_randfitting, bool isplot)
{

    cout << "Fitting ... " << endl;
    cout << "Valeurs A1 B1 C1 D1 : " << A1 << ", " << B1 << ", " << C1 << ", " << D1 << endl;
    cout << "Valeurs A2 B2 C2 D2 : " << A2 << ", " << B2 << ", " << C2 << ", " << D2 << endl;
  

    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

    std::vector<double> bornesKp;


    std::string name_fichierQ2 = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_analysis/acceptance/binningPKpKm/bornesKp.txt";
    std::ifstream fichierQ2(name_fichierQ2);

    double valeur;
    while (fichierQ2 >> valeur) {
         bornesKp.push_back(valeur);
     }
 
    fichierQ2.close();

    std::vector<std::vector<double>> liste_de_listes_Km(bornesKp.size()-1);

    size_t sizeKp = bornesKp.size()-1; // 5 nombre de le fichier donc 4 bins

    for (size_t i = 0; i < sizeKp; ++i) {
        
        std::string name_fichiert = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_analysis/acceptance/binningPKpKm/bornesKm_" + std::to_string(i) + ".txt";
        std::ifstream fichiert(name_fichiert);

        if (fichiert) {
            double valeur;
            while (fichiert >> valeur) {
                liste_de_listes_Km[i].push_back(valeur);
            }
            fichiert.close();
            //std::cout << "Fichier " << name_fichiert << " lu avec succès." << std::endl;
        } else {
            std::cerr << "Erreur lors de l'ouverture du fichier : " << name_fichiert << std::endl;
        }
    }

    //std::cout << "Liste des bornes en Kp (phase de FIT): "; // Affichage du contenu du vecteur
    //for (float val : bornesKp) std::cout << val << " ";
    //std::cout << std::endl;

    //for (size_t i = 0; i < liste_de_listes_Km.size(); ++i) { // Vérification : afficher le contenu
        //std::cout << "liste des bornes bins en Km pour le bin en Kp numero (phase de FIT) " << i << std::endl;
        //for (double val : liste_de_listes_Km[i]) {
            //std::cout << val << " ";
        //}
        //std::cout << std::endl;
    //}



    //liste des integral fit : 

    std::vector<std::vector<double>> liste_de_listes_sigmaDATA(bornesKp.size()-1);
    //std::vector<std::vector<double>> liste_de_listes_error(bornesQ2.size()-1);

    //std::vector<std::vector<double>> liste_de_listes_Nrec(bornesQ2.size()-1);
    //std::vector<std::vector<double>> liste_de_listes_Ngen(bornesQ2.size()-1);
    std::vector<std::vector<double>> liste_de_listes_sigmaMC(bornesKp.size()-1);
    //std::vector<std::vector<double>> liste_de_listes_Nrandom(bornesQ2.size()-1);

    //std::vector<std::vector<double>> liste_de_listes_Acc(bornesQ2.size()-1);



    

    // -----------------------
    // Ouvrir le fichier ROOT
    // -----------------------
    TFile *f = TFile::Open(chemin);
    if (!f || f->IsZombie()) {
        cout << "Erreur ouverture fichier" << endl;
    }

    TFile *fMC = TFile::Open(chemin2);

    TFile *f_Nfit = new TFile("NfitHisto.root", "RECREATE");


    //ouvrir le PDF de sortie 
    TCanvas *cfit = nullptr;
    
    if(isplot){

        cfit = new TCanvas("cfit", "Fit Minv", 800, 600);
        cfit->Print("Fits_Minv_Smearing.pdf[");

    }



    double Nphi_tot = 0;


    for (int i = 0; i < sizeKp; i++) {

        double Kp_lim1 = bornesKp[i];
        double Kp_lim2 = bornesKp[i+1];

        for (int j = 0; j < liste_de_listes_Km[i].size()-1; j++) {

            double Km_lim1 = liste_de_listes_Km[i][j];
            double Km_lim2 = liste_de_listes_Km[i][j+1];



            //std::cout << "\n=====================\n"
            //<< Form("-- new bin (i=%d, j=%d) --\n", i, j)
            //<< "=====================\n\n";


            // Nom de l'histogramme
            TString hname = Form("Histogramme_%d_%d", i, j);
            TH1F *h = (TH1F*)f->Get(hname);

            if (!h) {
                cout << "Histogramme non trouvé : " << hname << endl;
            
            }

            // -----------------------
            // Variable observable
            // -----------------------
            RooRealVar mass("mass", "M(K^{+}K^{-}) [GeV]", 0.992, 1.1);


            // -----------------------
            // Données (binnées)
            // -----------------------
            RooDataHist data("data", "data", mass, Import(*h));

            // -----------------------
            // Signal (Gaussian)
            // -----------------------
            RooRealVar mean("mean", "mean", 1.02, 1.01, 1.03);
            RooRealVar sigma("sigma", "sigma", 0.006, 0.002, 0.009);

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
            latex.DrawLatex(0.65, 0.85, Form("P_{Kp} #in [%.2f, %.2f] GeV^{2}", Kp_lim1, Kp_lim2));
            latex.DrawLatex(0.65, 0.80,Form("P_{Km} #in [%.2f, %.2f] GeV^{2}", Km_lim1, Km_lim2));

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

            TLatex latex4;
            latex4.SetNDC();
            latex4.SetTextSize(0.035);
            latex4.SetTextColor(kGreen);
            latex4.DrawLatex(0.65, 0.65, Form("#sigma = %.6f", sigma.getVal()));


            if(isplot){
                cfit->Print("Fits_Minv_Smearing.pdf");
            }
            //cout << "Signal yield = " << nsig.getVal()
            //<< " +/- " << nsig.getError() << endl;

            //cout << "Background yield = " << nbkg.getVal()
            //<< " +/- " << nbkg.getError() << endl;


            liste_de_listes_sigmaDATA[i].push_back(sigma.getVal());
            Nphi_tot += nsig.getVal();
            
            //liste_de_listes_error[i].push_back(nsig.getError());

            delete frame;



            // MC// MC// MC //FIT

            double m_kaons = 0.4937;

            
            int nSignToGenerate = (int) nsig.getVal();
            TString tname = Form("Tree_%d_%d", i, j);
            TTree* tMC = (TTree*)fMC->Get(tname);

            TLorentzVector* Kp = nullptr;
            TLorentzVector* Km = nullptr;

            tMC->SetBranchAddress("Kp", &Kp);
            tMC->SetBranchAddress("Km", &Km);

            Long64_t nTreeEntries = tMC->GetEntries();

            std::vector<TLorentzVector> vKp, vKm;

            for (Long64_t e = 0; e < nTreeEntries; ++e) {
                 tMC->GetEntry(e);
                 vKp.push_back(*Kp);
                 vKm.push_back(*Km);
            }

            double nbfill = 25000;
            //double scale = 1/nbfill;
            
            TH1F *stocker_Nfit = new TH1F(Form("stocker_Nfit_%d_%d", i, j), "Nf", 75, 0, 0.01);

            for(int k = 0; k < nb_randfitting; k++){


                //TH1F *hMC_signalonly_smeared = (TH1F*)h->Clone(Form("%s_withBkg", hname.Data())); // lui on va le smearer 
                //hMC_signalonly_smeared->Reset();
                //hMC_signalonly_smeared->Sumw2();
                TH1F *hSignal = (TH1F*)h->Clone(Form("hSignal_%d_%d_%d", i, j, k));
                TH1F *hBkg    = (TH1F*)h->Clone(Form("hBkg_%d_%d_%d", i, j, k));
                TH1F *hTot    = (TH1F*)h->Clone(Form("hTot_%d_%d_%d", i, j, k));

                hSignal->Reset();
                hBkg->Reset();
                hTot->Reset();

                hSignal->Sumw2();
                hBkg->Sumw2();
                hTot->Sumw2();

                ULong_t seed = 12345 + 100000*i + 1000*j + k;

                gRandom->SetSeed(seed);
                TRandom3 rng(seed);


                for (int n = 0; n < nbfill; ++n) {

                    Long64_t entry = rng.Integer(vKp.size());

                    const TLorentzVector& Kp = vKp[entry];
                    const TLorentzVector& Km = vKm[entry];
                    
                    
                    double sigmaSmearKp = Sigma(Kp.P(), 0, 1.2, 1.6, 2.2, 6.0, A1, B1, C1, D1);
                    double sigmaSmearKm = Sigma(Km.P(), 0, 1.2, 1.6, 2.2, 6.0, A2, B2, C2, D2);

                    double Pkp_smear = Kp.P() * rng.Gaus(1.0, sigmaSmearKp);
                    double Pkm_smear = Km.P() * rng.Gaus(1.0, sigmaSmearKm);

                    TVector3 Kp_smear_vect, Km_smear_vect;
                        
                    Kp_smear_vect.SetMagThetaPhi(Pkp_smear, Kp.Theta(), Kp.Phi());
                    Km_smear_vect.SetMagThetaPhi(Pkm_smear, Km.Theta(), Km.Phi());

                    TLorentzVector Kp_smear, Km_smear;

                    Kp_smear.SetVectM(Kp_smear_vect, m_kaons);
                    Km_smear.SetVectM(Km_smear_vect, m_kaons);        

                    double m_smear = (Kp_smear + Km_smear).M();

                    hSignal->Fill(m_smear);  // fill avec un nb egale au data integral
    
                }


                int nBkgToGenerate = (int) nbkg.getVal();

                TF1 fbkg("fbkg", "sqrt(x*x - 4*0.493677*0.493677) * exp([0]*x + [1]*x*x)", 0.992, 1.1);

                fbkg.SetParameters(a0.getVal(), a1.getVal());

                hBkg->FillRandom("fbkg", nbfill);// la on a un histo avec un background fill avec autant de donnée que l'integral background des DATAs

                //hMC_signalonly_smeared->Scale(scale);
                // on utilise hMC pour fill randome notre signal car hMC n'est pas exactement une gaussiene donc on introduit pas un biais (biais si on prenant la fonction signal des data qui est une pur gauss par ex)
                
                hSignal->Scale(nSignToGenerate / hSignal->Integral());
                hBkg->Scale(nBkgToGenerate / hBkg->Integral());

                hTot->Add(hSignal);
                hTot->Add(hBkg);

                if(k == 10){

                    if(isplot){
                        hSignal->Draw();
                        cfit->Print("Fits_Minv_Smearing.pdf");
                        hBkg->Draw();
                        cfit->Print("Fits_Minv_Smearing.pdf");
                        hTot->Draw();
                        cfit->Print("Fits_Minv_Smearing.pdf");
                    }

                }

                RooRealVar mass_mc("mass_mc", "M(K^{+}K^{-}) [GeV]", 0.992, 1.1);

                RooDataHist data_mc("data_mc", "data_mc", mass_mc, Import(*hTot));

                // -----------------------
                // Signal (Gaussian)
                // -----------------------
                RooRealVar mean_mc("mean_mc", "mean_mc", 1.02, 1.01, 1.03);
                RooRealVar sigma_mc("sigma_mc", "sigma_mc", 0.0055, 0.0020, 0.016); //debut , min , max

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
                RooRealVar nsig_mc("nsig_mc", "signal_mc yield", hTot->Integral()*0.5, 0, 2*hTot->Integral());
                RooRealVar nbkg_mc("nbkg_mc", "background_mc yield", hTot->Integral()*0.5, 0, 2*hTot->Integral());

                // -----------------------
                // Modèle total
                // -----------------------
                RooAddPdf model_mc("model_mc", "signal_mc + background_mc",
                    RooArgList(signal_mc, background_mc),
                    RooArgList(nsig_mc, nbkg_mc));

                // -----------------------
                // Fit
                // -----------------------
                RooFitResult* res_mc = model_mc.fitTo(data_mc, Extended(true), Save(), PrintLevel(-1), SumW2Error(true));
                //model.fitTo(data, Extended(true), PrintLevel(-1));

                // -----------------------
                // Plot
                // -----------------------
                RooPlot *frame_mc = mass_mc.frame(Title(" "));
                data_mc.plotOn(frame_mc, Name("myData_mc"));
                model_mc.plotOn(frame_mc, Name("myModel_mc"));
                model_mc.plotOn(frame_mc, Components(background_mc), LineStyle(kDashed), LineColor(kRed));
                model_mc.plotOn(frame_mc, Components(signal_mc), LineStyle(kDashed), LineColor(kBlue));


                if(k == 10){

                
                    frame_mc->Draw();
                    TLatex latex_mc;
                    latex_mc.SetNDC();
                    latex_mc.SetTextSize(0.035);
                    latex_mc.DrawLatex(0.65, 0.85, Form("P_{Kp} #in [%.2f, %.2f] GeV^{2}", Kp_lim1, Kp_lim2));
                    latex_mc.DrawLatex(0.65, 0.80,Form("P_{Km} #in [%.2f, %.2f] GeV^{2}", Km_lim1, Km_lim2));

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

                    TLatex latex4_mc;
                    latex4_mc.SetNDC();
                    latex4_mc.SetTextSize(0.035);
                    latex4_mc.SetTextColor(kGreen);
                    latex4_mc.DrawLatex(0.65, 0.65, Form("#sigma MC = %.6f", sigma_mc.getVal()));

                    if(isplot){
                        cfit->Print("Fits_Minv_Smearing.pdf");
                    }

                }

                stocker_Nfit->Fill(sigma_mc.getVal()); // fill l'histo des Nfit
                //cout << "Sigma MC : " << sigma_mc <<  " et sigma DATA " << sigma << endl;

                delete frame_mc;
                delete hSignal;
                delete hBkg;
                delete hTot;
                delete res_mc;

                }

                f_Nfit->cd();
                stocker_Nfit->Write();


                liste_de_listes_sigmaMC[i].push_back(stocker_Nfit->GetMean());
                delete stocker_Nfit;

        }
    }

    f_Nfit->Close();
    delete f_Nfit;

    f->Close();
    fMC->Close();

    delete f;
    delete fMC;
    //delete cfit;

    
    if(isplot){
        cfit->Print("Fits_Minv_Smearing.pdf]");
    }

    //cout << "le nombre total de phi est " << Nphi_tot << endl;

    return {liste_de_listes_sigmaDATA, liste_de_listes_sigmaMC};


}

#endif


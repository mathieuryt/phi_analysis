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


void acceptance2() {

    std::vector<double> liste = {2, 3.05333, 3.334, 3.566, 3.788, 4.006, 4.23333, 4.49, 4.76933, 5.1, 5.494, 6.076, 7.07933};
    double bin = 12;
    double bin2 = 200;

    double elim_a = 4;
    double elim_b = 5;

    double Q_in = 35.667;
    double Q_out = 32.451;
    double Q_s = 45.9;

    // Charger le fichier ROOT et l'histogramme
    gEnv->SetValue("Browser.Name", "TRootBrowser");
    TFile *file = TFile::Open("outputTCS_Phi_inbending.root");
    TTree *t_rec_in = (TTree*)file->Get("tree;1");
    TTree *t_gen_in = (TTree*)file->Get("tree_Gen;3");

    gEnv->SetValue("Browser.Name", "TRootBrowser");
    TFile *file2 = TFile::Open("outputTCS_Phi_outbending.root");
    TTree *t_rec_out = (TTree*)file2->Get("tree;1");
    TTree *t_gen_out = (TTree*)file2->Get("tree_Gen;3");

    //INBENDING

    TLorentzVector *electron_rec_in = nullptr;
    TLorentzVector *positron_rec_in = nullptr;
    TLorentzVector *proton_rec_in = nullptr;
 
    Float_t Q2_rec_in, MMassBeam_rec_in, Ephovar_rec_in, weightvar_rec_in;

    t_rec_in->SetBranchAddress("Electron", &electron_rec_in);
    t_rec_in->SetBranchAddress("Positron", &positron_rec_in);
    t_rec_in->SetBranchAddress("Proton", &proton_rec_in);
    t_rec_in->SetBranchAddress("Q2", &Q2_rec_in);
    t_rec_in->SetBranchAddress("MMassBeam", &MMassBeam_rec_in);
    t_rec_in->SetBranchAddress("Epho", &Ephovar_rec_in);
    t_rec_in->SetBranchAddress("weight", &weightvar_rec_in);


    // Initialisation des histogrammes

    std::vector<TH1F*> histograms;

    for (int i = 0; i < bin; i++) {
        std::cout << "Initialisation de l'histogramme numéro : " << i << std::endl;
        histograms.push_back(new TH1F(Form("Histogramme numero %d", i), "", bin2, 0.5, 1.5));
    }

    std::cout << "Taille du vecteur d'histogrammes : " << histograms.size() << std::endl;


    for (int i = 0; i < bin; i++) {

        elim_a = liste[i];
        elim_b = liste[i+1];

        std::cout << "INBENDING: Pour Epho compris entre " << elim_a << " et " << elim_b << std::endl;

        Long64_t nentries = t_rec_in->GetEntries();
 
        for (Long64_t j = 0; j < nentries; j++) {
            t_rec_in->GetEntry(j);

            if (Ephovar_rec_in > elim_a && Ephovar_rec_in < elim_b) {

                double M = (*electron_rec_in + *positron_rec_in).M();
                histograms[i]->Fill(M, weightvar_rec_in*Q_in);

            }
        }
    }

    //OUTBENDING

    TLorentzVector *electron_rec_out = nullptr;
    TLorentzVector *positron_rec_out = nullptr;
    TLorentzVector *proton_rec_out = nullptr;
 
    Float_t Q2_rec_out, MMassBeam_rec_out, Ephovar_rec_out, weightvar_rec_out;

    t_rec_out->SetBranchAddress("Electron", &electron_rec_out);
    t_rec_out->SetBranchAddress("Positron", &positron_rec_out);
    t_rec_out->SetBranchAddress("Proton", &proton_rec_out);
    t_rec_out->SetBranchAddress("Q2", &Q2_rec_out);
    t_rec_out->SetBranchAddress("MMassBeam", &MMassBeam_rec_out);
    t_rec_out->SetBranchAddress("Epho", &Ephovar_rec_out);
    t_rec_out->SetBranchAddress("weight", &weightvar_rec_out);

    for (int i = 0; i < bin; i++) {

        elim_a = liste[i];
        elim_b = liste[i+1];

        std::cout << "OUTBENDING : Pour Epho compris entre " << elim_a << " et " << elim_b << std::endl;

        Long64_t nentries = t_rec_out->GetEntries();
 
        for (Long64_t j = 0; j < nentries; j++) {
            t_rec_out->GetEntry(j);

            if (Ephovar_rec_out > elim_a && Ephovar_rec_out < elim_b) {

                double M = (*electron_rec_out + *positron_rec_out).M();
                histograms[i]->Fill(M, weightvar_rec_out*(Q_s+Q_out));

            }
        }
    }



    for (int k = 0; k < bin; k++) {

        // Créer un nom dynamique pour le TCanvas avec l'index k
        std::string canvasName = "c" + std::to_string(k);  
        TCanvas *c = new TCanvas(canvasName.c_str(), "Masse Invariante", 800, 600);
        
        histograms[k]->SetLineColor(kBlue);
        histograms[k]->GetXaxis()->SetTitle("Invariant mass [GeV]");
        histograms[k]->GetYaxis()->SetTitle("Number of events");
        histograms[k]->Draw("hist");

    }


    // Création et ouverture d’un fichier ROOT en écriture
    TFile *outputFile = new TFile("Histos_M2_REC.root", "RECREATE");

    for (int i = 0; i < bin; ++i) {

        if (histograms[i]) {
            histograms[i]->Write();  // Sauvegarde de chaque histogramme
        }
    }

    outputFile->Close();  // Ferme proprement le fichier
    std::cout << "Tous les histogrammes ont été sauvegardés dans Histos_M2Q.root" << std::endl;



}

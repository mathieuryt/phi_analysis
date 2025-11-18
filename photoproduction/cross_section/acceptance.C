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


void acceptance() {

    std::vector<double> liste = {2, 3.05333, 3.334, 3.566, 3.788, 4.006, 4.23333, 4.49, 4.76933, 5.1, 5.494, 6.076, 7.07933};
    double bin = 12;
    double bin2 = 200;

    double borne_plot_a = 0;
    double borne_plot_b = 11;

    const double m_elec = 0.000511;
    const double m_positron = 0.000511;

    double Q_in = 35.667;
    double Q_out = 32.451;
    double Q_s = 45.9;

    double Elim_histo_a = 1.58;
    double Elim_histo_b = 11;

    double echelle_plot_couleur = 0.00003;
    double echelle_plot_couleur2 = 0.00005;
    double hp_max = 10;
    double hangle_max = 70;

    // Charger le fichier ROOT et l'histogramme
    gEnv->SetValue("Browser.Name", "TRootBrowser");
    TFile *file = TFile::Open("outputTCS_Phi_inbending.root");
    TTree *t_rec_in = (TTree*)file->Get("tree;1");
    TTree *t_gen_in = (TTree*)file->Get("tree_Gen;3");

    gEnv->SetValue("Browser.Name", "TRootBrowser");
    TFile *file2 = TFile::Open("outputTCS_Phi_outbending.root");
    TTree *t_rec_out = (TTree*)file2->Get("tree;1");
    TTree *t_gen_out = (TTree*)file2->Get("tree_Gen;3");

    TH1F *Epho_gen = new TH1F("Epho_gen", "Histo nb ev gen par bin de Egamma", bin, liste.data());
    Epho_gen->GetXaxis()->SetRangeUser(borne_plot_a, borne_plot_b);
    Epho_gen->GetYaxis()->SetRangeUser(0, 290);
    Epho_gen->GetXaxis()->SetTitle("E_gamma [GeV]");
    Epho_gen->GetYaxis()->SetTitle("N_GEN");

 
    Float_t Ephovar_gen_in, weightvar_gen_in;

    t_gen_in->SetBranchAddress("Epho_Gen", &Ephovar_gen_in);
    t_gen_in->SetBranchAddress("weight", &weightvar_gen_in);

    Long64_t nentries_gen_in = t_gen_in->GetEntries();

    std::cout << "le nombre d'entries GEN est = " << nentries_gen_in << std::endl;

    for (Long64_t i = 0; i < nentries_gen_in; i++) {

        t_gen_in->GetEntry(i);

        Epho_gen->Fill(Ephovar_gen_in, weightvar_gen_in*Q_in);
 
    }

    Float_t Ephovar_gen_out, weightvar_gen_out;

    t_gen_out->SetBranchAddress("Epho_Gen", &Ephovar_gen_out);
    t_gen_out->SetBranchAddress("weight", &weightvar_gen_out);

    Long64_t nentries_gen_out = t_gen_out->GetEntries();

    std::cout << "le nombre d'entries GEN est = " << nentries_gen_out << std::endl;

    for (Long64_t i = 0; i < nentries_gen_out; i++) {

        t_gen_out->GetEntry(i);

        Epho_gen->Fill(Ephovar_gen_out, weightvar_gen_out*(Q_out+Q_s));
 
    }


    TCanvas *c2 = new TCanvas("c2", "nb ev GEN par bin Egamma" , 800, 600);
    Epho_gen->Draw("hist");  

    std::vector<float> liste_phigen;



    for (int i = 1; i <= bin; ++i) { 
        double phigen = Epho_gen->GetBinContent(i);
        liste_phigen.push_back(phigen);
    }



    std::ofstream fichier2("fichier_phigen.txt");

    if (fichier2) {
        for (float valeur : liste_phigen) {
            fichier2 << valeur << std::endl;
        }
        fichier2.close();
        std::cout << "Fichier des phi GEN écrit " << std::endl;
    }


}

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include <fstream>
#include <iostream>
#include <string>
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
#include <fstream>
#include <cmath>
#include <iomanip>
#include <algorithm>
using namespace std;


double Lambda(double x, double y, double z) {
    return (x - y - z)*(x - y - z) - 4 * y*z;
}

//From Byukling Kayanti Formula (5.14) Page 86

double T_min(double E) // arguments are squares of masses of particles in the reaction a+b->1+2, and s is the square of the total c.m. energy i.e. (a+b)^2
{
    const double Mp = 0.938;
    const double M_phi = 1.019;
    double ma_2 = 0;
    double mb_2 = Mp*Mp;

    double m1_2 = 1.019*1.019;
    double m2_2 = Mp*Mp;
    double s = Mp*Mp + 2*Mp*E;
    return abs(ma_2 + m1_2 - (1 / (2 * s))*((s + ma_2 - mb_2)*(s + m1_2 - m2_2) - sqrt(Lambda(s, ma_2, mb_2) * Lambda(s, m1_2, m2_2))));
}

//From Byukling Kayanti Formula (5.14) page 86

double T_max(double E) {

    const double Mp = 0.938;
    const double M_phi = 1.019;

    double ma_2 = 0;
    double mb_2 = Mp*Mp;

    double m1_2 = 1.019*1.019;
    double m2_2 = Mp*Mp;
    double s = Mp*Mp + 2*Mp*E;

    return abs(ma_2 + m1_2 - (1 / (2 * s))*((s + ma_2 - mb_2)*(s + m1_2 - m2_2) + sqrt(Lambda(s, ma_2, mb_2) * Lambda(s, m1_2, m2_2))));
}

double T_min_wrap(const double *x, const double *p) {
    return T_min(x[0]);
}

double T_max_wrap(const double *x, const double *p) {
    return T_max(x[0]);
}


double Evstheta_elastique(double theta_elec) {
    double Eb = 10.6;
    double Mp = 0.938;
    return Eb/(1+(2*Eb/Mp)*sin((theta_elec*3.14)/(2*180))*sin((theta_elec*3.14)/(2*180)));
}

double Evstheta_elastique_wrap(double *x, double *p) {
    return Evstheta_elastique(x[0]);
}



// histo data |p| vs theta et |t| vs E_gamma

void histo_data()
{

   // Créer les histogrammes
   double bins = 300;
   TH2F *ptheta_proton = new TH2F("ptheta_proton", "|p| vs theta for proton spring2019 inbending", bins, 0, 80, bins, -0.1, 1); 
   TH2F *ptheta_electron = new TH2F("ptheta_electron", "|p| vs theta for electron spring2019 inbending", bins, 0, 40, bins, -0.1, 11); 
   TH2F *ptheta_positron = new TH2F("ptheta_positron", "|p| vs theta for positron spring2019 inbending", bins, 0, 40, bins, -0.1, 11);
   TH2F *tvsE = new TH2F("tvs_E", "|t| vs E_gamma spring2019 inbending", bins, 0 , 11, bins, 0, 11); 

   gEnv->SetValue("Browser.Name", "TRootBrowser");

   TChain *t1 = new TChain("tree");  // "tree" doit être le nom EXACT du TTree dans les fichiers

   // Ajoute tous tes fichiers ROOT
   t1->Add("Data_pass2_fall2018_outbending.root");
   t1->Add("Data_pass2_spring2019_inbending.root");
   t1->Add("Data_pass2_fall2018_inbending.root");

   TLorentzVector *electron = nullptr;
   TLorentzVector *positron = nullptr;
   TLorentzVector *proton = nullptr;

   Float_t Q2, MMassBeam, Ephovar;

   const double m_elec = 0.000511;
   const double m_positron = 0.000511;
   const double Mp = 0.938;
   const double M_phi = 1.019;


   t1->SetBranchAddress("Proton", &proton);
   t1->SetBranchAddress("Electron", &electron);
   t1->SetBranchAddress("Positron", &positron);
   t1->SetBranchAddress("Q2", &Q2);
   t1->SetBranchAddress("MMassBeam", &MMassBeam);
   t1->SetBranchAddress("Epho", &Ephovar);

   Long64_t nentries = t1->GetEntries();
   for (Long64_t i = 0; i < nentries; i++) {
    
       t1->GetEntry(i);

       if (!electron || !positron) continue;

       double M = (*electron + *positron).M();

        // calcul du quadrivecteur k = p -p'

        Double_t E = Mp - proton->E();
        Double_t px = 0 - proton->Px();
        Double_t py = 0 - proton->Py();
        Double_t pz = 0 - proton->Pz();

        // calcul de |t| = abs(k²)
        Double_t t = abs(E*E - px*px - py*py - pz*pz);

       if (M<1.10 && M>0.90 && abs(MMassBeam) < 0.4 && Q2 < 0.4){


            ptheta_proton->Fill(proton->Theta()*180 / 3.14, proton->P());
            ptheta_electron->Fill(electron->Theta()*180 / 3.14, electron->P());
            ptheta_positron->Fill(positron->Theta()*180 / 3.14, positron->P());

            // calcul du quadrivecteur k = p -p'

            Double_t E = Mp - proton->E();
            Double_t px = 0 - proton->Px();
            Double_t py = 0 - proton->Py();
            Double_t pz = 0 - proton->Pz();

            // calcul de |t| = abs(k²)
            Double_t t = abs(E*E - px*px - py*py - pz*pz);

            tvsE->Fill(Ephovar, t);
        }

   }





    TF1 *fTmin = new TF1("fTmin", T_min_wrap, 1.528, 11, 0);
    fTmin->SetLineColor(kRed);
    fTmin->SetLineWidth(2);

    TF1 *fTmax = new TF1("fTmax", T_max_wrap, 1.528, 11, 0);
    fTmax->SetLineColor(kBlue);
    fTmax->SetLineWidth(2);

    TF1 *Evstheta_elastique_fct = new TF1("Evstheta_elastique_fct", Evstheta_elastique_wrap, 2, 15);

    Evstheta_elastique_fct->SetLineWidth(2);





   std::cout << "Histo 2D fini" << std::endl;


   ptheta_proton->GetXaxis()->SetTitle("#theta [degree]");
   ptheta_proton->GetYaxis()->SetTitle("|p| [GeV]");

   ptheta_electron->GetXaxis()->SetTitle("#theta [degree]");
   ptheta_electron->GetYaxis()->SetTitle("|p| [GeV]");

   ptheta_positron->GetXaxis()->SetTitle("#theta [degree]");
   ptheta_positron->GetYaxis()->SetTitle("|p| [GeV]");

   tvsE->GetXaxis()->SetTitle("E_gamma [GeV]");
   tvsE->GetYaxis()->SetTitle("|t| [GeV^{2}]");

   double echelle_plot_couleur = 5;
   double echelle_plot_couleur2 = 80;

   TCanvas *c1 = new TCanvas("c1", "|p| vs #theta for proton fall2018_outbending", 800, 600); 

   ptheta_proton->SetMinimum(0);  
   ptheta_proton->SetMaximum(echelle_plot_couleur);
   ptheta_proton->Draw("COLZ");
   c1->Update();

   TCanvas *c2 = new TCanvas("c2", "|p| vs #theta for electron fall2018_outbending", 800, 600); 

   ptheta_electron->SetMinimum(0);  
   ptheta_electron->SetMaximum(echelle_plot_couleur);
   ptheta_electron->Draw("COLZ");
   Evstheta_elastique_fct->Draw("SAME");
   c2->Update();



   TCanvas *c3 = new TCanvas("c3", "|p| vs #theta for positron fall2018_outbending", 800, 600);

   ptheta_positron->SetMinimum(0);  
   ptheta_positron->SetMaximum(echelle_plot_couleur);
   ptheta_positron->Draw("COLZ");
   c3->Update();

   TCanvas *c4 = new TCanvas("c4", "|t| vs E_gamma", 800, 600);

   tvsE->SetMinimum(0);  
   tvsE->SetMaximum(echelle_plot_couleur2);
   tvsE->Draw("COLZ");

   fTmax->Draw("SAME");
   fTmin->Draw("SAME");

   // Charger les bornes verticales (axe des X : E_gamma)
   std::vector<double> bornes_x;
   std::ifstream fichier_x("/local/home/mr282803/Documents/irfu/diff_cs/diff_bins/fichier_bornes_bins.txt");
   double val;
   while (fichier_x >> val) {
       bornes_x.push_back(val);
   }

   // Tracer les droites verticales
   for (double x : bornes_x) {
       TLine *line_x = new TLine(x, 0, x, 11); // t va de 0 à 11
       line_x->SetLineColor(kGreen);
       line_x->SetLineStyle(2); // dashed
       line_x->SetLineWidth(6); // ou 3 ou plus
       line_x->Draw("SAME");
    }

   // Tracer les droites horizontales entre chaque paire de bornes_x
    for (size_t i = 0; i + 1 < bornes_x.size(); ++i) {
        std::string nom_fichier = "/local/home/mr282803/Documents/irfu/diff_cs/diff_bins/fichier_bornes_bins_" + std::to_string(i) + ".txt";
        std::ifstream fichier_t(nom_fichier);
        std::vector<double> bornes_y;

        while (fichier_t >> val) {
            bornes_y.push_back(val);
        }

        for (double y : bornes_y) {
            TLine *line_y = new TLine(bornes_x[i], y, bornes_x[i + 1], y);
            line_y->SetLineColor(kBlack); // couleur plus claire pour ne pas masquer les données
            line_y->SetLineStyle(2); // pointillés
            line_y->SetLineWidth(5); // ou 3 ou plus
            line_y->Draw("SAME");
        }
    }

   c4->Update();
}


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


void FluxQ2() {

    
    gROOT->Reset();
    gSystem->Load("libPhysics.so");
    gSystem->Load("libCore.so");
    gSystem->Load("libTree.so");
    gSystem->Load("libRIO.so");
    gSystem->Load("libHist.so");
    gStyle->SetPalette(kBird); 


    //TRAITEMENT BORNES EN Q2
    std::vector<double> binQ2moy;
    std::vector<double> binxbmoy; //per bin of Q2
    std::vector<double> FluxQ2;
    double valeur;

    std::string name_fichierQ2 = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_analysis/bornes_bins/bornes_Q2_moy.txt";
    std::ifstream fichierQ2(name_fichierQ2);

    while (fichierQ2 >> valeur) {
         binQ2moy.push_back(valeur);
     }
 
    fichierQ2.close();



    std::string name_fichierxb = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_analysis/bornes_bins/bornes_xb_moy.txt";
    std::ifstream fichierxb(name_fichierxb);

    while (fichierxb >> valeur) {
         binxbmoy.push_back(valeur);
     }
 
    fichierxb.close();

    for (size_t i = 0; i < binQ2moy.size(); ++i) {

        //double W = 2.5; // valeur moyenne pour l'instant 
        double Q2 = binQ2moy[i];
        double xb = binxbmoy[i];
        double E = 10.4; // ATTENTION a changer en fonction du data set

        double alpha = 1/137.036;
        double Mp = 0.9395;
        double nu = Q2/(2*Mp*xb);
        //double nu = (W*W - Mp*Mp + Q2)/(2*Mp);
        double y = nu/E;
        //double xb = Q2/(2*Mp*nu);
        double epsilon = (1 - y - (Q2/(4*E*E)) )/(1- y + (y*y)/2 + (Q2/(4*E*E)) );
        double tau = (alpha * Q2 * (1-xb))/(8*3.14*Mp*Mp*E*E*xb*xb*xb*(1-epsilon));


        FluxQ2.push_back(tau);

    }




    std::string chemin = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_crosssection/F/";

    
    // Nom du fichier : integral_t_0.txt, integral_t_1.txt, ...
    std::string filename = chemin + Form("FluxQ2.txt");
    std::ofstream outfile(filename);

    if (!outfile) {
        std::cerr << "Erreur ouverture fichier " << filename << std::endl;

    }

    for (size_t i = 0; i < binQ2moy.size(); ++i) {
    // Écrire tous les nsig pour ce Q² en une seule ligne

    outfile << FluxQ2[i] << std::endl;

    }

    outfile.close();
    std::cout << "File flux per bins of Q2 write : " << filename << std::endl;



    }
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
#include "TLine.h"
#include "CreateRoot.h"
#include "FillBinning.h"
#include "FitforSigma.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <memory>

double best_chi2 = 1e99;
std::vector<std::vector<double>> bestSigmaDATA;
std::vector<std::vector<double>> bestSigmaMC;

double bestA1, bestB1, bestC1, bestD1;
double bestA2, bestB2, bestC2, bestD2;

bool changebest = false;
int nBestPlot = 0;

double ComputeChi2(double A1, double B1, double C1, double D1,
                   double A2, double B2, double C2, double D2, bool isplot)
{

    auto result = FitforSigma("/Users/mr282803/Documents/analysisRGB/neutron/neutron_analysis/acceptance/HistoMinvbinP_data.root", "/Users/mr282803/Documents/analysisRGB/neutron/neutron_analysis/acceptance/TreeKpKm_MC.root", A1, B1, C1, D1, A2, B2, C2, D2, 40, isplot);
    auto sigmaDATA = result.first;
    auto sigmaMC   = result.second;

    //for (size_t i = 0; i < sigmaDATA.size(); ++i) { // Vérification : afficher le contenu
        //std::cout << "liste des sigma DATA en Km pour le bin en Kp numero " << i << std::endl;
        //for (double val : sigmaDATA[i]) {
            //std::cout << val << " ";
        //}
        //std::cout << std::endl;
    //}

    //for (size_t i = 0; i < sigmaMC.size(); ++i) { // Vérification : afficher le contenu
        //std::cout << "liste des sigma MC en Km pour le bin en Kp numero " << i << std::endl;
        //for (double val : sigmaMC[i]) {
            //std::cout << val << " ";
        //}
        //std::cout << std::endl;
   // }

    //CALCUL DE CHI2
    double chi2_tot = 0.0;
    double err = 0.0001; 

    for (size_t i = 0; i < sigmaDATA.size(); ++i) { // Vérification : afficher le contenu

        for (int j = 0; j < sigmaDATA[i].size(); j++) {

            chi2_tot += ((sigmaDATA[i][j] - sigmaMC[i][j])*(sigmaDATA[i][j] - sigmaMC[i][j]))/(err*err);
            
        }
        //std::cout << std::endl;
    }
    cout << "chi2 " << chi2_tot << endl;

    return chi2_tot;
}


void FindSmearing_midnight()
{
    fillbinning("/Users/mr282803/Documents/analysisRGB/neutron/neutron_data2/fall2019_outbending_neutronKpKm_data_v2.root", false);
    fillbinning("/Users/mr282803/Documents/analysisRGB/neutron/neutron_data2/fall2019_outbending_neutronKpKm_mc_v6.root", true);

    auto chi2Functor = [](const double *p) {

        double A1 = p[0];
        double B1 = p[1];
        double C1 = p[2];
        double D1 = p[3];

        // si tu veux imposer Kp = Km comme dans ta grille :
        double A2 = p[4];
        double B2 = p[5];
        double C2 = p[6];
        double D2 = p[7];

        return ComputeChi2(A1, B1, C1, D1,
                           A2, B2, C2, D2, false);
    };

    ROOT::Math::Functor fcn(chi2Functor, 8);

    std::unique_ptr<ROOT::Math::Minimizer> min(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex")
    );

    min->SetMaxFunctionCalls(5000); //nb d'apelle de chi2
    min->SetMaxIterations(1000); 
    min->SetTolerance(1e-4); //tolerence sur chi2
    min->SetPrintLevel(1); // cb d"info on print

    min->SetFunction(fcn);
    min->SetLimitedVariable(0, "A1", 0.031, 0.001, 0.028, 0.034);
    min->SetLimitedVariable(1, "B1", 0.1021, 0.001, 0.098, 0.110);
    min->SetLimitedVariable(2, "C1", 0.085, 0.001, 0.080, 0.090);
    min->SetLimitedVariable(3, "D1", 0.125, 0.001, 0.115, 0.130);

    min->SetLimitedVariable(4, "A2", 0.031, 0.001, 0.028, 0.034);
    min->SetLimitedVariable(5, "B2", 0.079, 0.001, 0.075, 0.087);
    min->SetLimitedVariable(6, "C2", 0.035, 0.001, 0.030, 0.040);
    min->SetLimitedVariable(7, "D2", 0.085, 0.001, 0.080, 0.090);

    cout << "Repetability ..." << endl;

    std::cout << "Chi2 param a 0.35 : " << ComputeChi2(0.035,0.035,0.035,0.035,
                         0.035,0.035,0.035,0.035,false) << std::endl;

    std::cout << "Chi2 param a 0.35 2eme fois : " << ComputeChi2(0.035,0.035,0.035,0.035,
                         0.035,0.035,0.035,0.035,false) << std::endl;

    std::cout << "Chi2 param a 0.35 3eme fois : " << ComputeChi2(0.035,0.035,0.035,0.035,
                         0.035,0.035,0.035,0.035,false) << std::endl;
    
    bool ok = min->Minimize();

    const double *xs = min->X();
    const double *err = min->Errors();

    std::cout << "\n===== RESULTAT MINUIT =====\n";
    std::cout << "status = " << min->Status() << "\n";
    std::cout << "ok     = " << ok << "\n";
    std::cout << "chi2   = " << min->MinValue() << "\n";

    std::cout << "A1 = " << xs[0] << " +/- " << err[0] << "\n";
    std::cout << "B1 = " << xs[1] << " +/- " << err[1] << "\n";
    std::cout << "C1 = " << xs[2] << " +/- " << err[2] << "\n";
    std::cout << "D1 = " << xs[3] << " +/- " << err[3] << "\n";

    std::cout << "A2 = " << xs[4] << " +/- " << err[4] << "\n";
    std::cout << "B2 = " << xs[5] << " +/- " << err[5] << "\n";
    std::cout << "C2 = " << xs[6] << " +/- " << err[6] << "\n";
    std::cout << "D2 = " << xs[7] << " +/- " << err[7] << "\n";

    ComputeChi2(xs[0], xs[1], xs[2], xs[3],
            xs[4], xs[5], xs[6], xs[7],
            true);


}
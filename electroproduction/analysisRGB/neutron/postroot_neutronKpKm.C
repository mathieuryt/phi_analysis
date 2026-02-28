#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TCut.h"
#include "TStyle.h"
#include "TLine.h"
#include "TString.h"
#include <vector>
#include <iostream>
#include "TLatex.h"

struct Plot1D {
    TString name;      // nom de l'histo root
    TString globaltitle; //titre du plot
    TString title;     // titre axe X
    TString var;       // variable TTree
    TString additionalplot; //utilse pour plot une variable differente genre nucleon.P() et missing_nucleon.P();
    int     nBins;
    double  xmin;
    double  xmax;
    double  vline;     // <0 → pas de ligne
    double  vline2;     // <0 → pas de ligne
    TString additionalcut; //cut additionel
    bool logy; //mettre en log
    bool plotMC; // superposé le mc
    bool plotMCcontamination;
    bool plotDATA;
};

struct Plot2D {
    TString name;          // nom de l'histo root
    TString globaltitle;
    TString titleX;        // titre axe X
    TString titleY;        // titre axe Y
    TString varX;          // variable X
    TString varY;          // variable Y
    int     nBinsX;
    double  xmin;
    double  xmax;
    int     nBinsY;
    double  ymin;
    double  ymax;
    TString additionalcut;
    bool plotMC;           // true → MC, false → Data
};

int postroot_neutronKpKm()
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetPalette(kBird); 
    //gStyle->SetOptStat(0);


    // ===============================
    // FICHIERS
    // ===============================

    TString data_adress;
    TString mc_adress;
    TString mc_contamination_adress;
    double integralMC;
    double integralMC_contamination;
    double integralDATA;


    bool fall2018_inbending = true;
    bool fall2018_outbending = false;
    bool spring2019 = false;

    if (fall2018_inbending){

        data_adress = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_data2/fall2019_outbending_neutronKpKm_bis_angle.root";
        mc_adress = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_data2/fall2019_outbending_neutronKpKm_mc_v2.root";
        mc_contamination_adress = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_data2/fall2019_outbending_contamination_proton_KpKm_mc_v1.root";

    }

    if (fall2018_outbending){

        data_adress = "";
        mc_adress = "";
        mc_contamination_adress = "";

    }

    if (spring2019){

        data_adress = "";
        mc_adress = "";
        mc_contamination_adress = "";
     
    }

    TFile *fData = new TFile(data_adress);
    TFile *fMC   = new TFile(mc_adress);
    TFile *fMC_contamination = new TFile(mc_contamination_adress);

    TTree *tData = (TTree*)fData->Get("tree");
    TTree *tMC   = (TTree*)fMC->Get("tree");
    TTree *tMC_contamination = (TTree*)fMC_contamination->Get("tree");


    // ===============================
    // CUTS 
    // ===============================
    TCut cut = "Q2>1.0 && Electron.P()>2 && Neutron.Theta()*180./3.141592 > 4 && Neutron.P()>0.25 && angle_neutron_missnucl*180./3.14159 < 5";


    // cut angle full : && angle_neutron_missnucl*180./3.14159 < 5
    // cut phi et theta : Neutron.Theta()*180./3.141592 - Missing_nucleon.Theta()*180./3.141592 < 15 && Neutron.Theta()*180./3.141592 - Missing_nucleon.Theta()*180./3.141592 > -15 && Neutron.Phi()*180./3.141592 - Missing_nucleon.Phi()*180./3.141592 > -15 && Neutron.Phi()*180./3.141592 - Missing_nucleon.Phi()*180./3.141592 < 15 
        //"positron_SF>0.15 && electron_SF>0.15 && "
        //"positron_HTCC_ECAL_match==1 && electron_HTCC_ECAL_match==1 && "
        //"abs(MMassBeam)<0.4 && Q2>0 && abs(Q2)<0.5";

    TCut weightMC = "weight_mc_correction2";


    TString status_cut ="Status_Kp >= 2000 && Status_Kp <= 2999 && Status_Km >= 2000 && Status_Km <= 2999";
    TString MM_cut = "Missing.M2() < 0.5 && Missing.M2() > -0.5";

    // ===============================
    // LISTE DES VARIABLES (LABELS)
    // ===============================
    std::vector<Plot1D> plots = {

        {"Q2"," ; ",   "Q^{2} [GeV^{2}]",     "Q2", "Nothing",       80, 0.0, 10.0, -1, -1   , "", false, true, false, true},
        {"Q2_2"," ; ",   "Q^{2} [GeV^{2}]",     "Q2", "Nothing",       80, 0.0, 10.0, -1, -1   , "Missing.M2() < 3.5 && Missing.M2() > -3.5", false, true, false, true},
        {"Q2_3"," ; ",   "Q^{2} [GeV^{2}]",     "Q2", "Nothing",       80, 0.0, 10.0, -1, -1   , "Missing.M2() < 0.5 && Missing.M2() > -0.5", false, true,false, true},
        {"Q2_4"," ; ",   "Q^{2} [GeV^{2}]",     "Q2", "Nothing",       80, 0.0, 10.0, -1, -1   , "angle_neutron_missnucl*180./3.14159 < 15", false, true, false,true},
        {"Q2_5"," ; ",   "Q^{2} [GeV^{2}]",     "Q2", "Nothing",       80, 0.0, 10.0, -1, -1   , status_cut + "&&" + "angle_neutron_missnucl*180./3.14159 < 15 && Missing.M2() < 0.5 && Missing.M2() > -0.5", false, true,false, true},
        
        {"t"," t without cuts; ",   "t [GeV^{2}]",     "t", "t_missing_nucleon",      80, -8.0, 0.0, -1,-1   , "",false, true,false, true},
        {"t_2"," t with status cut ; ",   "t [GeV^{2}]",     "t", "t_missing_nucleon",      80, -8.0, 0.0, -1,-1   , status_cut ,false, true,false, true},
        {"t_3"," t with status cut and MM cut ; ",   "t [GeV^{2}]",     "t", "t_missing_nucleon",      80, -8.0, 0.0, -1,-1   , status_cut + "&&" + MM_cut,false, true, false,true},
        {"t_4"," t missing with status cut and MM cut ; ",   "t [GeV^{2}]",     "t_missing_nucleon", "Nothing",      80, -8.0, 0.0, -1,-1   , status_cut + "&&" + MM_cut,false, true, false,true},
        {"W"," W with status cut; ",  "W [GeV]",     "W", "Nothing",        80, 0.0, 6.0, -1,-1, status_cut,false, true,false, true},
        {"W2"," W^{2} with status cut; ",  "W^{2} [GeV^{2}]",     "W*W", "Nothing",        80, 0.0, 15.0, -1,-1, status_cut,false, true,false, true},

        {"Electron_momentum","Electron momemtum ; ",  "p [GeV]",     "Electron.P()",  "Nothing",       80, 0.0, 12.0, -1,-1, status_cut,false, true,false, true},
        {"Electron_theta","Electron theta ; ",  "#theta [deg]",     "Electron.Theta()*180./3.141592",  "Nothing",       80, 0.0, 50.0, -1,-1, status_cut,false, true,false, true},
        {"Electron_phi","Electron phi ; ",  "#phi [deg]",     "Electron.Phi()*180./3.141592",  "Nothing",       80, -200.0, 200.0, -1,-1, status_cut,false, true,false, true},

        {"Kp_momentum","Kp momemtum ; ",  "p [GeV]",     "Kp.P()", "Nothing",        80, 0.0, 12.0, -1,-1, status_cut,false, true, false,true},
        {"Kp_theta","Kp theta ; ",  "#theta [deg]",     "Kp.Theta()*180./3.141592",  "Nothing",       80, 0.0, 50.0, -1,-1, status_cut,false, true, false,true},
        {"Kp_phi","Kp phi ; ",  "#phi [deg]",     "Kp.Phi()*180./3.141592",  "Nothing",       80, -200.0, 200.0, -1,-1, status_cut,false, true, false,true},

        {"Km_momentum","Km momemtum ; ",  "p [GeV]",     "Km.P()", "Nothing",        80, 0.0, 12.0, -1,-1, status_cut,false, true, false,true},
        {"Km_theta","Km theta ; ",  "#theta [deg]",     "Km.Theta()*180./3.141592","Nothing",         80, 0.0, 50.0, -1,-1, status_cut,false, true,false, true},
        {"Km_phi","Km phi ; ",  "#phi [deg]",     "Km.Phi()*180./3.141592",  "Nothing",       80, -200.0, 200.0, -1,-1, status_cut,false, true,false, true},

        {"Neutron_momentum","Neutron momemtum ; ",  "p [GeV]",     "Neutron.P()",  "Missing_nucleon.P()",       80, 0.0, 12.0, -1,-1, status_cut,false, true, false, true},
        {"Neutron_momentum_2","Missing nucleon momemtum ; ",  "p [GeV]",     "Missing_nucleon.P()",  "Nothing",       80, 0.0, 12.0, -1,-1, status_cut,false, true, false, true},
        {"Neutron_momentum_3","Missing nucleon momemtum ; ",  "p [GeV]",     "Neutron_gen_associated.P()",  "Nothing",       80, 0.0, 12.0, -1,-1, status_cut,false, true, false, false},
        {"Neutron_theta","Neutron theta ; ",  "#theta [deg]",     "Neutron.Theta()*180./3.141592", "Missing_nucleon.Theta()*180./3.141592",       80, 0.0, 70.0, -1,-1, status_cut,false, true,false, true},
        {"Neutron_theta_2","Neutron theta with FC ; ",  "#theta [deg]",     "Neutron.Theta()*180./3.141592", "Nothing",       80, 0.0, 70.0, -1,-1, status_cut + " && Status_n >= 2000 && Status_n <= 2999 && ((Lv_PCAL_n > 9.0 && Lw_PCAL_n > 9.0 && Lu_PCAL_n < 404.0) || x_PCAL_n == 0) && ((Lv_ECIN_n > 9.0 && Lw_ECIN_n > 9.0 && Lu_ECIN_n < 404.0) || x_ECIN_n == 0) && ((Lv_ECOUT_n > 9.0 && Lw_ECOUT_n > 9.0 && Lu_ECOUT_n < 404.0) || x_ECOUT_n == 0)",false, true,false, false},
        {"Neutron_phi"," ; ",  "#phi [deg]",     "Neutron.Phi()*180./3.141592",  "Missing_nucleon.Phi()*180./3.141592",       80, -200.0, 200.0, -1,-1, status_cut,false, true, false,true},
        {"Neutron_phi_3"," ; ",  "#phi [deg]",     "Neutron.Phi()*180./3.141592",  "Nothing",       80, -200.0, 200.0, -1,-1, status_cut ,false, true, false,false},
        {"Neutron_phi_2"," ; ",  "#phi [deg]",     "Neutron.Phi()*180./3.141592",  "Nothing",       350, -200.0, 200.0, -1,-1, status_cut + " && Status_n >= 2000 && Status_n <= 2999 && ((Lv_PCAL_n > 9.0 && Lw_PCAL_n > 9.0 && Lu_PCAL_n < 404.0) || x_PCAL_n == 0) && ((Lv_ECIN_n > 9.0 && Lw_ECIN_n > 9.0 && Lu_ECIN_n < 404.0) || x_ECIN_n == 0) && ((Lv_ECOUT_n > 9.0 && Lw_ECOUT_n > 9.0 && Lu_ECOUT_n < 404.0) || x_ECOUT_n == 0)",false, true, false,false},
        {"Neutron_phi_5"," ; ",  "#phi [deg]",     "Neutron.Phi()*180./3.141592",  "Nothing",       80, -200.0, 200.0, -1,-1, status_cut + " && Status_n >= 2000 && Status_n <= 2999 && ((Lv_PCAL_n > 9.0 && Lw_PCAL_n > 9.0 && Lu_PCAL_n < 404.0) || x_PCAL_n == 0) && ((Lv_ECIN_n > 9.0 && Lw_ECIN_n > 9.0 && Lu_ECIN_n < 404.0) || x_ECIN_n == 0) && ((Lv_ECOUT_n > 9.0 && Lw_ECOUT_n > 9.0 && Lu_ECOUT_n < 404.0) || x_ECOUT_n == 0) && (x_PCAL_n != 0 || x_ECIN_n != 0 || x_ECOUT_n != 0) && e_vx < 0.1 && e_vx > -0.1 && e_vy < 0.1 && e_vy > -0.1 && n_vx < 1 && n_vx > -1 && n_vy < 1 && n_vy > -1",false, true, false,false},
        {"Neutron_phi_4"," ; ",  "#phi [deg]",     "Neutron.Phi()*180./3.141592",  "Nothing",       80, -200.0, 200.0, -1,-1, status_cut + " && Status_n >= 2000 && Status_n <= 2999 && x_PCAL_n >= 0 && x_ECIN_n >= 0 && x_ECOUT_n >= 0",false, true, false,false},
        {"Energy_cal"," sum E_pcal Ecin Ecout for x == 0 in all the 3 cal; ",  "Energy [GeV]",     "E_PCAL_n + E_ECIN_n + E_ECOUT_n",  "Nothing",       80, -20.0, 20.0, -1,-1, status_cut + "&& x_PCAL_n == 0 && x_ECIN_n == 0 && x_ECOUT_n == 0",false, true, false,false},
        {"Energy_cal"," sum E_pcal Ecin Ecout; ",  "Energy [GeV]",     "E_PCAL_n + E_ECIN_n + E_ECOUT_n",  "Nothing",       300, -1.0, 3.0, -1,-1, status_cut + " && E_PCAL_n + E_ECIN_n + E_ECOUT_n == 0" ,false, true, false,false},

        {"Delta_neutron_momentum","p_{neutron} - p_{missing_nucleon} ; ",  "p [GeV]",     "Neutron.P() - Missing_nucleon.P()",  "Nothing",       80, -10, 10.0, -1,-1, status_cut,false, true, false,true},
        {"Delta_neutron_theta","#theta_{neutron} - #theta_{missing_nucleon} ; ",  "#theta [deg]",     "Neutron.Theta()*180./3.141592 - Missing_nucleon.Theta()*180./3.141592", "Nothing",       80, -70, 70.0, -1,-1, status_cut,false, true, false,true},
        {"Delta_neutron_phi","#phi_{neutron} - #phi_{missing_nucleon} ; ",  "#phi [deg]",     "Neutron.Phi()*180./3.141592 - Missing_nucleon.Phi()*180./3.141592",  "Nothing",       80, -200.0, 200.0, -1,-1, status_cut,false, true,false, true},

        {"Delta_neutron_fullangle","angle between neutron and missing_nucleon ; ",  "angle [deg]",     "angle_neutron_missnucl*180./3.141592",  "Nothing",       80, -1.0, 80.0, -1 ,5.0, status_cut,false, true, true,true},


        //{"MinvPiPlusMoinsD","test meson rho 0 cut MM<0.1 ; ",   "M_{inv}(#pi^{-} #pi^{+})",     "MinvPiPlusMoinsD",         80, 0.0, 1.4, 0.775 , -1, "MissingMassD < 0.1", false, false, true},
        //{"MinvKsKlD", " test meson w cut sur MM pi0; ",   "M_{inv}(pi+ pi- pi0) [GeV^{2}]",     "MinvKsKlD",         80, 0.7, 1.5, 0.782, -1, "MissingMassD < 0.25 && MissingMassD > 0.08 ",false, false, true},
        //{"vertex_eks_cut"," vertex e- Ks after cuts on MM and Minv; ",   "||e^{-} - Ks||_{xyz} [cm]",     "((Ks_vx - e_vx)*(Ks_vx - e_vx) + (Ks_vy - e_vy)*(Ks_vy - e_vy) + (Ks_vz - e_vz)*(Ks_vz - e_vz))**(1/2)",         80, -0.4, 20.0, 1, 7, "MinvPiPlusMoinsD > 0.4 && MinvPiPlusMoinsD < 0.6 && MissingMassD > 0.4 && MissingMassD < 0.6",true, true, true},
        //{"vertex_pippim_cut"," vertex pi+ pi- after cuts on MM and Minv; ",   "||#pi^{-} - #pi^{+}||_{xyz} [cm]",     "((pip_vx - pim_vx)*(pip_vx - pim_vx) + (pip_vy - pim_vy)*(pip_vy - pim_vy) + (pip_vz - pim_vz)*(pip_vz - pim_vz))**(1/2)",         80, -0.2, 15.0, 0, 3.5, "MinvPiPlusMoinsD > 0.4 && MinvPiPlusMoinsD < 0.6 && MissingMassD > 0.4 && MissingMassD < 0.6", true, true, true},
        //{"vertex_eks_cut_2"," vertex e- Ks after cuts on MM and Minv; ",   "||e^{-} - Ks||_{xyz} [cm]",     "((Ks_vx - e_vx)*(Ks_vx - e_vx) + (Ks_vy - e_vy)*(Ks_vy - e_vy) + (Ks_vz - e_vz)*(Ks_vz - e_vz))**(1/2)",         80, -0.4, 20.0, 1, 7, "MinvPiPlusMoinsD > 0.4 && MinvPiPlusMoinsD < 0.6 && MissingMassD > 0.4 && MissingMassD < 0.6",false, true, false},
        //{"vertex_pippim_cut_2"," vertex pi+ pi- after cuts on MM and Minv; ",   "||#pi^{-} - #pi^{+}||_{xyz} [cm]",     "((pip_vx - pim_vx)*(pip_vx - pim_vx) + (pip_vy - pim_vy)*(pip_vy - pim_vy) + (pip_vz - pim_vz)*(pip_vz - pim_vz))**(1/2)",         80, -0.2, 15.0, 0, 3.5, "MinvPiPlusMoinsD > 0.4 && MinvPiPlusMoinsD < 0.6 && MissingMassD > 0.4 && MissingMassD < 0.6", false, true, false},
        //{"MM", "Missing Mass",   "MissingMass [GeV]",     "MM",         80, 0.0, 1.5, -2,-1, "", false, true, false},
        //{"MM", "MM ;",   "MM[GeV]",  "MissingMass",         80, -1.7, 1.7, 1.019,-1, "", false, true, true},


        {"MM_2", "MM_2 ;",   "MissingMass^{2}[GeV^{2}]",  "Missing.M2()",   "Nothing",      160, -2, 2, 0.5, 0.5, status_cut, false, true,false, true},
        {"MinvKpKm_hist_0", "MinvKpKm without cut ;",   "Minv (K+ K-) [GeV]",  "MinvKpKm",   "Nothing",     100, 0.95, 1.12, 1.019,-1, "", false, true, true,true},
        {"MinvKpKm_hist_1", "MinvKpKm with cut on status ;",   "Minv (K+ K-) [GeV]",  "MinvKpKm",    "Nothing",     100, 0.95, 1.12, 1.019,-1, status_cut, false, true, true,true},
        {"MinvKpKm_hist_2", "Minv(K+ K-) with cut on status and MM^{2};",   "Minv (K+ K-) [GeV]",  "MinvKpKm",  "Nothing",     100, 0.95, 1.12, 1.019,-1, status_cut + "&&" + MM_cut, false, true,true, true},
        {"MinvKpKm_hist_2bis", "Minv(K+ K-) with cut on status and MM^{2} bin en t;",   "Minv (K+ K-) [GeV]",  "MinvKpKm",  "Nothing",     100, 0.95, 1.12, 1.019,-1, "t < -0.5 && t > -0.90 && " + status_cut + "&&" + MM_cut, false, true,true, true},
        {"MinvKpKm_hist_3", "MinvKpKm with cut on status, MM^{2}, and missE, missPt.. ;",   "Minv (K+ K-) [GeV]",  "MinvKpKm", "Nothing",        100, 0.95, 1.12, 1.019,-1, status_cut + "&&" + MM_cut + "&&" + "Missing.E() < 1.0 && Missing.Pt() < 2.0 && Kp.P() > 0.2 && Km.P() > 0.2", false, true,true, true}

        /*

        {"vertex1"," ; ",   "#Delta v_{z} (e - #pi^{+}) [cm]",     "e_vz - pip_vz",         80, -15.0, 15.0, -1,-1, "",false, false, true},
        {"vertex2"," ; ",   "#Delta v_{z} (e - #pi^{-}) [cm]",     "e_vz - pim_vz",         80, -15.0, 15.0, -1,-1, "",false, false, true},
        {"vertex3"," ; ",   "#Delta v_{z} (e - Ks) [cm]",     "e_vz - Ks_vz",         80, -15.0, 15.0, -1,-1, "",false, false, true},
        {"vertex4"," ; ",   "||e^{-} - Ks||_{xyz} [cm]",     "((Ks_vx - e_vx)*(Ks_vx - e_vx) + (Ks_vy - e_vy)*(Ks_vy - e_vy) + (Ks_vz - e_vz)*(Ks_vz - e_vz))**(1/2)",         80, -0.4, 20.0, 1.5, 5, "",false, false, true},
        {"vertex5"," ; ",  "#Delta v_{z} (#pi^{-} - #pi^{+}) [cm]",     "pim_vz - pip_vz",         80, -15.0, 15.0, -1, -1, "",false, false, true},
        {"vertex6"," ; ",   "||#pi^{-} - #pi^{+}||_{xyz} [cm]",     "((pip_vx - pim_vx)*(pip_vx - pim_vx) + (pip_vy - pim_vy)*(pip_vy - pim_vy) + (pip_vz - pim_vz)*(pip_vz - pim_vz))**(1/2)",         80, -0.2, 15.0, 0, 2.5, "", false, false, true},

        {"status_pim","Status #pi^{-} ; ",   "status",     "Status_pim",         80, -5000, 5000.0, -1,-1, "", false, false, true},
        {"status_pip","Status #pi^{+} ; ",   "status",     "Status_pip",         80, -5000, 5000.0, -1,-1, "", false, false, true},
        {"status_pr","Status proton; ",   "status",     "Status_pr",         80, -5000, 5000.0, -1,-1, "", false, false, true},
        {"status_el","Status electron; ",   "status",     "Status_el",         80, -5000, 5000.0, -1,-1, "", false, false, true},


        {"MinvPiPlusMoinsD"," ; ",   "M_{inv}(#pi^{-} #pi^{+})",     "MinvPiPlusMoinsD",         80, 0.0, 1.4, 0.4 , 0.6, "", false, false, true},
        {"MissingMassD"," ; ",   "MissingMass [GeV]",     "MissingMassD",         80, 0.0, 1.4, 0.4, 0.6, "", false, false, true},
        {"MinvKsKlD", " ; ",   "M_{inv}(Ks Kl) [GeV^{2}]",     "MinvKsKlD",         80, 0.7, 1.5, 1.019, -1, "",false, false, true},

        {"MinvPiPlusMoinsD_meeting", "Minv with cut on MM, status, R1, R2 ; ",   "M_{inv}(#pi^{-} #pi^{+})",     "MinvPiPlusMoinsD",         80, 0.0, 1.4, 0.496,-1, MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,false, false, true},
        {"MissingMassD_meeting", "MM with cut on Minv, status, R1, R2 ; ",  "MissingMass [GeV]",     "MissingMassD",         80, 0.0, 1.4, 0.496,-1, Minv_cut + status_cut + " && " + R1_cut + " && " + R2_cut,false, false, true},
        {"MinvKsKlD_meeting", "MinvKsKl with cut on Minv, MM, status, R1, R2 ; ",   "Minv Ks Kl [GeV^{2}]",     "MinvKsKlD",         80, 0.7, 1.5, 1.019,-1, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,false, false, true},

        {"MinvPiPlusMoinsD_mconly", "Minv MC with cut on MM, status, R1, R2 ; ",   "M_{inv}(#pi^{-} #pi^{+})",     "MinvPiPlusMoinsD",         80, 0.0, 1.4, 0.496,-1, MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,false, true, false},
        {"MissingMassD_mconly", "MM MC with cut on Minv, status, R1, R2 ; ",  "MissingMass [GeV]",     "MissingMassD",         80, 0.0, 1.4, 0.496,-1, Minv_cut + status_cut + " && " + R1_cut + " && " + R2_cut,false, true, false},
        {"MinvKsKlD_mconly", "MinvKsKl MC with cut on Minv, MM, status, R1, R2 ; ",   "Minv Ks Kl [GeV^{2}]",     "MinvKsKlD",         80, 0.7, 1.5, 1.019,-1, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,false, true, false},
        

        {"MinvPiPlusMoinsD_2", "Minv with cut on MM, status, R1, R2 ; ",   "M_{inv}(#pi^{-} #pi^{+})",     "MinvPiPlusMoinsD",         80, 0.0, 1.4, 0.496,-1, MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,false, true, true},
        {"MinvPiPlusMoinsD_3", "Minv with cut on MM, status, R1, R2 ; ",   "M_{inv}(#pi^{-} #pi^{+})",     "MinvPiPlusMoinsD",         80, 0.0, 1.4, 0.496,-1, MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,true, true, true},
        {"MissingMassD_2", "MM with cut on Minv, status, R1, R2 ; ",  "MissingMass [GeV]",     "MissingMassD",         80, 0.0, 1.4, 0.496,-1, Minv_cut + status_cut + " && " + R1_cut + " && " + R2_cut,false, true, true},
        {"MissingMassD_3", "MM with cut on Minv, status, R1, R2 ; ",  "MissingMass [GeV]",     "MissingMassD",         80, 0.0, 1.4, 0.496,-1, Minv_cut + status_cut + " && " + R1_cut + " && " + R2_cut,true, true, true},
        {"MinvKsKlD_2", "MinvKsKl with cut on Minv, MM, status, R1, R2 ; ",   "Minv Ks Kl [GeV^{2}]",     "MinvKsKlD",         80, 0.7, 1.5, 1.019,-1, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,false, true, true},
        {"MinvKsKlD_3", "MinvKsKl with cut on Minv, MM, status, R1, R2 ; ",   "Minv Ks Kl [GeV^{2}]",     "MinvKsKlD",         80, 0.7, 1.5, 1.019,-1, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,true, true, true},

        {"p_electron", "Electron.P() with cut on Minv, MM, M_{#phi} status, R1, R2 ; ",   "Electron.P() [GeV]",     "Electron.P()",         80, 0.0, 10.0, -1,-1, Minv_MM_cut + Mphi_cut + status_cut + " && " + R1_cut + " && " + R2_cut,true, true, true},
        {"p_proton", "Proton.P() with cut on Minv, MM, M_{#phi}, status, R1, R2 ; ",   "Proton.P() [GeV]",     "Proton.P()",         80, 0.0, 10.0, -1,-1, Minv_MM_cut + Mphi_cut + status_cut + " && " + R1_cut + " && " + R2_cut,true, true, true},
        {"p_pip", "PiPlus.P() with cut on Minv, MM, M_{#phi}, status, R1, R2 ; ",   "PiPlus.P() [GeV]",     "PiPlusD.P()",         80, 0.0, 8.0, -1,-1, Minv_MM_cut + Mphi_cut + status_cut + " && " + R1_cut + " && " + R2_cut,true, true, true},
        {"p_pim", "PiMinus.P() with cut on Minv, MM, M_{#phi}, status, R1, R2 ; ",   "PiMinus.P() [GeV]",     "PiMinusD.P()",         80, 0.0, 8.0, -1,-1, Minv_MM_cut + Mphi_cut + status_cut + " && " + R1_cut + " && " + R2_cut,true, true, true},

        {"theta_electron", "angle #theta for e- with cut on Minv, MM, M_{#phi}, status, R1, R2 ; ",   "#theta (e-) [GeV]",     "Electron.Theta()*180./3.141592",         80, 0.0, 45.0, -1,-1, Minv_MM_cut + Mphi_cut + status_cut + " && " + R1_cut + " && " + R2_cut, true, true, true},
        {"theta_proton", "angle #theta for p with cut on Minv, MM, M_{#phi}, status, R1, R2 ; ",   "#theta (p) [GeV]",     "Proton.Theta()*180./3.141592",         80, 0.0, 80.0, -1,-1, Minv_MM_cut + Mphi_cut + status_cut + " && " + R1_cut + " && " + R2_cut, true, true, true},
        {"theta_pip", "angle #theta for #pi^{+} with cut on Minv, MM, M_{#phi}, status, R1, R2 ; ",   "#theta (#pi^{+}) [GeV]",     "PiPlusD.Theta()*180./3.141592",         80, 0.0, 120.0, -1,-1, Minv_MM_cut + Mphi_cut + status_cut + " && " + R1_cut + " && " + R2_cut, true, true, true},
        {"theta_pim", "angle #theta for #pi^{-} with cut on Minv, MM, M_{#phi}, status, R1, R2 ; ",   "#theta (#pi^{-}) [GeV]",     "PiMinusD.Theta()*180./3.141592",         80, 0.0, 120.0, -1,-1, Minv_MM_cut + Mphi_cut + status_cut + " && " + R1_cut + " && " + R2_cut, true, true, true},

        */


        //{"Q2_2"," Q2 with cut on Minv, MM, status, R1, R2 ; ",   "Q^{2} [GeV^{2}]",     "Q2",         80, 0.0, 10.0, -1, -1   , Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut, true, true, true},
        //{"t_2"," t with cut on Minv, MM, status, R1, R2 ; ",   "t [GeV^{2}]",     "t",         80, -8.0, 0.0, -1,-1   , Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,true, true, true},
        //{"W_2"," W with cut on Minv, MM, status, R1, R2 ; ",  "W [GeV]",     "W",         80, 0.0, 6.0, -1,-1, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,true, true, true}
        
    };

    std::vector<Plot2D> plots2 = {

    { "t_vs_Q2_data", "t vs Q^{2}" , "t [GeV^{2}]" , "Q^{2} [GeV^{2}]" , "t" , "Q2" ,       60, -8.0, 0.0,     60, 0.0, 10.0, "",   false },
    { "t_vs_Q2_MC", "t vs Q^{2}" , "t [GeV^{2}]" , "Q^{2} [GeV^{2}]" , "t" , "Q2" ,       60, -8.0, 0.0,     60, 0.0, 10.0, "",   true },
    { "Q2_vs_xbj_data", "Q^{2} vs x_{bj}" , "x_{bj}" , "Q^{2} [GeV^{2}]" , "Q2/(2*0.939*(10.4-Electron.E()))" , "Q2" ,       60, 0.0, 1.0,     60, 0.0, 10.0, "Status_n >= 2000 && Status_n <= 2999 && Status_Kp >= 2000 && Status_Kp <= 2999 && Status_Km >= 2000 && Status_Km <= 2999",   false },

    //theta vs theta
    { "thetatheta_neutron_km_data", "#theta_{neutron} vs #theta_{K^{-}} with cut status" , "#theta_{neutron} [degree]" , "#theta_{K^{-}} [degree]" , "Neutron.Theta()*180./3.141592" , "Km.Theta()*180./3.141592" ,       60, 0.0, 120.0,     60, 0.0, 120.0, status_cut,   false },
    { "thetatheta_neutron_km_mc", "#theta_{neutron} vs #theta_{K^{-}} with cut status" , "#theta_{neutron} [degree]" , "#theta_{K^{-}} [degree]" , "Neutron.Theta()*180./3.141592" , "Km.Theta()*180./3.141592" ,       60, 0.0, 120.0,     60, 0.0, 120.0, status_cut,   true },

    { "thetatheta_neutron_kp_data", "#theta_{neutron} vs #theta_{K^{+}} with cut status" , "#theta_{neutron} [degree]" , "#theta_{K^{+}} [degree]" , "Neutron.Theta()*180./3.141592" , "Kp.Theta()*180./3.141592" ,       60, 0.0, 120.0,     60, 0.0, 120.0, status_cut,   false },
    { "thetatheta_neutron_kp_mc", "#theta_{neutron} vs #theta_{K^{+}} with cut status" , "#theta_{neutron} [degree]" , "#theta_{K^{+}} [degree]" , "Neutron.Theta()*180./3.141592" , "Kp.Theta()*180./3.141592" ,       60, 0.0, 120.0,     60, 0.0, 120.0, status_cut,   true },

    //theta vs p
    { "thetavsp_neutron_data", "#theta_{neutron} vs p_{neutron} with all cut" , "#theta_{neutron} [degree]" , "p_{neutron} [GeV]" , "Neutron.Theta()*180./3.141592" , "Neutron.P()" ,       60, 0.0, 120.0,     60, 0.0, 8.0, status_cut + " && " + MM_cut,   false },
    { "thetavsp_neutron_mc", "#theta_{neutron} vs p_{neutron} with cut status" , "#theta_{neutron} [degree]" , "p_{neutron} [GeV]" , "Neutron.Theta()*180./3.141592" , "Neutron.P()" ,       60, 0.0, 120.0,     60, 0.0, 8.0, status_cut,   true },

    { "thetavsp_kp_data", "#theta_{K^{+}} vs p_{K^{+}} with cut status" , "#theta_{K^{+}} [degree]" , "p_{K^{+}} [GeV]" , "Kp.Theta()*180./3.141592" , "Kp.P()" ,       60, 0.0, 120.0,     60, 0.0, 8.0, status_cut,   false },
    { "thetavsp_kp_mc", "#theta_{K^{+}} vs p_{K^{+}} with cut status" , "#theta_{K^{+}} [degree]" , "p_{K^{+}} [GeV]" , "Kp.Theta()*180./3.141592" , "Kp.P()" ,       60, 0.0, 120.0,     60, 0.0, 8.0, status_cut,   true },

    { "thetavsp_km_data", "#theta_{K^{-}} vs p_{K^{-}} with cut status" , "#theta_{K^{-}} [degree]" , "p_{K^{-}} [GeV]" , "Km.Theta()*180./3.141592" , "Km.P()" ,       60, 0.0, 120.0,     60, 0.0, 8.0, status_cut,   false },
    { "thetavsp_km_mc", "#theta_{K^{-}} vs p_{K^{-}} with cut status" , "#theta_{K^{-}} [degree]" , "p_{K^{-}} [GeV]" , "Km.Theta()*180./3.141592" , "Km.P()" ,       60, 0.0, 120.0,     60, 0.0, 8.0, status_cut,   true },

    { "xvsy_PCAL_n_mc", "x vs y neutron PCAL without FC" , "x [cm]" , "y [cm]" , "x_PCAL_n" , "y_PCAL_n" ,       200, -450.0, 450.0,     200, -450.0, 450.0, "Status_n >= 2000 && Status_n <= 2999 && x_PCAL_n != 0",   true },
    { "xvsy_PCAL_n_mc_2", "x vs y neutron PCAL do not pass FC" , "x [cm]" , "y [cm]" , "x_PCAL_n" , "y_PCAL_n" ,       200, -450.0, 450.0,     200, -450.0, 450.0, "(Lv_PCAL_n < 14.0 || Lw_PCAL_n < 14.0 || Lu_PCAL_n > 403.0) && Status_n >= 2000 && Status_n <= 2999 && x_PCAL_n != 0",   true },

    { "xvsy_ECIN_n_mc", "x vs y neutron ECIN without FC" , "x [cm]" , "y [cm]" , "x_ECIN_n" , "y_ECIN_n" ,       200, -450.0, 450.0,     200, -450.0, 450.0, "Status_n >= 2000 && Status_n <= 2999 && x_ECIN_n != 0",   true },
    { "xvsy_ECIN_n_mc_2", "x vs y neutron ECIN do not pass FC" , "x [cm]" , "y [cm]" , "x_ECIN_n" , "y_ECIN_n" ,       200, -450.0, 450.0,     200, -450.0, 450.0, "(Lv_ECIN_n < 14.0 || Lw_ECIN_n < 14.0 || Lu_ECIN_n > 403.0) && Status_n >= 2000 && Status_n <= 2999 && x_ECIN_n != 0",   true },

    { "xvsy_ECOUT_n_mc", "x vs y neutron ECOUT without FC" , "x [cm]" , "y [cm]" , "x_ECOUT_n" , "y_ECOUT_n" ,       200, -450.0, 450.0,     200, -450.0, 450.0, "Status_n >= 2000 && Status_n <= 2999 && x_ECOUT_n != 0",   true },
    { "xvsy_ECOUT_n_mc_2", "x vs y neutron ECOUT do not pass FC" , "x [cm]" , "y [cm]" , "x_ECOUT_n" , "y_ECOUT_n" ,       200, -450.0, 450.0,     200, -450.0, 450.0, "(Lv_ECOUT_n < 14.0 || Lw_ECOUT_n < 14.0 || Lu_ECOUT_n > 403.0) && Status_n >= 2000 && Status_n <= 2999 && x_ECOUT_n != 0",   true },

    { "thetavsphi_neutron", "theta vs phi neutron" , "phi [degree]" , "theta [degree]" , "Neutron.Phi()*180./3.141592" , "Neutron.Theta()*180./3.141592" ,       200, -200.0, 200.0,     200, -10.0, 45.0, "",   true },
    { "thetavsphi_neutron2", "theta vs phi neutro with fiducial cut" , "phi [degree]" , "theta [degree]" , "Neutron.Phi()*180./3.141592" , "Neutron.Theta()*180./3.141592" ,       200, -200.0, 200.0,     200, -10.0, 45.0, status_cut + " && Status_n >= 2000 && Status_n <= 2999 && ((Lv_PCAL_n > 9.0 && Lw_PCAL_n > 9.0 && Lu_PCAL_n < 404.0) || x_PCAL_n == 0) && ((Lv_ECIN_n > 9.0 && Lw_ECIN_n > 9.0 && Lu_ECIN_n < 404.0) || x_ECIN_n == 0) && ((Lv_ECOUT_n > 9.0 && Lw_ECOUT_n > 9.0 && Lu_ECOUT_n < 404.0) || x_ECOUT_n == 0)",   true },
    { "thetavsphi_neutron3", "theta vs phi neutron with fiducial cut and vertex vx vy cut" , "phi [degree]" , "theta [degree]" , "Neutron.Phi()*180./3.141592" , "Neutron.Theta()*180./3.141592" ,       100, -200.0, 200.0,     100, -10.0, 45.0, status_cut + " && Status_n >= 2000 && Status_n <= 2999 && ((Lv_PCAL_n > 9.0 && Lw_PCAL_n > 9.0 && Lu_PCAL_n < 404.0) || x_PCAL_n == 0) && ((Lv_ECIN_n > 9.0 && Lw_ECIN_n > 9.0 && Lu_ECIN_n < 404.0) || x_ECIN_n == 0) && ((Lv_ECOUT_n > 9.0 && Lw_ECOUT_n > 9.0 && Lu_ECOUT_n < 404.0) || x_ECOUT_n == 0) && e_vx < 1 && e_vx > -1 && e_vy < 1 && e_vy > -1 && n_vx < 1 && n_vx > -1 && n_vy < 1 && n_vy > -1",   true },





    /*
    //p vs theta 
    { "ptheta_e_data", " p vs #theta electron" , "#theta [degree]" , "p [GeV]" , "Electron.Theta()*180./3.141592" , "Electron.P()" ,       60, 0.0, 45.0,     60, 0.0, 10.0, "",   false },
    { "ptheta_e_MC", "p vs #theta electron" , "#theta [degree]" , "p [GeV]" , "Electron.Theta()*180./3.141592" , "Electron.P()" ,       60, 0.0, 45.0,     60, 0.0, 10.0, "",   true },

    { "ptheta_p_data", " p vs #theta proton" , "#theta [degree]" , "p [GeV]" , "Proton.Theta()*180./3.141592" , "Proton.P()" ,       60, 0.0, 70.0,     60, 0.0, 10.0, "",   false },
    { "ptheta_p_MC", "p vs #theta proton" , "#theta [degree]" , "p [GeV]" , "Proton.Theta()*180./3.141592" , "Proton.P()" ,       60, 0.0, 70.0,     60, 0.0, 10.0, "",   true },

    { "ptheta_pip_data", " p vs #theta #pi^{+} " , "#theta [degree]" , "p [GeV]" , "PiPlusD.Theta()*180./3.141592" , "PiPlusD.P()" ,       60, 0.0, 120.0,     60, 0.0, 10.0, "",   false },
    { "ptheta_pip_MC", "p vs #theta #pi^{+}" , "#theta [degree]" , "p [GeV]" , "PiPlusD.Theta()*180./3.141592" , "PiPlusD.P()" ,       60, 0.0, 120.0,     60, 0.0, 10.0, "",   true },

    { "ptheta_pim_data", " p vs #theta #pi^{-}" , "#theta [degree]" , "p [GeV]" , "PiMinusD.Theta()*180./3.141592" , "PiMinusD.P()" ,       60, 0.0, 120.0,     60, 0.0, 10.0, "",   false },
    { "ptheta_pim_MC", "p vs #theta #pi^{-}" , "#theta [degree]" , "p [GeV]" , "PiMinusD.Theta()*180./3.141592" , "PiMinusD.P()" ,       60, 0.0, 120.0,     60, 0.0, 10.0, "",   true },
    
    { "ptheta_e_data_2", " p vs #theta electron with cut on Minv, MM, status, R1, R2" , "#theta [degree]" , "p [GeV]" , "Electron.Theta()*180./3.141592" , "Electron.P()" ,       60, 0.0, 45.0,     60, 0.0, 10.0, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,   false },
    { "ptheta_e_MC_2", "p vs #theta electron with cut on Minv, MM, status, R1, R2" , "#theta [degree]" , "p [GeV]" , "Electron.Theta()*180./3.141592" , "Electron.P()" ,       60, 0.0, 45.0,     60, 0.0, 10.0, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,   true },

    { "ptheta_p_data_2", " p vs #theta proton with cut on Minv, MM, status, R1, R2" , "#theta [degree]" , "p [GeV]" , "Proton.Theta()*180./3.141592" , "Proton.P()" ,       60, 0.0, 70.0,     60, 0.0, 10.0, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,   false },
    { "ptheta_p_MC_2", "p vs #theta proton with cut on Minv, MM, status, R1, R2" , "#theta [degree]" , "p [GeV]" , "Proton.Theta()*180./3.141592" , "Proton.P()" ,       60, 0.0, 70.0,     60, 0.0, 10.0, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,   true },

    { "ptheta_pip_data_2", " p vs #theta #pi^{+} with cut on Minv, MM, status, R1, R2" , "#theta [degree]" , "p [GeV]" , "PiPlusD.Theta()*180./3.141592" , "PiPlusD.P()" ,       60, 0.0, 120.0,     60, 0.0, 10.0, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,   false },
    { "ptheta_pip_MC_2", "p vs #theta #pi^{+} with cut on Minv, MM, status, R1, R2" , "#theta [degree]" , "p [GeV]" , "PiPlusD.Theta()*180./3.141592" , "PiPlusD.P()" ,       60, 0.0, 120.0,     60, 0.0, 10.0, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,   true },

    { "ptheta_pim_data_2", " p vs #theta #pi^{-} with cut on Minv, MM, status, R1, R2" , "#theta [degree]" , "p [GeV]" , "PiMinusD.Theta()*180./3.141592" , "PiMinusD.P()" ,       60, 0.0, 120.0,     60, 0.0, 10.0, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,   false },
    { "ptheta_pim_MC_2", "p vs #theta #pi^{-} with cut on Minv, MM, status, R1, R2" , "#theta [degree]" , "p [GeV]" , "PiMinusD.Theta()*180./3.141592" , "PiMinusD.P()" ,       60, 0.0, 120.0,     60, 0.0, 10.0, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,   true },

    { "ptheta_Ks_data_2", " p vs #theta Ks with cut on Minv, MM, status, R1, R2" , "#theta [degree]" , "p [GeV]" , "Ks.Theta()*180./3.141592" , "Ks.P()" ,       60, 0.0, 120.0,     60, 0.0, 10.0, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,   false },
    { "ptheta_Ks_MC_2", "p vs #theta Ks with cut on Minv, MM, status, R1, R2" , "#theta [degree]" , "p [GeV]" , "Ks.Theta()*180./3.141592" , "Ks.P()" ,       60, 0.0, 120.0,     60, 0.0, 10.0, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,   true },

    { "ptheta_Kl_data_2", " p vs #theta Kl with cut on Minv, MM, status, R1, R2" , "#theta [degree]" , "p [GeV]" , "MissingD.Theta()*180./3.141592" , "MissingD.P()" ,       60, 0.0, 120.0,     60, 0.0, 10.0, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,   false },
    { "ptheta_Kl_MC_2", "p vs #theta Kl with cut on Minv, MM, status, R1, R2" , "#theta [degree]" , "p [GeV]" , "MissingD.Theta()*180./3.141592" , "MissingD.P()" ,       60, 0.0, 120.0,     60, 0.0, 10.0, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,   true },
    
    //theta vs theta
    { "thetatheta_prot_pim_MC", "#theta_{proton} vs #theta_{#pi^{-}} with cut on Minv, MM, status, R1, R2" , "#theta_{proton} [degree]" , "#theta_{#pi^{-}} [degree]" , "Proton.Theta()*180./3.141592" , "PiMinusD.Theta()*180./3.141592" ,       60, 0.0, 120.0,     60, 0.0, 120.0, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,   true },
    { "thetatheta_prot_pip_MC", "#theta_{proton} vs #theta_{#pi^{+}} with cut on Minv, MM, status, R1, R2" , "#theta_{proton} [degree]" , "#theta_{#pi^{+}} [degree]" , "Proton.Theta()*180./3.141592" , "PiPlusD.Theta()*180./3.141592" ,       60, 0.0, 120.0,     60, 0.0, 120.0, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,   true },

    //minv vs
    { "MinvKsKl_W_data", " M_{inv}(K_{s}K_{L}) vs W with cut on Minv, MM, status, R1, R2" , "W [GeV]" , "M_{inv}(K_{s}K_{L}) [GeV]" , "W" , "MinvKsKlD" ,       60, 1.0, 5.0,     60, 0.0, 3.5, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,   false },
    { "MinvKsKl_W_MC", " M_{inv}(K_{s}K_{L}) vs W with cut on Minv, MM, status, R1, R2" , "W [GeV]" , "M_{inv}(K_{s}K_{L}) [GeV]" , "W" , "MinvKsKlD" ,       60, 1.0, 5.0,     60, 0.0, 2, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,   true },

    { "MinvKsKl_t_data", " M_{inv}(K_{s}K_{L}) vs t with cut on Minv, MM, status, R1, R2" , "t [GeV^{2}]" , "M_{inv}(K_{s}K_{L}) [GeV]" , "t" , "MinvKsKlD" ,       60, -8.0, 5.0,     60, 0.0, 3.5, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,   false },
    { "MinvKsKl_t_MC", " M_{inv}(K_{s}K_{L}) vs t with cut on Minv, MM, status, R1, R2" , "t [GeV^{2}]" , "M_{inv}(K_{s}K_{L}) [GeV]" , "t" , "MinvKsKlD" ,       60, -8.0, 5.0,     60, 0.0, 2.0, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,   true },

    { "MinvKsKl_Q2_data", " M_{inv}(K_{s}K_{L}) vs Q^{2} with cut on Minv, MM, status, R1, R2" , "Q^{2} [GeV^{2}]" , "M_{inv}(K_{s}K_{L}) [GeV]" , "Q2" , "MinvKsKlD" ,       60, 0.0, 10.0,     60, 0.0, 3.5, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,   false },
    { "MinvKsKl_Q2_MC", " M_{inv}(K_{s}K_{L}) vs Q^{2} with cut on Minv, MM, status, R1, R2" , "Q^{2} [GeV^{2}]" , "M_{inv}(K_{s}K_{L}) [GeV]" , "Q2" , "MinvKsKlD" ,       60, 0.0, 8.0,     60, 0.0, 2.0, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,   true },

    { "MinvKsKl_pe_data", " M_{inv}(K_{s}K_{L}) vs Electron.P() with cut on Minv, MM, status, R1, R2" , "p [GeV]" , "M_{inv}(K_{s}K_{L}) [GeV]" , "Electron.P()" , "MinvKsKlD" ,       60, 1.0, 10.0,     60, 0.0, 3.5, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,   false },
    { "MinvKsKl_pe_MC", " M_{inv}(K_{s}K_{L}) vs Elelcron.P() with cut on Minv, MM, status, R1, R2" , "p [GeV]" , "M_{inv}(K_{s}K_{L}) [GeV]" , "Electron.P()" , "MinvKsKlD" ,       60, 1.0, 10.0,     60, 0.0, 2.0, Minv_MM_cut + status_cut + " && " + R1_cut + " && " + R2_cut,   true }
    
    */

    //{ "Minvpions_MM_data", " M_{inv}(#pi^{-} #pi^{+}) vs MissingMass" , "Missing Mass[GeV]" , "M_{inv}(#pi^{-} #pi^{+}) [GeV]" , "MissingMassD" , "MinvPiPlusMoinsD" ,       80, 0.0, 1.4,     80, 0.0, 1.4,"",   false }

    };

    // ===============================
    // CANVAS + PDF
    // ===============================
    TCanvas *c = new TCanvas("c", "", 900, 700);
    TString pdf = "Data_vs_MC_allVars_fall2019outbending_neutronKpKpm_bis2.pdf";

    c->Clear();
    c->SetFillColor(0);

    TLatex latex;
    latex.SetTextAlign(13);   // align left, top
    latex.SetTextSize(0.045);

    latex.DrawLatexNDC(0.15, 0.85, "Plot analysis DATA and MC KpKm");
    latex.DrawLatexNDC(0.15, 0.75, "DATA : inbending Fall 2018");
    latex.DrawLatexNDC(0.15, 0.65, "Number of data files : ");
    latex.DrawLatexNDC(0.15, 0.55, "MC simulation : N_{gen} = 2 000 000");

    c->Print(pdf + "(");

    // ===============================
    // BOUCLE SUR LES PLOTS
    // ===============================
    for (const auto &p : plots) {

        if (p.logy)
            c->SetLogy();
        else
            c->SetLogy(0);

        // histos
        TH1D *hData = new TH1D("hData_"+p.name, "", p.nBins, p.xmin, p.xmax);
        TH1D *hData_bis = new TH1D("hData_bis_"+p.name, "", p.nBins, p.xmin, p.xmax); //pour superposer des plots de data type missing nucleon par exemple
        TH1D *hMC   = new TH1D("hMC_"+p.name,   "", p.nBins, p.xmin, p.xmax);
        TH1D *hMC_contamination   = new TH1D("hMC_contamination_"+p.name,   "", p.nBins, p.xmin, p.xmax);

        hData->Sumw2();
        hData_bis->Sumw2();
        hMC->Sumw2();
        hMC_contamination->Sumw2();

        TCut additionalcut_applied = p.additionalcut.Data();

        // fill
        tData->Draw(p.var + ">>hData_"+p.name, cut * additionalcut_applied, "goff");
        tMC->Draw(p.var + ">>hMC_"+p.name, cut * additionalcut_applied * weightMC, "goff");

        if(p.additionalplot != "Nothing"){

            tData->Draw(p.additionalplot + ">>hData_bis_"+p.name, cut * additionalcut_applied, "goff");

        }

        if(p.plotMCcontamination == true){

            tMC_contamination->Draw(p.var + ">>hMC_contamination_"+p.name, cut * additionalcut_applied * weightMC, "goff");

        }

        // normalisation MC → Data
        if (hMC->Integral() > 0){
            hMC->Scale((18.3/2)*1.28);
            integralMC = hMC->Integral();
        }

        if (hMC_contamination->Integral() > 0){
            hMC_contamination->Scale(18.3/2);
            integralMC_contamination = hMC_contamination->Integral();
        }

        int binMin = hData->FindBin(1.01);
        int binMax = hData->FindBin(1.035);

        double integralDATA = hData->Integral(binMin, binMax);

        

        // style
        hData->SetMinimum(0.01);
        hData->SetMaximum(hData->GetMaximum()*1.2);
        hData->SetMarkerStyle(20);
        hData->SetMarkerSize(0.6);   // points plus petits
        hData->SetMarkerColor(kBlue+1);  // points bleus
        hData->SetLineColor(kBlue+1);    // barres d’erreur bleues

        hData_bis->SetMarkerColor(kGreen+1);  // points vert 
        hData_bis->SetLineColor(kGreen+1);


        hMC->SetFillColorAlpha(kRed+1, 0.35);
        hMC->SetMarkerStyle(0);        // pas de points
        hMC->SetLineColor(0);    // couleur des barres d’erreur

        hMC_contamination->SetFillColorAlpha(kMagenta+1, 0.35);
        hMC_contamination->SetMarkerStyle(0);        // pas de points
        hMC_contamination->SetLineColor(0);    // couleur des barres d’erreur

        hData->SetTitle(p.globaltitle + p.title + ";Events");

        // draw
        c->Clear();
        if (p.plotDATA){

            hData->Draw("E");

        }

        if (hMC->Integral() > 0 && p.plotMC) {

            if (p.plotDATA == false){

                hMC->Draw("HIST");

            }


            hMC->Draw("HIST SAME");

            TH1D* hMCerr = (TH1D*)hMC->Clone();
            hMCerr->SetFillStyle(0);       // pas de remplissage
            hMCerr->SetLineColor(kRed+1);  // barres d'erreur rouges
            hMCerr->SetMarkerStyle(0);
            hMCerr->Draw("E SAME");             
    
        }

        if(p.plotDATA){

            hData->Draw("E SAME");

        }

        if(p.additionalplot != "Nothing"){

            hData_bis->Draw("E SAME");
        
        }

        if(p.plotMCcontamination){

            hMC_contamination->Draw("HIST SAME");
        }
        

        // ligne verticale optionnelle
        if (p.vline > 0) {

            TLine *l = new TLine(-p.vline, 0, -p.vline, hData->GetMaximum());
            l->SetLineColor(kGreen+2);
            l->SetLineStyle(2);
            l->SetLineWidth(2);
            l->Draw("same");

        }

        if (p.vline2 > 0) {

            TLine *l = new TLine(p.vline2, 0, p.vline2, hData->GetMaximum());
            l->SetLineColor(kGreen+2);
            l->SetLineStyle(2);
            l->SetLineWidth(2);
            l->Draw("same");

        }

        // legend
        //TLegend *leg = new TLegend(0.60, 0.52, 0.84, 0.72);
        TLegend *leg = new TLegend(0.60, 0.72, 0.84, 0.92);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);
        TLegendEntry *eData = leg->AddEntry(hData, "Detected neutron", "lep");

        if (p.additionalplot != "Nothing") {

            leg->AddEntry(hData_bis, "Missing_nucleon", "lep");


        }
        eData->SetMarkerSize(0.2); 
        if (hMC->Integral() > 0 && p.plotMC){
            leg->AddEntry(hMC, "MC neutron target", "f");

            TString mcText;
            mcText.Form("Int. MC = %.1f", integralMC);
            leg->AddEntry((TObject*)nullptr, mcText, "");
        }

            if (hMC_contamination->Integral() > 0 && p.plotMCcontamination){
            leg->AddEntry(hMC_contamination, "MC proton target", "f");

            TString mcText_contamination;
            mcText_contamination.Form("Int. MC = %.1f", integralMC_contamination);
            leg->AddEntry((TObject*)nullptr, mcText_contamination, "");
        }

        TString mcText2;
        mcText2.Form("Int. DATA = %.1f", integralDATA);
        //leg->AddEntry((TObject*)nullptr, mcText2, "");

        leg->Draw();

        c->Print(pdf);
    }

    // ===============================
    // BOUCLE SUR LES PLOTS 2D
    // ===============================

    c->SetLogy(0);


    for (const auto &p : plots2) {


    TH2D *h2 = new TH2D("h2_" + p.name, "", p.nBinsX, p.xmin, p.xmax, p.nBinsY, p.ymin, p.ymax);

    TCut additionalcut_applied = p.additionalcut.Data();

    // Remplissage
    if (p.plotMC) {

        tMC->Draw(p.varY + ":" + p.varX + ">>h2_" + p.name, cut * additionalcut_applied * weightMC, "goff");


    } else {
        tData->Draw(p.varY + ":" + p.varX + ">>h2_" + p.name, cut * additionalcut_applied, "goff");
    }

    // Titres

    //h2->SetTitle(p.globaltitle);
    h2->GetXaxis()->SetTitle(p.titleX);
    h2->GetYaxis()->SetTitle(p.titleY);
   
    c->SetLeftMargin(0.14);
    c->SetRightMargin(0.16);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);

    // Draw
    c->Clear();
    h2->Draw("COLZ");

    // Label Data / MC
    TLatex lat;
    lat.SetTextSize(0.04);
    lat.SetNDC();

    TString label = p.plotMC ? "Monte Carlo : " : "Data : ";
    label += p.globaltitle;

    lat.DrawLatex(0.15, 0.95, label);



    c->Print(pdf);

    }   


    c->Clear();
    c->Print(pdf + ")");

    return 0;
}

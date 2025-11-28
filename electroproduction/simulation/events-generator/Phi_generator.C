#include <iostream>
#include <cmath>
#include <random>
#include "TRandom3.h"
#include "TLegend.h"
#include <vector>
#include <fstream>
#include "TLatex.h"
#include "TF2.h"
#include "TF1.h"
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include <Math/Integrator.h>
#include "TH3D.h"
#include "TH3.h"
#include <memory>
#include <thread>
#include <atomic>
#include <filesystem> // C++17
namespace fs = std::filesystem;

// In this Class : cross section (clas12 parametrisation from proposal) and useful functions like Lambda(), t_min() and t_max()

class Physics {


   public:

       double Lambda(double x, double y, double z) {

         return x*x + y*y + z*z - 2*x*y - 2*x*z - 2*y*z;

       }

       
       double sigmaT(double W, double Q2){

        // simga T (Q2, W). gamma* + p -> phi + p 
        double W_th = 1.957;
        double Mp = 0.938;
        double m_phi = 1.019;
        double alpha_1 = 400;
        double alpha_2 = 1.0;
        double alpha_3 = 0.32;
        double c_T = alpha_1*pow(1 - (W_th*W_th)/(W*W), alpha_2)*pow(W, alpha_3);
        double vt = 3.0;
        double sigmaT = c_T/pow(1 + Q2/(m_phi*m_phi), vt);

        return sigmaT;

       }

       double R(double Q2){

        // simga L = R*sigmaT (Q2, W). gamma* + p -> phi + p 

        double m_phi = 1.019;
        double cR = 0.4;
        double R = (cR*Q2)/(m_phi*m_phi);
        return R;

       }

       double t_min_calcul(double Q2, double W){

        double Mp = 0.938;
        double m_phi = 1.019;

        double p_cm = pow(Lambda(W*W, Mp*Mp, -Q2), 0.5)/(2*W);
        double p_cm2 = pow(Lambda(W*W, Mp*Mp, m_phi*m_phi), 0.5)/(2*W);
        double e_cm = sqrt(p_cm*p_cm + Mp*Mp);
        double e_cm2 = sqrt(p_cm2*p_cm2 + Mp*Mp);

        double t_min_calcul = 2*(Mp*Mp-e_cm*e_cm2 + p_cm*p_cm2);
        return t_min_calcul;

       }

       double t_max_calcul(double Q2, double W){

        double Mp = 0.938;
        double m_phi = 1.019;

        double p_cm = pow(Lambda(W*W, Mp*Mp, -Q2), 0.5)/(2*W);
        double p_cm2 = pow(Lambda(W*W, Mp*Mp, m_phi*m_phi), 0.5)/(2*W);
        double e_cm = sqrt(p_cm*p_cm + Mp*Mp);
        double e_cm2 = sqrt(p_cm2*p_cm2 + Mp*Mp);

        double t_max_calcul = 2*(Mp*Mp-e_cm*e_cm2 - p_cm*p_cm2);
        return t_max_calcul;

       }

       double dsumdt(double t, double Q2, double W, double E){

        // t-depence dsigmaL/dt and dsigmaT/dt 
        // Then sum with episolon wich depend of E -> we obtain dsigma tot/dt for this processus gamma* + p -> phi + p : 

        double mg = 1.0;
        double Mp = 0.938;
        double m_phi = 1.019;
   

        double t_min = t_min_calcul(Q2, W);

        double F_int = mg/(3*pow(mg-t_min, 3));
        double Ft = mg/pow(mg-t, 4);

        double dsigmaTdt = sigmaT(W, Q2)*(Ft/F_int);
        double sigmaL = sigmaT(W, Q2)*R(Q2);
        double dsigmaLdt = sigmaL*(Ft/F_int);

        double nu = (W*W - Mp*Mp + Q2)/(2*Mp);
        double y = nu/E;
        double epsilon = (1 - y - (Q2/(4*E*E)) )/(1- y + (y*y)/2 + (Q2/(4*E*E)) );

        double dsumdt = dsigmaTdt + epsilon*dsigmaLdt;

        return dsumdt;
        

       }

       double dsigmadt_tot(double t, double Q2, double W, double E){

        // cross section total dsigma / dQ^2 dt xb 
        // it's like : photon flux (QED) * gamma* + p -> phi + p (QCD) wich gives ep->e'p'phi

        double alpha = 1/137.036;
        double Mp = 0.938;
        double nu = (W*W - Mp*Mp + Q2)/(2*Mp);
        double y = nu/E;
        double xb = Q2/(2*Mp*nu);
        double epsilon = (1 - y - (Q2/(4*E*E)) )/(1- y + (y*y)/2 + (Q2/(4*E*E)) );
        double tau = (alpha * Q2 * (1-xb))/(8*3.14*Mp*Mp*E*E*xb*xb*xb*(1-epsilon));
        double dsigmadt_tot = tau*dsumdt(t, Q2, W, E);

        return dsigmadt_tot;

       }


   };


// Generator : in root .x Phi_generator(1000, 2000) for generate 1000 lund file with 2000 events
// This code generate also a pdf "phi_generator_plots.pdf" with histograms to check

void Phi_generator(double nb_fichier, double nb_event) { // arg1 : number of lund file you want 
                                                         // arg2 : number of event you want per Lund file

   // Initialisation of histograms for the pdf file                                               
   double bin1 = 200;
   double bin2 = 100;

   TH2F *pvstheta_el = new TH2F("pvstheta_el", "p vs theta for electron", bin1, 0, 90, bin1, 0.0, 11);
   TH2F *pvstheta_pr = new TH2F("pvstheta_pr", "p vs theta for proton", bin1, 0, 90, bin1, 0.0, 11);
   TH2F *pvstheta_kl = new TH2F("pvstheta_kl", "p vs theta for kl", bin1, 0, 90, bin1, 0.0, 11); 
   TH2F *pvstheta_ks = new TH2F("pvstheta_ks", "p vs theta for ks", bin1, 0, 90, bin1, 0.0, 11); 
   TH2F *pvstheta_pip = new TH2F("pvstheta_pip", "p vs theta for pi+", bin1, 0, 90, bin1, 0.0, 11); 
   TH2F *pvstheta_pim = new TH2F("pvstheta_pim", "p vs theta for pi-", bin1, 0, 90, bin1, 0.0, 11); 
   TH2F *tvsQ2 = new TH2F("tvsQ2", "t vs Q2 with kinematics cuts", 300, -20, 0, 300, 0, 10);
   TH2F *Q2vsW = new TH2F("Q2vsW", "Q2 vs W", 100, 0, 10, 100, 0, 10); 

   TH1F *Minv_pip_pim = new TH1F("Minv_pip_pim", "Invariant mass of pi+ pi-", bin2, 0, 1.0);
   TH1F *Minv_pip_pim_kl = new TH1F("Minv_pip_pim_kl", "Invariant mass of (pi+ pi- + Kl)", 400, 0.8, 1.2);
   TH1F *hist_vz = new TH1F("hist_vz", "vz of e-", bin2, -6, 1.0);
   TH1F *hist_Ks_vx = new TH1F("hist_Ks_vx", "vertex x of Ks", bin2, -3, 3.0);
   TH1F *hist_Ks_vy = new TH1F("hist_Ks_vy", "vertex y of Ks", bin2, -3, 3.0);
   TH1F *hist_Ks_vz = new TH1F("hist_Ks_vz", "vertex z of Ks", bin2, -4, 4.0);
   TH1F *hist_dist_vevks = new TH1F("hist_dist_vevks", "distance (x, y, z) between vertex e and vertex Ks", bin2, 2, 4);
   TH1F *histt = new TH1F("histt", "t", bin2, -8, 0);
   TH1F *hist_Q2 = new TH1F("hist_Q2", "Q2", bin2, 0, 8);
   TH1F *hist_xb = new TH1F("hist_xb", "xb", bin2, 0, 1);
   TH1F *hist_W = new TH1F("hist_W", "W", bin2, 1.5, 4.5);
   TH1F *hist_p_electron = new TH1F("hist_p_electron", "p electron", bin2, 0, 12);
   TH1F *hist_p_proton = new TH1F("hist_p_proton", "p proton", bin2, 0, 2);
   TH1F *hist_p_Ks = new TH1F("hist_p_Ks", "p Ks", bin2, 0, 4.2);
   TH1F *hist_p_Kl = new TH1F("hist_p_Kl", "p Kl", bin2, 0, 4.2);
   TH1F *hist_theta_electron = new TH1F("hist_theta_electron", "theta electron", bin2, 0, 1);
   TH1F *hist_theta_proton = new TH1F("hist_theta_proton", "theta proton", bin2, 0, 2);
   TH1F *hist_theta_Ks = new TH1F("hist_theta_Ks", "theta Ks", bin2, 0, 3);
   TH1F *hist_theta_Kl = new TH1F("hist_theta_Kl", "theta Kl", bin2, 0, 3);
   TH1F *hist_phi_electron = new TH1F("hist_phi_electron", "phi electron", bin2, -4, 4);
   TH1F *hist_phi_proton = new TH1F("hist_phi_proton", "phi proton", bin2, -4, 3.18);
   TH1F *hist_phi_Ks = new TH1F("hist_phi_Ks", "phi Ks", bin2, -4, 4);
   TH1F *hist_phi_Kl = new TH1F("hist_phi_Kl", "phi Kl", bin2, -4, 4);
   
   auto setTitles = [](TH1 *h, const char *xtitle, const char *ytitle) {
               h->GetXaxis()->SetTitle(xtitle);
               h->GetYaxis()->SetTitle(ytitle);
               h->GetXaxis()->CenterTitle();
               h->GetYaxis()->CenterTitle();
   };

   setTitles(pvstheta_el, "#theta [degree]", "p [GeV]");
   setTitles(pvstheta_pr, "#theta [degree]", "p [GeV]");
   setTitles(pvstheta_kl, "#theta [degree]", "p [GeV]");
   setTitles(pvstheta_ks, "#theta [degree]", "p [GeV]");
   setTitles(pvstheta_pip, "#theta [degree]", "p [GeV]");
   setTitles(pvstheta_pim, "#theta [degree]", "p [GeV]");
   setTitles(tvsQ2, "t [GeV^2]", "Q2 [GeV^2]");
   setTitles(Q2vsW, "Q2 [GeV^2]", "W [GeV^2]");

   Minv_pip_pim->GetXaxis()->SetTitle("Minv pi+ pi- [GeV]");
   Minv_pip_pim_kl->GetXaxis()->SetTitle("Minv pi+ pi- + kl [GeV]");
   hist_vz->GetXaxis()->SetTitle("vz e- [cm]");
   hist_Ks_vx->GetXaxis()->SetTitle("vx Ks [cm]");
   hist_Ks_vy->GetXaxis()->SetTitle("vy ks [cm]");
   hist_Ks_vz->GetXaxis()->SetTitle("vz Ks [cm]");
   hist_dist_vevks->GetXaxis()->SetTitle("distance vertex e- and vertex Ks [cm]");
   histt->GetXaxis()->SetTitle("t [GeV^2]");
   hist_Q2->GetXaxis()->SetTitle("Q2 [GeV^2]");
   hist_xb->GetXaxis()->SetTitle("xb");
   hist_W->GetXaxis()->SetTitle("W [GeV]");
   hist_p_electron->GetXaxis()->SetTitle("p [GeV]");
   hist_p_proton->GetXaxis()->SetTitle("p [GeV]");
   hist_p_Ks->GetXaxis()->SetTitle("p [GeV]");
   hist_p_Kl->GetXaxis()->SetTitle("p [GeV]");
   hist_theta_electron->GetXaxis()->SetTitle("theta [rad]");
   hist_theta_proton->GetXaxis()->SetTitle("theta [rad]");
   hist_theta_Ks->GetXaxis()->SetTitle("theta [rad]");
   hist_theta_Kl->GetXaxis()->SetTitle("theta [rad]");
   hist_phi_electron->GetXaxis()->SetTitle("theta [rad]");
   hist_phi_electron->GetXaxis()->SetTitle("theta [rad]");
   hist_phi_electron->GetXaxis()->SetTitle("theta [rad]");
   hist_phi_electron->GetXaxis()->SetTitle("theta [rad]");

   //initialisation of lund files
   std::vector<std::ofstream> files(nb_fichier);

   //some counter usefule
   double sum_weight_phasespace = 0; //sum of weight associated with the size of phasespace
   double cpt_nan = 0; //counter of NaN 

   //LOOP in LUND FILE
   for (int i = 0; i < nb_fichier; ++i) {
       std::string filename = "PhiGen_" + std::to_string(i) + ".txt";
       files[i].open(filename);

   //LOOP in event for a lundfile fixed
   for (Int_t n = 0; n < nb_event; n++) {

   double indice_evenement = n; // index of event

   //Cross section and some useful functions
   Physics phys;

   //Mass
   double m_electron = 0.00051;
   double Mp = 0.938;
   double m_kaons = 0.4937;
   double m_phi = 1.019;
   double m_pions = 0.1395;

   //Kinamitics varibale
   double Q2_min = 1;
   double Q2_max = 6.5;
   double Eb = 10.6;

   double s_23 = m_electron*m_electron + 2*Mp*Eb + Mp*Mp; // The mandelstam variable s for the process 2->3 (ep->e'p'phi) 

   double Q2, W, W2, W2_min, W2_max, nu, nu_min, nu_max, xb, xb_min, xb_max, theta_electron, phi_electron, E_electron, 
   t, t_min, t_min2, t_max2, t_max, E_proton, s_22, u_22, theta_proton_gamma, phi_proton_gamma, theta_ks, phi_ks, theta_pi, phi_pi;

   double E_sum, Px_sum, Py_sum, Pz_sum; // useful variable to verify the conservation of energy and momemtum at the end

   //Vertex variable
   double Ks_px, Ks_py, Ks_pz, norme, Ks_px_norme, Ks_py_norme, Ks_pz_norme;
   double vz0, vx, vy, vz, Ks_vx, Ks_vy, Ks_vz;
   double cT_Ks = 2.8;


   //Weight 
   double weight; // weight of phase space
   double weight_crosssection; // weight of cross section
   double BR_kaons = 0.34; //branching ratio
   double BR_pions = 0.692; // branching ratio
   double finalweight; 

   //Initial quadrivector
   TLorentzVector Target(0.0, 0.0, 0.0, Mp);  
   TLorentzVector Beam(0.0, 0.0, Eb, sqrt(Eb*Eb + m_electron*m_electron));

   TLorentzVector Proton, Electron, q, Phi;

   // random generator
   std::random_device rd;
   std::mt19937 gen(rd()); 

   //Q2 generation
   std::uniform_real_distribution<double> distQ2(Q2_min, Q2_max);
   Q2 = distQ2(gen);
   
   //W2 phase space
   //find the ref for the next formula in Byckling, E., and Kajantie, K. (1973b). Particle kinematics. Wiley-Interscience. (chapter 2->3 scattering)
   W2_min = pow((m_phi + Mp), 2);
   W2_max = s_23 + m_electron*m_electron - (1/(2*m_electron*m_electron))*(( s_23 + m_electron*m_electron - Mp*Mp )*(2*m_electron*m_electron + Q2) - sqrt(phys.Lambda(s_23, m_electron*m_electron, Mp*Mp))*sqrt(phys.Lambda(-Q2, m_electron*m_electron, m_electron*m_electron)));

   //nu phase space
   nu_min = (W2_min + Q2 - Mp*Mp)/(2*Mp);
   nu_max = (W2_max + Q2 - Mp*Mp)/(2*Mp);

   //xb phase space 
   xb_min = Q2/(2*Mp*nu_max);
   xb_max = Q2/(2*Mp*nu_min);

   //xb generation
   std::uniform_real_distribution<double> distxb(xb_min, xb_max);
   xb = distxb(gen);
   nu = Q2/(2*Mp*xb);


   //Kinematics of scatering electron
   E_electron = Eb - nu; //Energy of scattering electron
   theta_electron = 2*TMath::ASin(sqrt(Q2/(Eb*E_electron*4))); //Theta of scattering electron

   std::uniform_real_distribution<double> distphi_electron(0, 2*TMath::Pi());
   phi_electron = distphi_electron(gen);

   Electron.SetPxPyPzE(E_electron*TMath::Sin(theta_electron),0,E_electron*TMath::Cos(theta_electron),E_electron);
   Electron.Rotate(phi_electron,Beam.Vect());

   //Kinematics of virtual photon
   q = Beam - Electron;

   //Kinematics of proton 
   W = sqrt(Mp*Mp + 2*Mp*nu - Q2);
   t_min = phys.t_min_calcul(Q2, W);
   t_max = phys.t_max_calcul(Q2, W);
   //test with an other formula for t_min/max -> it's ok same value, t_min=t_min2 and t_max=t_max2
   //and find the ref for the formula in Byckling, E., and Kajantie, K. (1973b). Particle kinematics. Wiley-Interscience. (chapter 2->2 scattering)
   t_min2 = -Q2 + m_phi*m_phi - (1/(2*W*W))*((W*W - Q2 - Mp*Mp)*(W*W + m_phi*m_phi - Mp*Mp) + sqrt(phys.Lambda(W*W, -Q2, Mp*Mp))*sqrt(phys.Lambda(W*W, m_phi*m_phi, Mp*Mp)));
   t_max2 = -Q2 + m_phi*m_phi - (1/(2*W*W))*((W*W - Q2 - Mp*Mp)*(W*W + m_phi*m_phi - Mp*Mp) - sqrt(phys.Lambda(W*W, -Q2, Mp*Mp))*sqrt(phys.Lambda(W*W, m_phi*m_phi, Mp*Mp)));

   std::uniform_real_distribution<double> dist_t_proton(t_min, t_max);
   t = dist_t_proton(gen); // Mandelstan variable t = (p - p')^2 wich is the same for the 2->2 process (q + p -> phi + p') the 2->3 process (ep->e'p'phi)
   E_proton = (-t + 2*Mp*Mp)/(2*Mp);
   s_22 = -Q2 + Mp*Mp + 2*q.E()*Mp; // Mandelstan variable s for the 2->2 process (q + p -> phi + p')
   u_22 = -Q2 + 2*Mp*Mp + m_phi*m_phi - t - s_22; // Mandelstan u variable for the 2->2 process (q + p -> phi + p'). we use u + t + s = sum(mass^2)


   //find the ref for the next formula in Byckling, E., and Kajantie, K. (1973b). Particle kinematics. Wiley-Interscience. (chapter 2->2 scattering)
   theta_proton_gamma = TMath::ACos(((s_22 + Q2 - Mp*Mp)*(Mp*Mp + Mp*Mp - t) + 2*Mp*Mp*(u_22 + Q2 - Mp*Mp))/(sqrt(phys.Lambda(s_22, -Q2, Mp*Mp))*sqrt(phys.Lambda(t, Mp*Mp, Mp*Mp))));

   //axis relative to the photon
   TVector3 ez = q.Vect().Unit();               
   TVector3 ex = ez.Orthogonal().Unit();        
   TVector3 ey = ez.Cross(ex).Unit();  

   std::uniform_real_distribution<double> dist_phi_proton_gamma(0, 2*TMath::Pi());

   phi_proton_gamma = dist_phi_proton_gamma(gen); // phi relative to the photon axis 

   //calculate the direction of the proton vector with theta and phi 
   //important to use the axis relative to the photon 
   TVector3 proton_vect =
    ex * (sin(theta_proton_gamma)*cos(phi_proton_gamma)) +
    ey * (sin(theta_proton_gamma)*sin(phi_proton_gamma)) +
    ez * cos(theta_proton_gamma);

   proton_vect *= sqrt(E_proton*E_proton - Mp*Mp); //norme of momemtum

   Proton.SetVectM(proton_vect, Mp); // defintion of the quadri-vector

   //Kinematics of Phi meson
   Phi = Target + q - Proton;

   //Kinematics of Ks and Kl;
   std::uniform_real_distribution<double> dist_costheta_ks(-1, 1); // important to generate cos(theta) instead of theta beacause d(solide_angle) = sin(theta)dtheta * dphi
   std::uniform_real_distribution<double> dist_phi_ks(0, 2*TMath::Pi());

   theta_ks = TMath::ACos(dist_costheta_ks(gen));
   phi_ks = dist_phi_ks(gen);

   TVector3 pKs_CMS;
      pKs_CMS.SetMagThetaPhi(sqrt((m_phi/2)*(m_phi/2) - m_kaons*m_kaons), theta_ks, phi_ks);

   TLorentzVector Ks(pKs_CMS, m_phi/2);
   TLorentzVector Kl(-pKs_CMS, m_phi/2);

   TVector3 boost_lab = Phi.BoostVector();
   Ks.Boost(boost_lab);
   Kl.Boost(boost_lab);

   //Kinematics of pi+ and pi-;
   std::uniform_real_distribution<double> dist_costheta_pi(-1, 1);
   std::uniform_real_distribution<double> dist_phi_pi(0, 2*TMath::Pi());

   theta_pi = TMath::ACos(dist_costheta_pi(gen));
   phi_pi = dist_phi_pi(gen);

   TVector3 pPi_CMS;
      pPi_CMS.SetMagThetaPhi(sqrt((m_kaons/2)*(m_kaons/2) - m_pions*m_pions), theta_pi, phi_pi);

   TLorentzVector Pip(pPi_CMS, m_kaons/2);
   TLorentzVector Pim(-pPi_CMS, m_kaons/2);

   TVector3 boost_lab2 = Ks.BoostVector();
   Pip.Boost(boost_lab2);
   Pim.Boost(boost_lab2);

   //Weight of PhaseSpace
   weight = abs(Q2_max - Q2_min)*abs(xb_max - xb_min)*abs(t_max-t_min);
   sum_weight_phasespace += weight;

   weight_crosssection = phys.dsigmadt_tot(t, Q2, W, Eb);
   finalweight = weight*weight_crosssection*BR_pions*BR_kaons;

   //Vertex of e- 
   std::uniform_real_distribution<double> dist_vz0(-5.5, -0.5);
   vz0 = dist_vz0(gen);
   vx=0.;
   vy=0.;
   vz=vz0;

   //Vertex of Ks

   Ks_px = Ks.Px();
   Ks_py = Ks.Py();
   Ks_pz = Ks.Pz();

   norme = sqrt(Ks_px*Ks_px + Ks_py*Ks_py + Ks_pz*Ks_pz); //magnitude of the momentum of the Ks

   //normalisation of Ks vector (just to have the direction of emission)
   Ks_px_norme = Ks_px/norme; 
   Ks_py_norme = Ks_py/norme;
   Ks_pz_norme = Ks_pz/norme; 

   //shifted the vertex in the direction calculate previously * 2.8 cm
   // Ks_vx, Ks_vy and Ks_vz are usefule to set the vertex of pi+ and pi- in Lund file
   Ks_vx = 0 + Ks_px_norme*cT_Ks;
   Ks_vy = 0 + Ks_py_norme*cT_Ks;
   Ks_vz = vz0 + Ks_pz_norme*cT_Ks;

   // Fill Lund 

   //header
   files[i] << 4 << setw(15) << 1 << setw(5) << 1 << setw(15) << 0 << setw(15) << 0 << setw(15) << 11 << setw(15) << Eb << setw(15) << 0 << setw(15) << indice_evenement << setw(15) << finalweight << endl;
  
   // proton scattering
   files[i] << 1 << setw(5) << 1 << setw(5) << 1 << setw(7) << 2212 << setw(5) << 0 << setw(5) << 0 << setw(15) << Proton.Px() << setw(15) << Proton.Py() << setw(15) << Proton.Pz();
   files[i] << setw(15) << Proton.E() << setw(15) << Mp << setw(15) << 0. << setw(15) << 0. << setw(15) << vz << endl;

   // electron scattering
   files[i] << 2 << setw(5) << 1 << setw(5) << 1 << setw(7) << 11 << setw(5) << 0 << setw(5) << 0 << setw(15) << Electron.Px() << setw(15) << Electron.Py() << setw(15) << Electron.Pz();
   files[i] << setw(15) << Electron.E() << setw(15) << m_electron << setw(15) << 0. << setw(15) << 0. << setw(15) << vz << endl;

   // PI+
   files[i] << 3 << setw(5) << 1 << setw(5) << 1 << setw(7) << 211 << setw(5) << 0 << setw(5) << 0 << setw(15) << Pip.Px() << setw(15) << Pip.Py() << setw(15) << Pip.Pz();
   files[i] << setw(15) << Pip.E() << setw(15) << m_pions << setw(15) << Ks_vx << setw(15) << Ks_vy << setw(15) << Ks_vz << endl;
  
   // PI-
   files[i] << 4 << setw(5) << 1 << setw(5) << 1 << setw(7) << -211 << setw(5) << 0 << setw(5) << 0 << setw(15) << Pim.Px() << setw(15) << Pim.Py() << setw(15) << Pim.Pz();
   files[i] << setw(15) << Pim.E() << setw(15) << m_pions << setw(15) << Ks_vx << setw(15) << Ks_vy << setw(15) << Ks_vz << endl;


   //Fill histograms

   pvstheta_el->Fill(Electron.Theta()*(180/3.14), Electron.P(), finalweight);
   pvstheta_pr->Fill(Proton.Theta()*(180/3.14), Proton.P(), finalweight);
   pvstheta_kl->Fill(Kl.Theta()*(180/3.14), Kl.P(), finalweight);
   pvstheta_ks->Fill(Ks.Theta()*(180/3.14), Ks.P(), finalweight);
   pvstheta_pip->Fill(Pip.Theta()*(180/3.14), Pip.P(), finalweight);
   pvstheta_pim->Fill(Pim.Theta()*(180/3.14), Pim.P(), finalweight);

   if (Electron.P() > 5.25 && Electron.P() < 8 && Proton.P() > 0.45){ //same cut than Bhawani plot 

      tvsQ2->Fill(t, Q2, finalweight);
      Q2vsW->Fill(Q2, W, finalweight);
      histt->Fill(t, finalweight);
      hist_Q2->Fill(Q2, finalweight);
      hist_xb->Fill(xb, finalweight);
      hist_W->Fill(W, finalweight);

      hist_p_electron->Fill(Electron.P(), finalweight);
      hist_p_proton->Fill(Proton.P(), finalweight);
      hist_p_Ks->Fill(Ks.P(), finalweight);
      hist_p_Kl->Fill(Kl.P(), finalweight);

      hist_theta_electron->Fill(Electron.Theta(), finalweight);
      hist_theta_proton->Fill(Proton.Theta(), finalweight);
      hist_theta_Ks->Fill(Ks.Theta(), finalweight);
      hist_theta_Kl->Fill(Kl.Theta(), finalweight);

      hist_phi_electron->Fill(Electron.Phi(), finalweight);
      hist_phi_proton->Fill(Proton.Phi(), finalweight);
      hist_phi_Ks->Fill(Ks.Phi(), finalweight);
      hist_phi_Kl->Fill(Kl.Phi(), finalweight);

   }

   Minv_pip_pim->Fill((Pip + Pim).M(), finalweight);
   Minv_pip_pim_kl->Fill((Pip + Pim + Kl).M(), finalweight);
   hist_vz->Fill(vz, finalweight);
   hist_Ks_vx->Fill(Ks_vx, finalweight);
   hist_Ks_vy->Fill(Ks_vy, finalweight);
   hist_Ks_vz->Fill(Ks_vz, finalweight);
   hist_dist_vevks->Fill(sqrt((Ks_vx - vx)*(Ks_vx - vx) + (Ks_vy - vy)*(Ks_vy - vy) + (Ks_vz - vz)*(Ks_vz - vz)), finalweight);


   //Verification -- Tests

   if (std::isnan(weight)){ // Verify if there is NaN
      cpt_nan += 1;
   }
   if (finalweight < 0){ 
      cpt_nan += 1;
   }

   //momentum conservation test 
   Px_sum = Pip.Px() + Pim.Px() + Kl.Px() + Proton.Px() + Electron.Px();
   Py_sum = Pip.Py() + Pim.Py() + Kl.Py() + Proton.Py() + Electron.Py();
   Pz_sum = Pip.Pz() + Pim.Pz() + Kl.Pz() + Proton.Pz() + Electron.Pz();

   //energy conservation test
   E_sum = Pip.E() + Pim.E() + Kl.E() + Proton.E() + Electron.E();

   if(n % 500000 == 0){

   cout << "\n--- Resume event---" << endl;
   cout << "s 2->3 : " << s_23 << endl;
   cout << "Q2 : " << Q2 << endl;
   cout << "W2_min : " << W2_min << endl;
   cout << "W2_max : " << W2_max << endl;
   cout << "xb : " << xb << endl;

   cout << "Scatering Electron : " << endl;
   cout << "E scatering electron : " << E_electron << endl;
   cout << "theta scatering electron : " << theta_electron << endl;
   cout << "ELECTRON Px : " << Electron.Px() << endl;
   cout << "ELECTRON Py : " << Electron.Py() << endl;
   cout << "ELECTRON Pz : " << Electron.Pz() << endl;
   cout << "ELECTRON E : " << Electron.E() << endl;

   cout << "q virtual photon : " << endl;
   cout << "photon Px : " << q.Px() << endl;
   cout << "photon Py : " << q.Py() << endl;
   cout << "photon Pz : " << q.Pz() << endl;
   cout << "photon E : " << q.E() << endl;

   cout << "W : " << W << endl;
   cout << "t_min : " << t_min << endl;
   cout << "t_max : " << t_max << endl;
   cout << "t_min2 : " << t_min2 << endl;
   cout << "t_max2 : " << t_max2 << endl;
   cout << "t : " << t << endl;
   cout << "E_proton : " << E_proton << endl;
   cout << "Theta Proton -- Virtual photon : " << theta_proton_gamma << endl;

   cout << "Scatering Proton : " << endl;
   cout << "Proton Px : " << Proton.Px() << endl;
   cout << "Proton Py : " << Proton.Py() << endl;
   cout << "Proton Pz : " << Proton.Pz() << endl;
   cout << "Proton E : " << Proton.E() << endl;

   cout << "Phi : " << endl;
   cout << "Phi Px : " << Phi.Px() << endl;
   cout << "Phi Py : " << Phi.Py() << endl;
   cout << "Phi Pz : " << Phi.Pz() << endl;
   cout << "Phi E : " << Phi.E() << endl;

   cout << "Ks : " << endl;
   cout << "Ks Px : " << Ks.Px() << endl;
   cout << "Ks Py : " << Ks.Py() << endl;
   cout << "Ks Pz : " << Ks.Pz() << endl;
   cout << "Ks E : " << Ks.E() << endl;

   cout << "Kl : " << endl;
   cout << "Kl Px : " << Kl.Px() << endl;
   cout << "Kl Py : " << Kl.Py() << endl;
   cout << "Kl Pz : " << Kl.Pz() << endl;
   cout << "Kl E : " << Kl.E() << endl;

   cout << "pi+ : " << endl;
   cout << "pi+ Px : " << Pip.Px() << endl;
   cout << "pi+ Py : " << Pip.Py() << endl;
   cout << "pi+ Pz : " << Pip.Pz() << endl;
   cout << "pi+ E : " << Pip.E() << endl;

   cout << "pi- : " << endl;
   cout << "pi- Px : " << Pim.Px() << endl;
   cout << "pi- Py : " << Pim.Py() << endl;
   cout << "pi- Pz : " << Pim.Pz() << endl;
   cout << "pi- E : " << Pim.E() << endl;

   // Px_sum and Py_sum expected to 0.0 beceause there no momentum in x or y initialy
   cout << "Sum of x final state momentum :  " << Px_sum << endl;
   cout << "Sum of y final state momentum :  " << Py_sum << endl;
   cout << "Sum of z final state momentum :  " << Pz_sum << endl;

   cout << "Sum of E final state momentum :  " << E_sum << endl;

   cout << "Weight phasespace event : " << weight << endl;
   cout << "finalweight phasespace event : " << finalweight << endl;
   cout << "----------------\n" << endl;
   }



   } // end of loop in events
   
   files[i].close();

   } // end of loop in lund file

   cout << "Writting pdf..." <<endl;
   TString pdfFile = "phi_generator_plots.pdf";
   TCanvas *c = new TCanvas("c", "Plots", 800, 600);

   pvstheta_el->Draw("COLZ"); c->Print(pdfFile + "("); 
   pvstheta_pr->Draw("COLZ"); c->Print(pdfFile);
   pvstheta_kl->Draw("COLZ"); c->Print(pdfFile);
   pvstheta_ks->Draw("COLZ"); c->Print(pdfFile);
   pvstheta_pip->Draw("COLZ"); c->Print(pdfFile);
   pvstheta_pim->Draw("COLZ"); c->Print(pdfFile);
   Minv_pip_pim->Draw("HIST"); c->Print(pdfFile);
   Minv_pip_pim_kl->Draw("HIST"); c->Print(pdfFile);
   hist_vz->Draw("HIST"); c->Print(pdfFile);
   hist_Ks_vx->Draw("HIST"); c->Print(pdfFile);
   hist_Ks_vy->Draw("HIST"); c->Print(pdfFile);
   hist_Ks_vz->Draw("HIST"); c->Print(pdfFile);
   hist_dist_vevks->Draw("HIST"); c->Print(pdfFile);
   hist_p_electron->Draw("HIST"); c->Print(pdfFile);
   hist_p_proton->Draw("HIST"); c->Print(pdfFile);
   hist_p_Ks->Draw("HIST"); c->Print(pdfFile);
   hist_p_Kl->Draw("HIST"); c->Print(pdfFile);
   hist_theta_electron->Draw("HIST"); c->Print(pdfFile);
   hist_theta_proton->Draw("HIST"); c->Print(pdfFile);
   hist_theta_Ks->Draw("HIST"); c->Print(pdfFile);
   hist_theta_Kl->Draw("HIST"); c->Print(pdfFile);
   hist_phi_electron->Draw("HIST"); c->Print(pdfFile);
   hist_phi_proton->Draw("HIST"); c->Print(pdfFile);
   hist_phi_Ks->Draw("HIST"); c->Print(pdfFile);
   hist_phi_Kl->Draw("HIST"); c->Print(pdfFile);
   hist_Q2->Draw("HIST"); c->Print(pdfFile);
   hist_xb->Draw("HIST"); c->Print(pdfFile);
   histt->Draw("HIST"); c->Print(pdfFile);
   hist_W->Draw("HIST"); c->Print(pdfFile);
   Q2vsW->Draw("HIST"); c->Print(pdfFile);
   tvsQ2->Draw("COLZ");

   c->Print(pdfFile + ")");

   delete c;

   //cout << "Sum phasespace weight : " << sum_weight_phasespace << endl;
   cout << "number of nan : " << cpt_nan << endl;

} //end of generator function

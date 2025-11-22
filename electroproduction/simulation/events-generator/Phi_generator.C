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


// Distribution E_gamma



// Classe section efficace

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

        // t-depence dsigmaL/dt and dsigmaT/dt + sum avec episolon qui depent de E -> dsigma tot/dt pour le processus gamma* + p -> phi + p : 

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

        // cross section total dsigma / dQ^2 dt xb partie emission du photon QED * gamma* + p -> phi + p (QCD) ce qui done bien ep->e'p'phi

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


// Generateur 
void Phi_generator() {


   double sum_weight_phasespace = 0;
   double cpt_nan = 0;


   for (int i = 0; i < 2000000; ++i) {

   //Cross section and some useful functions
   Physics phys;

   //Mass
   double m_electron = 0.00051;
   double Mp = 0.938;
   double m_kaons = 0.4937;
   double m_phi = 1.019;
   double m_pions = 0.1395;

   //Kinamitics varibale
   double Q2_min = 0.1;
   double Q2_max = 7;
   double Eb = 10.6;

   double s_23 = m_electron*m_electron + 2*Mp*Eb + Mp*Mp;
   double Q2, W, W2, W2_min, W2_max, nu, nu_min, nu_max, xb, xb_min, xb_max, theta_electron, phi_electron, E_electron, 
   t, t_min, t_max, E_proton, s_22, u_22, theta_proton_gamma, phi_proton_gamma, theta_ks, phi_ks, theta_pi, phi_pi;

   //Weight of PhaseSpace
   double weight;


   //Initial quadrivector
   TLorentzVector Target(0.0, 0.0, 0.0, Mp);  
   TLorentzVector Beam(0.0, 0.0, Eb, sqrt(Eb*Eb + m_electron*m_electron));

   TLorentzVector Proton, Electron, q, Phi;

   // Generator
   std::random_device rd;
   std::mt19937 gen(rd()); 

   //Q2 generation
   std::uniform_real_distribution<double> distQ2(Q2_min, Q2_max);
   Q2 = distQ2(gen);
   
   //W2 phase space
   W2_min = pow((m_phi + Mp), 2);
   W2_max = s_23 + m_electron*m_electron - (1/(2*m_electron*m_electron))*(( s_23 + m_electron*m_electron - Mp*Mp )*(2*m_electron*m_electron + Q2) - sqrt(phys.Lambda(s_23, m_electron*m_electron, Mp*Mp))*sqrt(phys.Lambda(-Q2, m_electron*m_electron, m_electron*m_electron)));

   //nu phase space
   nu_min = (W2_min + Q2 - Mp*Mp)/(2*Mp);
   nu_max = (W2_max + Q2 - Mp*Mp)/(2*Mp);

   //xb phase space 
   xb_min = Q2/(2*Mp*nu_max);
   xb_max = Q2/(2*Mp*nu_min);

   std::uniform_real_distribution<double> distxb(xb_min, xb_max);
   xb = distxb(gen);
   nu = Q2/(2*Mp*xb);


   //Kinematics of scatering electron
   E_electron = Eb - nu; //Energy
   theta_electron = 2*TMath::ASin(sqrt(Q2/(Eb*E_electron*4))); //Theta

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

   std::uniform_real_distribution<double> dist_t_proton(t_min, t_max);
   t = dist_t_proton(gen);
   E_proton = (-t + 2*Mp*Mp)/(2*Mp);
   s_22 = -Q2 + Mp*Mp + 2*q.E()*Mp;
   u_22 = -Q2 + 2*Mp*Mp + m_phi*m_phi - t - s_22;

   theta_proton_gamma = TMath::ACos(((s_22 + Q2 - Mp*Mp)*(Mp*Mp + Mp*Mp - t) + 2*Mp*Mp*(u_22 + Q2 - Mp*Mp))/(sqrt(phys.Lambda(s_22, -Q2, Mp*Mp))*sqrt(phys.Lambda(t, Mp*Mp, Mp*Mp))));

   TVector3 ez = q.Vect().Unit();               
   TVector3 ex = ez.Orthogonal().Unit();        
   TVector3 ey = ez.Cross(ex).Unit();  

   std::uniform_real_distribution<double> dist_phi_proton_gamma(0, 2*TMath::Pi());

   phi_proton_gamma = dist_phi_proton_gamma(gen);

   TVector3 proton_vect =
    ex * (sin(theta_proton_gamma)*cos(phi_proton_gamma)) +
    ey * (sin(theta_proton_gamma)*sin(phi_proton_gamma)) +
    ez * cos(theta_proton_gamma);

   proton_vect *= sqrt(E_proton*E_proton - Mp*Mp);

   Proton.SetVectM(proton_vect, Mp);

   //Kinematics of Phi meson
   Phi = Target + q - Proton;

   //Kinematics of Ks and Kl;
   std::uniform_real_distribution<double> dist_costheta_ks(-1, 1);
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

   if (std::isnan(weight)){
   } else {
      sum_weight_phasespace += weight;
   }

   if (std::isnan(weight)){
      cpt_nan += 1;
   }

   if(i % 2000000 == 0){

   cout << "\n--- Résumé evenement---" << endl;
   cout << "s 23 : " << s_23 << endl;
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

   cout << "Weight event : " << weight << endl;
   cout << "----------------\n" << endl;

   }

   }
   cout << "Somme weight event : " << sum_weight_phasespace << endl;
   cout << "number of nan : " << cpt_nan << endl;
   

}

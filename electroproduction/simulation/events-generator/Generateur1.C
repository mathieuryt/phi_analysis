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
#include <iostream>
#include <cmath>
using namespace std;



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
void Generateur1(double nb_fichier, double nb_event) {

   //paramètres 
   double Eb = 10.6;    

   // plage de calcul de la distrib
   
   int points = 1000;

   //masse particules
   double me = 0.00051;
   double m_kaons = 0.4937;
   double m_phi = 1.019;
   double m_pi = 0.1395;
   double Mp = 0.938;

   // compteur de particules
   double counter = 0;
   double N_retenus = 0;
   double strange_weight = 0;

   Physics phys;

   // paramettres section efficace


   //ini fichiers
   std::vector<std::ofstream> files(nb_fichier);

   // Définition du proton cible au repos
   TLorentzVector target(0.0, 0.0, 0.0, Mp);  
   TLorentzVector beam(0.0, 0.0, Eb, Eb);
  
   Double_t masses[3] = {Mp, m_phi, me}; // masse proton et phi et e-

   TGenPhaseSpace event;

   // initialisation des histo

   double bin1 = 200;
   double bin2 = 100;

   TH2F *pvstheta_el = new TH2F("pvstheta_el", "p vs theta for electron", bin1, 0, 90, bin1, 0.0, 11);
   TH2F *pvstheta_pr = new TH2F("pvstheta_pr", "p vs theta for proton", bin1, 0, 90, bin1, 0.0, 11);
   TH2F *pvstheta_kl = new TH2F("pvstheta_kl", "p vs theta for kl", bin1, 0, 90, bin1, 0.0, 11); 
   TH2F *pvstheta_pip = new TH2F("pvstheta_pip", "p vs theta for pi+", bin1, 0, 90, bin1, 0.0, 11); 
   TH2F *pvstheta_pim = new TH2F("pvstheta_pim", "p vs theta for pi-", bin1, 0, 90, bin1, 0.0, 11); 

   TH1F *Minv_pip_pim = new TH1F("Minv_pip_pim", "Invariant mass of pi+ pi-", bin2, 0, 1.0);
   TH1F *Minv_pip_pim_kl = new TH1F("Minv_pip_pim_kl", "Invariant mass of (pi+ pi- + Kl)", 400, 0.8, 1.2);
   TH1F *hist_vz = new TH1F("hist_vz", "vz of e-", bin2, -6, 1.0);
   TH1F *hist_Ks_vx = new TH1F("hist_Ks_vx", "vertex x of Ks", bin2, -3, 3.0);
   TH1F *hist_Ks_vy = new TH1F("hist_Ks_vy", "vertex y of Ks", bin2, -3, 3.0);
   TH1F *hist_Ks_vz = new TH1F("hist_Ks_vz", "vertex z of Ks", bin2, -4, 4.0);
   TH1F *hist_dist_vevks = new TH1F("hist_dist_vevks", "distance (x, y, z) between vertex e and vertex Ks", bin2, 2, 4);
   TH1F *histt = new TH1F("histt", "t", bin2, -8, 0);
   TH1F *hist_Q2 = new TH1F("hist_Q2", "Q2", bin2, 0, 8);

   auto setTitles = [](TH1 *h, const char *xtitle, const char *ytitle) {
               h->GetXaxis()->SetTitle(xtitle);
               h->GetYaxis()->SetTitle(ytitle);
               h->GetXaxis()->CenterTitle();
               h->GetYaxis()->CenterTitle();
   };

   setTitles(pvstheta_el, "#theta [degree]", "p [GeV]");
   setTitles(pvstheta_pr, "#theta [degree]", "p [GeV]");
   setTitles(pvstheta_kl, "#theta [degree]", "p [GeV]");
   setTitles(pvstheta_pip, "#theta [degree]", "p [GeV]");
   setTitles(pvstheta_pim, "#theta [degree]", "p [GeV]");

   Minv_pip_pim->GetXaxis()->SetTitle("Minv pi+ pi- [GeV]");
   Minv_pip_pim_kl->GetXaxis()->SetTitle("Minv pi+ pi- + kl [GeV]");
   hist_vz->GetXaxis()->SetTitle("vz e- [cm]");
   hist_Ks_vx->GetXaxis()->SetTitle("vx Ks [cm]");
   hist_Ks_vy->GetXaxis()->SetTitle("vy ks [cm]");
   hist_Ks_vz->GetXaxis()->SetTitle("vz Ks [cm]");
   hist_dist_vevks->GetXaxis()->SetTitle("distance vertex e- and vertex Ks [cm]");
   histt->GetXaxis()->SetTitle("t [GeV^2]");
   hist_Q2->GetXaxis()->SetTitle("Q2 [GeV^2]");
   


   // Création et ouverture des fichiers
   for (int i = 0; i < nb_fichier; ++i) {
       std::string filename = "PhiGen_" + std::to_string(i) + ".txt";
       files[i].open(filename);

       if (!files[i]) {
           std::cerr << "Erreur lors de l'ouverture du fichier : " << filename << std::endl;
       }


       for (Int_t n = 0; n < nb_event; n++) {

         TLorentzVector Z = beam + target;   
         event.SetDecay(Z, 3, masses); // decay e- et proton ---> phi + proton + e-
         Double_t weight = event.Generate();
  
         TLorentzVector *pProton = event.GetDecay(0); // Proton final
         TLorentzVector *pPhi = event.GetDecay(1); //  phi
         TLorentzVector *pEl = event.GetDecay(2); //  electron
         
  
         Double_t t = 2*(Mp*Mp - Mp*pProton->E()); // calcul de t d'une autre maniere (renvoi bien t)

         TLorentzVector q = beam - *pEl;


         double Q2 = -q.M2();
         double W = (target + q).M();
  
         double weight_crosssection = phys.dsigmadt_tot(t, Q2, W, Eb);
  
  
         // verification de la conservation de l'énergie
         Double_t E_ini = Eb + 0.938;
         Double_t E_fin = pPhi->E() + pProton->E() + pEl->E();
  
         // verification que somme impulsion ini = somme impulsion final
         // normalement tout les P_diff valent 0
  
         Double_t Px_diff = 0 + 0 - pProton->Px() - pPhi->Px() - pEl->Px();
         Double_t Py_diff = 0 + 0 - pProton->Py() - pPhi->Py() - pEl->Py();
         Double_t Pz_diff = Eb + 0 - pProton->Pz() - pPhi->Pz() -pEl->Pz();


  
         // histo fill en pondérant par Poids (celui des E_gamma) et weigh_crosssection (qui vient du modèle de section efficace)
  
         // j'ai pas mis les poids "weight" car pour l'instant ils valent 1.
  
  
  
  
         double s = Mp*Mp + 2*Mp*Eb;
  
         //double t_min = phys.T_min(0., Mp*Mp, 0., Mp*Mp, s);
         //double t_max = phys.T_max(0., Mp*Mp, 0., Mp*Mp, s);
  
         //double correction_weight = abs(t_min - t_max);
   
  
         double tot_weight = weight*weight_crosssection;
  
  
  
         // DEUXIEME DECAY : désintégration du PHI --> Ks+Kl

  
         Double_t masses_K[2] = { m_kaons, m_kaons };
         TGenPhaseSpace decay;
         decay.SetDecay(*pPhi, 2, masses_K);
  
         Double_t weight2 = decay.Generate();
  
         TLorentzVector *pKs = decay.GetDecay(0);
         TLorentzVector *pKl = decay.GetDecay(1);

         double BR_kaons = 0.34;

         double tot_weight2 = weight2*BR_kaons;


         double angle = pKs->Vect().Angle(pKs->Vect());
         double angle_deg = angle * 180.0 / 3.14;

         double delta_theta = abs(pKs->Theta() - pKl->Theta());
         double delta_theta_deg = delta_theta * 180.0 / 3.14;

         double delta_phi = pKs->Phi() - pKl->Phi();
         double delta_phi_deg = delta_phi * 180.0 / 3.14;

         if (delta_phi_deg>180){

            delta_phi_deg = 360 - delta_phi_deg;

         }

         if (delta_phi_deg<-180){
            

            delta_phi_deg = -(delta_phi_deg + 360);

         }
  
         // histo p et theta
  
         Double_t theta_phi = pPhi->Theta();
         Double_t theta_Ks = pKs->Theta();
         Double_t theta_Kl = pKl->Theta();
         Double_t theta_proton = pProton->Theta(); // proton difuse
         Double_t theta_electron = pEl->Theta(); // electron difuse
  
  
         Double_t P_phi = pPhi->P();
         Double_t P_Ks = pKs->P();
         Double_t P_Kl = pKl->P();
         Double_t P_proton = pProton->P(); // proton difuse
         Double_t P_electron = pEl->P(); // electron difuse

         
         // 3EME DECAYS Ks --> pi+ + pi-
         
         Double_t masses_pi[2] = { m_pi, m_pi };
         TGenPhaseSpace decay3;
         decay3.SetDecay(*pKs, 2, masses_pi);
  
         Double_t weight3 = decay3.Generate();
  
         TLorentzVector *pPip = decay3.GetDecay(0);
         TLorentzVector *pPim = decay3.GetDecay(1);

         double BR_pions = 0.692;

         double tot_weight3 = weight3*BR_pions;

         Double_t theta_pip = pPip->Theta();
         Double_t theta_pim = pPim->Theta();

         Double_t P_pip = pPip->P();
         Double_t P_pim = pPim->P();


         // calcul du decalage a cause du temps de vol du Ks

         Double_t Ks_px = pKs->Px();
         Double_t Ks_py = pKs->Py(); 
         Double_t Ks_pz = pKs->Pz(); 

         double norme = sqrt(Ks_px*Ks_px + Ks_py*Ks_py + Ks_pz*Ks_pz);

         double Ks_px_norme = Ks_px/norme;
         double Ks_py_norme = Ks_py/norme;
         double Ks_pz_norme = Ks_pz/norme; 

         // le vecteur (Ks_px_norme, Ks_py_norme, Ks_pz_norme) est de norme 1 et pointe dans la direction de p_Ks

         double cT_Ks = 2.8; // longueur de vol du Ks si unite bien en [cm]

         // Générateur aléatoire basé sur le matériel (seed)
         std::random_device rd;
         std::mt19937 gen(rd());  // Mersenne Twister (générateur)
         std::uniform_real_distribution<> dist(-5.5, -0.5); // pour un réel
         // std::uniform_int_distribution<> dist(-5, 0);   // pour un entier

         double vz0 = dist(gen);
   

         Double_t Ks_vx, Ks_vy, Ks_vz; // pour Ks (peut etre Kl)

         Ks_vx = 0 + Ks_px_norme*cT_Ks;
         Ks_vy = 0 + Ks_py_norme*cT_Ks;
         Ks_vz = vz0 + Ks_pz_norme*cT_Ks;

         Double_t vx, vy, vz; // pour e- et p difusé
         vx=0.;
         vy=0.;
         vz=vz0;

         double finalweight = tot_weight*tot_weight2*tot_weight3;

         double nu_test1 = (W*W - Mp*Mp + Q2)/(2*Mp);
         Double_t nu_test2 = Eb - pEl->E();



         if (n % 1000000000 == 0) {

          std::cout << "\n--- Résumé evenement---" << std::endl;
          std::cout << "Événement " << n << " : Eb = " << Eb << std::endl;
          std::cout << "Événement " << n << " : masse du phi = " << pPhi->M() << std::endl;
          std::cout << "Événement " << n << " : masse p diffuse = " << pProton->M() << std::endl;
          std::cout << "Événement " << n << " : masse e diffuse = " << pEl->M() << std::endl;

          std::cout << "t et Q2" << std::endl;
          std::cout << "Événement " << n << " : t = " << t << std::endl;
          std::cout << "Événement " << n << " : Q2 = " << Q2 << std::endl;
          std::cout << "Événement " << n << " : W = " << W << std::endl;
          std::cout << "Événement " << n << " : nu1 = " << nu_test1 << std::endl;
          std::cout << "Événement " << n << " : n2 = " << nu_test2 << std::endl;

          std::cout << "weight tgenphase" << std::endl;
          std::cout << "Événement " << n << " : weight = " << weight << std::endl;

          std::cout << "Verification conservation de l'énergie " << std::endl;
          std::cout << "Événement " << n << " : E ini = " << E_ini << std::endl;
          std::cout << "Événement " << n << " : E fini = " << E_fin << std::endl;

          std::cout << "Verification somme impulsion ini = somme impulsion fini " << std::endl;
          std::cout << "Événement " << n << " : Px_diff = " << Px_diff << std::endl;
          std::cout << "Événement " << n << " : Py_diff = " << Py_diff << std::endl;
          std::cout << "Événement " << n << " : Pz_diff = " << Pz_diff << std::endl;

          std::cout << "pPhi: (E=" << pPhi->E() 
          
            << ", px=" << pPhi->Px() 
            << ", py=" << pPhi->Py() 
            << ", pz=" << pPhi->Pz() 
            << ", M=" << pPhi->M() << ")" << std::endl;

         std::cout << "Événement " << n << " : poids premiere decay = " << weight << std::endl;
         cout << "Section efficace : " << weight_crosssection << " nb " << endl;
         std::cout << "Événement " << n << " : poids deuxieme decay = " << weight2 << std::endl;
         std::cout << "Événement " << n << " : poids troisieme decay = " << weight3 << std::endl;
         std::cout << "Événement " << n << " : t_min_calcul = " << phys.t_min_calcul(Q2, W) << std::endl;

         std::cout << "----------------\n" << std::endl;
  
         }


         if (finalweight < 0){
            strange_weight += 1;

            std::cout << "\n--- Résumé evenement ETRANGE---" << std::endl;
            std::cout << "Événement " << n << " : Eb = " << Eb << std::endl;
            std::cout << "Événement " << n << " : masse du phi = " << pPhi->M() << std::endl;
            std::cout << "Événement " << n << " : masse p diffuse = " << pProton->M() << std::endl;
            std::cout << "Événement " << n << " : masse e diffuse = " << pEl->M() << std::endl;

            std::cout << "t et Q2" << std::endl;
            std::cout << "Événement " << n << " : t = " << t << std::endl;
            std::cout << "Événement " << n << " : Q2 = " << Q2 << std::endl;
            std::cout << "Événement " << n << " : W = " << W << std::endl;

            std::cout << "weight tgenphase" << std::endl;
            std::cout << "Événement " << n << " : weight = " << weight << std::endl;

            std::cout << "Verification conservation de l'énergie " << std::endl;
            std::cout << "Événement " << n << " : E ini = " << E_ini << std::endl;
            std::cout << "Événement " << n << " : E fini = " << E_fin << std::endl;

            std::cout << "Verification somme impulsion ini = somme impulsion fini " << std::endl;
            std::cout << "Événement " << n << " : Px_diff = " << Px_diff << std::endl;
            std::cout << "Événement " << n << " : Py_diff = " << Py_diff << std::endl;
            std::cout << "Événement " << n << " : Pz_diff = " << Pz_diff << std::endl;

            std::cout << "pPhi: (E=" << pPhi->E() 
          
               << ", px=" << pPhi->Px() 
               << ", py=" << pPhi->Py() 
               << ", pz=" << pPhi->Pz() 
               << ", M=" << pPhi->M() << ")" << std::endl;

            std::cout << "Événement " << n << " : poids premiere decay = " << weight << std::endl;
            cout << "Section efficace : " << weight_crosssection << " nb " << endl;
            std::cout << "Événement " << n << " : poids deuxieme decay = " << weight2 << std::endl;
            std::cout << "Événement " << n << " : poids troisieme decay = " << weight3 << std::endl;
            std::cout << "Événement " << n << " : final weight = " << finalweight << std::endl;
            std::cout << "Événement " << n << " : t_min_calcul = " << phys.t_min_calcul(Q2, W) << std::endl;

            std::cout << "Événement " << n << " : sigmaT = " << phys.sigmaT(W, Q2) << std::endl;
            std::cout << "Événement " << n << " : R = " << phys.R(Q2) << std::endl;
            std::cout << "Événement " << n << " : dsumdt = " << phys.dsumdt(t, Q2, W, Eb) << std::endl;
            std::cout << "Événement " << n << " : dsigmadt_tot = " << phys.dsigmadt_tot(t, Q2, W, Eb) << std::endl;
            

            std::cout << "----------------\n" << std::endl;

         }
         
         //histo

         if(t < -0.1 && t > -8 && Q2 > 0.5 && Q2 < 8){

            N_retenus += 1;

            pvstheta_el->Fill(theta_electron*(180/3.14), P_electron, finalweight);
            pvstheta_pr->Fill(theta_proton*(180/3.14), P_proton, finalweight);
            pvstheta_kl->Fill(theta_Kl*(180/3.14), P_Kl, finalweight);
            pvstheta_pip->Fill(theta_pip*(180/3.14), P_pip, finalweight);
            pvstheta_pim->Fill(theta_pim*(180/3.14), P_pim, finalweight);

            Minv_pip_pim->Fill((*pPip + *pPim).M(), finalweight);
            Minv_pip_pim_kl->Fill((*pPip + *pPim + *pKl).M(), finalweight);
            hist_vz->Fill(vz, finalweight);
            hist_Ks_vx->Fill(Ks_vx, finalweight);
            hist_Ks_vy->Fill(Ks_vy, finalweight);
            hist_Ks_vz->Fill(Ks_vz, finalweight);
            hist_dist_vevks->Fill(sqrt((Ks_vx - vx)*(Ks_vx - vx) + (Ks_vy - vy)*(Ks_vy - vy) + (Ks_vz - vz)*(Ks_vz - vz)), finalweight);
            histt->Fill(t, finalweight);
            hist_Q2->Fill(Q2, finalweight);


            double indice_evenement = n; 
  
            //header
            files[i] << 4 << setw(15) << 1 << setw(5) << 1 << setw(15) << 0 << setw(15) << 0 << setw(15) << 11 << setw(15) << Eb << setw(15) << 0 << setw(15) << indice_evenement << setw(15) << finalweight << endl;
  
            // proton difusé
            files[i] << 1 << setw(5) << 1 << setw(5) << 1 << setw(7) << 2212 << setw(5) << 0 << setw(5) << 0 << setw(15) << pProton->Px() << setw(15) << pProton->Py() << setw(15) << pProton->Pz();
            files[i] << setw(15) << pProton->E() << setw(15) << Mp << setw(15) << 0. << setw(15) << 0. << setw(15) << vz << endl;

            // electron difusé
            files[i] << 2 << setw(5) << 1 << setw(5) << 1 << setw(7) << 11 << setw(5) << 0 << setw(5) << 0 << setw(15) << pEl->Px() << setw(15) << pEl->Py() << setw(15) << pEl->Pz();
            files[i] << setw(15) << pEl->E() << setw(15) << me << setw(15) << 0. << setw(15) << 0. << setw(15) << vz << endl;

            // PI+
            files[i] << 3 << setw(5) << 1 << setw(5) << 1 << setw(7) << 211 << setw(5) << 0 << setw(5) << 0 << setw(15) << pPip->Px() << setw(15) << pPip->Py() << setw(15) << pPip->Pz();
            files[i] << setw(15) << pPip->E() << setw(15) << m_pi << setw(15) << Ks_vx << setw(15) << Ks_vy << setw(15) << Ks_vz << endl;
  
            // PI-
            files[i] << 4 << setw(5) << 1 << setw(5) << 1 << setw(7) << -211 << setw(5) << 0 << setw(5) << 0 << setw(15) << pPim->Px() << setw(15) << pPim->Py() << setw(15) << pPim->Pz();
            files[i] << setw(15) << pPim->E() << setw(15) << m_pi << setw(15) << Ks_vx << setw(15) << Ks_vy << setw(15) << Ks_vz << endl;


         }
  
  
     }

     files[i].close();

   }

   cout << "Ecriture du pdf..." <<endl;

   double echelle_plot_couleur = 2.0;

   pvstheta_el->SetMinimum(0);  // échelle de couleur
   pvstheta_el->SetMaximum(echelle_plot_couleur);

   pvstheta_pr->SetMinimum(0);  // échelle de couleur
   pvstheta_pr->SetMaximum(echelle_plot_couleur);

   pvstheta_kl->SetMinimum(0);  // échelle de couleur
   pvstheta_kl->SetMaximum(echelle_plot_couleur);

   pvstheta_pip->SetMinimum(0);  // échelle de couleur
   pvstheta_pip->SetMaximum(echelle_plot_couleur);

   pvstheta_pim->SetMinimum(0);  // échelle de couleur
   pvstheta_pim->SetMaximum(echelle_plot_couleur);

   TString pdfFile = "generateur_plots.pdf";

   // Premier canvas (ouvre le PDF)
   TCanvas *c = new TCanvas("c", "Plots", 800, 600);

   pvstheta_el->Draw("COLZ");

   c->Print(pdfFile + "("); // ouvre le fichier PDF

   pvstheta_pr->Draw("COLZ"); 
   c->Print(pdfFile);

   pvstheta_kl->Draw("COLZ"); 
   c->Print(pdfFile);

   pvstheta_pip->Draw("COLZ"); 
   c->Print(pdfFile);

   pvstheta_pim->Draw("COLZ");
   c->Print(pdfFile);

   Minv_pip_pim->Draw("HIST"); c->Print(pdfFile);
   Minv_pip_pim_kl->Draw("HIST"); c->Print(pdfFile);
   hist_vz->Draw("HIST"); c->Print(pdfFile);
   hist_Ks_vx->Draw("HIST"); c->Print(pdfFile);
   hist_Ks_vy->Draw("HIST"); c->Print(pdfFile);
   hist_Ks_vz->Draw("HIST"); c->Print(pdfFile);
   hist_dist_vevks->Draw("HIST"); c->Print(pdfFile);
   histt->Draw("HIST"); c->Print(pdfFile);

   hist_Q2->Draw("HIST");
   c->Print(pdfFile + ")");

   delete c;

   std::cout << nb_fichier << "fichiers de " << nb_event << "events generés" << std::endl;
   std::cout << "apres les cuts sur t et Q2 le nombre d'événements retenus est : " << N_retenus << std::endl;
   std::cout << "Le nombre de weight bizarre (< 0) est : " << strange_weight << std::endl;
}

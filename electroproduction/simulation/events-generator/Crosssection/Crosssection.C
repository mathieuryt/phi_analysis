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


        double TAU(double theta){

        // test sut TAU

        double E = 2.344;
        double alpha = 1/137.036;
        double Mp = 0.938;
        double nu = 1.5;
        double y = nu/E;
        double Q2 = 4*E*(E-nu)*sin(theta/2)*sin(theta/2);
        double xb = Q2/(2*Mp*nu);
        double epsilon = (1 - y - (Q2/(4*E*E)) )/(1- y + (y*y)/2 + (Q2/(4*E*E)) );
        double tau = (alpha * Q2 * (1-xb))/(8*3.14*Mp*Mp*E*E*xb*xb*xb*(1-epsilon));

        return tau;

       }


   };





// Generateur 
void Crosssection() {

    Physics class1;

    double Eb = 4.0;

    // SIGMA TRANSVERSE EN FONCTION DE Q2 à W fixe 

    std::ofstream file("/Users/mr282803/test/cs_test.txt");

    for (double Q2 = 0.1; Q2 <= 30.0; Q2 += 0.05) {
        double W = 2.5;
        file << Q2 << " " << class1.sigmaT(W, Q2) << "\n";
    }

    file.close();
    std::cout << "Données écrites dans cs_test.txt" << std::endl;




    // SIGMA TRANSVERSE EN FONCTION DE W à Q2 fixe 

    std::ofstream file2("/Users/mr282803/test/cs_test2.txt");

    for (double W = 2; W <= 100.0; W += 0.1) {
        double Q2 = 2.5;
        file2 << W << " " << class1.sigmaT(W, Q2) << "\n";
    }

    file2.close();
    std::cout << "Données écrites dans cs_test2.txt" << std::endl;



    // ratio SIGMA TRANSVERSE SIGMA LONGITUDINALE en fonction de Q2

    std::ofstream file3("/Users/mr282803/test/cs_test3.txt");

    for (double Q2 = 0.1; Q2 <= 40.0; Q2 += 0.1) {

        file3 << Q2 << " " << class1.R(Q2) << "\n";
    }

    file3.close();
    std::cout << "Données écrites dans cs_test3.txt" << std::endl;




    // d sigma dt (somme des deux differentiel transverse et longitudinale) dsigmadt (gamma* p -> phi + p)

    double tmintest = class1.t_min_calcul(2.5, 2.5);
    std::cout << "pour W = 2.5 Gev et Q2 = 2.5 GeV2, t_min = " << tmintest << std::endl;


    std::ofstream file4("/Users/mr282803/test/cs_test4.txt");

    for (double t = -3.0; t <= tmintest; t += 0.01) {

        double t_print = tmintest - t;

        file4 << t_print << " " << class1.dsumdt(t, 2.5, 2.5, Eb) << "\n";
    }

    file4.close();
    std::cout << "Données écrites dans cs_test4.txt" << std::endl;





    // flux de photon * dsigam(gamma* p -> phi + p)/dt = dsigma(e p->phi + p)/dtdQ2dxb

    std::ofstream file5("/Users/mr282803/test/cs_test5.txt");

    for (double t = -3; t <= tmintest; t += 0.01) {

        double t_print = tmintest - t;

        file5 << t_print << " " << class1.dsigmadt_tot(t, 2.5, 2.5, Eb) << "\n";
    }

    file5.close();
    std::cout << "Données écrites dans cs_test5.txt" << std::endl;



    // TAU

    std::ofstream file6("/Users/mr282803/test/cs_test6.txt");

    for (double theta = 0; theta <= 0.2; theta += 0.001) {

        double theta_conv = theta*1000; //mrad

        file6 << theta_conv << " " << class1.TAU(theta) << "\n";
    }

    file6.close();
    std::cout << "Données écrites dans cs_test6.txt" << std::endl;


}

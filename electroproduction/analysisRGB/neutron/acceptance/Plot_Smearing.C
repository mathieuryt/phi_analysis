void Plot_Smearing()
{
    gStyle->SetOptStat(0);
    TString file_data_name = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_data2/fall2019_outbending_neutronKpKm_data_v2.root";
    TString file_mc_name   = "/Users/mr282803/Documents/analysisRGB/neutron/neutron_data2/fall2019_outbending_neutronKpKm_mc_v7.root";

    TFile *f_data = TFile::Open(file_data_name, "READ");
    TFile *f_mc   = TFile::Open(file_mc_name, "READ");

    TTree *t_data = (TTree*)f_data->Get("tree");
    TTree *t_mc   = (TTree*)f_mc->Get("tree");

    // -------------------------
    // Histograms
    // -------------------------
    TH1F *h_data = new TH1F("h_data", "t comparison; t [GeV^{2}]; Counts", 100, -10.0, 0.5);
    TH1F *h_mc   = new TH1F("h_mc",   "t comparison; t [GeV^{2}]; Counts", 100, -10.0, 0.5);
    TH1F *h_mcs  = new TH1F("h_mcs",  "t comparison; t [GeV^{2}]; Counts", 100, -10.0, 0.5);

    TH1F *h_data_minv = new TH1F("h_data_minv", "Invariant Mass K+ K-; Minv (K+K-) [GeV]; Counts", 100, 0.98, 1.1);
    TH1F *h_mc_minv   = new TH1F("h_mc_minv",   "Invariant Mass K+ K-; Minv (K+K-) [GeV]; Counts", 100, 0.98, 1.1);
    TH1F *h_mcs_minv  = new TH1F("h_mcs_minv",  "Invariant Mass K+ K-; Minv (K+K-) [GeV]; Counts", 100, 0.98, 1.1);

    h_data->SetLineColor(kBlack);
    h_mc->SetLineColor(kBlue+1);
    h_mcs->SetLineColor(kRed+1);

    h_data_minv->SetLineColor(kBlue);
    h_mc_minv->SetLineColor(kRed+1);
    h_mcs_minv->SetLineColor(kGreen+1);

    h_data->SetLineWidth(2);
    h_mc->SetLineWidth(2);
    h_mcs->SetLineWidth(2);

    h_data_minv->SetLineWidth(2);
    h_mc_minv->SetLineWidth(2);
    h_mcs_minv->SetLineWidth(2);

    h_mc_minv->Sumw2();
    h_mcs_minv->Sumw2();

    // -------------------------
    // DATA branches
    // -------------------------
    TLorentzVector *Missing_data = nullptr;
    TLorentzVector *Electron_data = nullptr;
    TLorentzVector *Neutron_data = nullptr;
    TLorentzVector *Kp = nullptr;
    TLorentzVector *Km = nullptr;

    double angle_data, Q2_data;
    double status_Kp_data, status_Km_data;
    double MinvKpKm_data;
    double t_data_var;

    t_data->SetBranchAddress("Missing", &Missing_data);
    t_data->SetBranchAddress("Electron", &Electron_data);
    t_data->SetBranchAddress("Neutron", &Neutron_data);
    t_data->SetBranchAddress("Kp", &Kp);
    t_data->SetBranchAddress("Km", &Km);
    t_data->SetBranchAddress("angle_neutron_missnucl", &angle_data);
    t_data->SetBranchAddress("Q2", &Q2_data);
    t_data->SetBranchAddress("Status_Kp", &status_Kp_data);
    t_data->SetBranchAddress("Status_Km", &status_Km_data);
    t_data->SetBranchAddress("MinvKpKm", &MinvKpKm_data);
    t_data->SetBranchAddress("t_missing_nucleon", &t_data_var);

    // -------------------------
    // MC branches
    // -------------------------
    TLorentzVector *Missing_mc = nullptr;
    TLorentzVector *Electron_mc = nullptr;
    TLorentzVector *Neutron_mc = nullptr;
    TLorentzVector *Kp_mc = nullptr;
    TLorentzVector *Km_mc = nullptr;
    TLorentzVector *Kp_mc_smear = nullptr;
    TLorentzVector *Km_mc_smear = nullptr;

    double angle_mc, Q2_mc;
    double status_Kp_mc, status_Km_mc;
    double MinvKpKm_mc;
    double t_missing_nucleon_mc;
    double t_missing_smeared_mc;
    double real_weight_SDME_effcorrection;

    t_mc->SetBranchAddress("Missing", &Missing_mc);
    t_mc->SetBranchAddress("Electron", &Electron_mc);
    t_mc->SetBranchAddress("Neutron", &Neutron_mc);
    t_mc->SetBranchAddress("Kp", &Kp_mc);
    t_mc->SetBranchAddress("Km", &Km_mc);
    t_mc->SetBranchAddress("Kp_smear", &Kp_mc_smear);
    t_mc->SetBranchAddress("Km_smear", &Km_mc_smear);
    t_mc->SetBranchAddress("angle_neutron_missnucl", &angle_mc);
    t_mc->SetBranchAddress("Q2", &Q2_mc);
    t_mc->SetBranchAddress("Status_Kp", &status_Kp_mc);
    t_mc->SetBranchAddress("Status_Km", &status_Km_mc);
    t_mc->SetBranchAddress("MinvKpKm", &MinvKpKm_mc);
    t_mc->SetBranchAddress("t_missing_nucleon", &t_missing_nucleon_mc);
    t_mc->SetBranchAddress("t_missing_smeared", &t_missing_smeared_mc);
    t_mc->SetBranchAddress("real_weight_SDME_effcorrection", &real_weight_SDME_effcorrection);

    // -------------------------
    // Fill DATA
    // -------------------------
    Long64_t n_data = t_data->GetEntries();

    for (Long64_t i = 0; i < n_data; i++) {
        t_data->GetEntry(i);

        if (!Missing_data || !Electron_data || !Neutron_data) continue;

        bool cut =
            angle_data * 180. / 3.14159 < 5 &&
            Q2_data > 1 &&
            Missing_data->M2() < 0.5 &&
            Missing_data->M2() > -0.5 &&
            Electron_data->P() > 2 &&
            Neutron_data->Theta() * 180. / 3.141592 > 4 &&
            Neutron_data->P() > 0.25 &&
            status_Kp_data >= 2000 &&
            status_Kp_data <= 2999 &&
            status_Km_data >= 2000 &&
            status_Km_data <= 2999 &&
            MinvKpKm_data < 1.04 &&
            MinvKpKm_data > 1.00;

        bool cut2 =
            angle_data * 180. / 3.14159 < 5 &&
            Q2_data > 1 &&
            Missing_data->M2() < 0.5 &&
            Missing_data->M2() > -0.5 &&
            Electron_data->P() > 2 &&
            Neutron_data->Theta() * 180. / 3.141592 > 4 &&
            Neutron_data->P() > 0.25 &&
            status_Kp_data >= 2000 &&
            status_Kp_data <= 2999 &&
            status_Km_data >= 2000 &&
            status_Km_data <= 2999;

        if (cut) {
            h_data->Fill(t_data_var);
        }

        if(cut2){
            double Minvdata = (*Kp + *Km).M();
            h_data_minv->Fill(Minvdata);
        }
    }

    // -------------------------
    // Fill MC
    // -------------------------
    Long64_t n_mc = t_mc->GetEntries();

    for (Long64_t i = 0; i < n_mc; i++) {
        t_mc->GetEntry(i);

        if (!Missing_mc || !Electron_mc || !Neutron_mc) continue;

        bool cut =
            angle_mc * 180. / 3.14159 < 5 &&
            Q2_mc > 1 &&
            Missing_mc->M2() < 0.5 &&
            Missing_mc->M2() > -0.5 &&
            Electron_mc->P() > 2 &&
            Neutron_mc->Theta() * 180. / 3.141592 > 4 &&
            Neutron_mc->P() > 0.25 &&
            status_Kp_mc >= 2000 &&
            status_Kp_mc <= 2999 &&
            status_Km_mc >= 2000 &&
            status_Km_mc <= 2999 &&
            MinvKpKm_mc < 1.04 &&
            MinvKpKm_mc > 1.00;


        bool cut2 =
            angle_mc * 180. / 3.14159 < 5 &&
            Q2_mc > 1 &&
            Missing_mc->M2() < 0.5 &&
            Missing_mc->M2() > -0.5 &&
            Electron_mc->P() > 2 &&
            Neutron_mc->Theta() * 180. / 3.141592 > 4 &&
            Neutron_mc->P() > 0.25 &&
            status_Kp_mc >= 2000 &&
            status_Kp_mc <= 2999 &&
            status_Km_mc >= 2000 &&
            status_Km_mc <= 2999;

        if (cut) {
            h_mc->Fill(t_missing_nucleon_mc, real_weight_SDME_effcorrection);
            h_mcs->Fill(t_missing_smeared_mc, real_weight_SDME_effcorrection);
        }

        if(cut2){

            double Minv_mc = (*Kp_mc + *Km_mc).M();
            double Minv_mc_smear = (*Kp_mc_smear + *Km_mc_smear).M();
            h_mc_minv->Fill(Minv_mc, real_weight_SDME_effcorrection);
            h_mcs_minv->Fill(Minv_mc_smear, real_weight_SDME_effcorrection);
        }
    }

    // Optional: normalize all to unit area
 

    // -------------------------
    // Draw
    // -------------------------
    TCanvas *c = new TCanvas("c", "Minv comparison", 900, 700);

    //h_data->Draw("HIST");
    h_data->SetMaximum(h_data->GetMaximum()*3.1);
    //h_mc->Draw("HIST SAME");
    //h_mcs->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.35, 0.68, 0.68, 0.88);
    leg->AddEntry(h_data, "Data: t", "l");
    leg->AddEntry(h_mc, "MC: t_missing_nucleon", "l");
    leg->AddEntry(h_mcs, "MC smeared: t_missing_smeared", "l");
    //leg->Draw();

    h_data_minv->Draw("E");
    h_mc_minv->Scale(1100.0/1700.0);
    h_mcs_minv->Scale(1100.0/1700.0);
    h_data_minv->SetMaximum(h_mc_minv->GetMaximum()*1.1);
    h_mc_minv->Draw("HIST SAME");
    //h_mcs_minv->Draw("HIST SAME");

    TLegend *leg2 = new TLegend(0.55, 0.68, 0.78, 0.88);
    leg2->AddEntry(h_data_minv, "Data RGB fall2019", "l");
    leg2->AddEntry(h_mc_minv, "MC ", "l");
    //leg2->AddEntry(h_mcs_minv, "MC smeared ", "l");
    leg2->Draw();

    cout << "integral MC : " << h_mc_minv->Integral() << endl;
    cout << "integral MC smeared : " << h_mcs_minv->Integral() << endl;

    c->SaveAs("comparison_Minv_data_mc_smeared_v1_withoutsmearing.pdf");
}
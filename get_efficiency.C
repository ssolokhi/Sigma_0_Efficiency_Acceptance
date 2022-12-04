#include <iostream>
#include <memory>

#include "TMath.h"
#include "TFile.h"
#include "TList.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraphErrors.h"

double background_fit(double *x, double *par) {
	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
}

double Gaussian_fit(double *x, double *par) {
	double argument = (x[0] - par[1])/par[2];
	return par[0]*TMath::Exp(-argument*argument/2)/(par[2]*TMath::Sqrt(2*TMath::Pi()));
}

double Gaussian_with_background_fit(double *x, double *par) {
	return Gaussian_fit(x, par) + background_fit(x, &par[3]);
}

void get_efficiency() {
	const double rapidity_cut = 0.5;
	const int n_pT_bins = 8;

	const double pT_edges[n_pT_bins+1] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 10.0};
	const double pT_bins_center[n_pT_bins] = {2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.5, 8.0};

	double Sigma_yield[n_pT_bins] = {0};
	double Sigma_mass[n_pT_bins] = {0};
	double Sigma_mass_width[n_pT_bins] = {0};

	TFile *file_data = TFile::Open("HISTOS-DT-4nov22.root", "READ");
	if (!file_data) return;
	TDirectoryFile *folder_name = dynamic_cast<TDirectoryFile*>(file_data->Get("MyTask"));
	if (!folder_name) return;
	TList *histos_data = dynamic_cast<TList*>(folder_name->Get("Histos"));
	if (!histos_data) return;

	TH2F *inv_mass_vs_pT_bins = dynamic_cast<TH2F*>(histos_data->FindObject("hLamGv0"));
	TFile *file_MC = TFile::Open("HISTOS-MC-4nov22.root", "READ");
	if (!file_MC) return;
	folder_name = dynamic_cast<TDirectoryFile*>(file_MC->Get("MyTask"));
	if (!folder_name) return;
	TList *histos_MC = dynamic_cast<TList*>(folder_name->Get("Histos"));
	if (!histos_MC) return;

	TH2F *generated_Sigma_spectrum = dynamic_cast<TH2F*>(histos_data->FindObject("hmc4piSig0Eta"));

	TH1F *inv_mass_in_pT_bin[n_pT_bins]= {0};
	const int n_bins_in_histo = 50;
	const double left_mass_histo_bound = 1.14;
	const double right_mass_histo_bound = 1.22;
	const double bin_width = (right_mass_histo_bound - left_mass_histo_bound)/n_bins_in_histo;
	const int bins_per_1GeV = 4;
	const int offset = 1;

	const double Sigma_0_mass_PDG = 1.192;
	const double deviation = 0.02;
	const int polynome_degree = 3;

	for (int i = 0; i < n_pT_bins; ++i) {
		inv_mass_in_pT_bin[i] = new TH1F(Form("inv_mass_in_pT_bin_%d",i), Form("inv_mass_in_pT_bin_%d",i), n_bins_in_histo, left_mass_histo_bound, right_mass_histo_bound);
		
		int bin_index = int(bins_per_1GeV*pT_edges[i] + offset);
		int next_bin_index = int(bins_per_1GeV*pT_edges[i+1] + offset);
		inv_mass_in_pT_bin[i] = reinterpret_cast<TH1F*>(inv_mass_vs_pT_bins->ProjectionX(Form("p_T_%d", i), bin_index, next_bin_index - 1));
	}

	for (int i = 0; i < n_pT_bins; ++i) {
		double left_fit_bound = Sigma_0_mass_PDG - deviation;
		double right_fit_bound = Sigma_0_mass_PDG + deviation;
		TF1 *background = new TF1("background", background_fit, left_fit_bound, right_fit_bound, polynome_degree + 1);
		inv_mass_in_pT_bin[i]->Fit("background", "IMEQ+", "", left_fit_bound, right_fit_bound);

		TF1 *full_fit = new TF1("full_fit", Gaussian_with_background_fit, left_fit_bound, right_fit_bound, polynome_degree + 1 + 3);
		full_fit->SetParLimits(0, 0, 10000);
		full_fit->SetParLimits(1, 1.182, 1.202);
		full_fit->SetParLimits(2, 0.001, 0.1);
		full_fit->SetLineColor(3);

		inv_mass_in_pT_bin[i]->Fit("full_fit", "IMEQ+", "", left_fit_bound, right_fit_bound);

		Sigma_yield[i] = 0;
		Sigma_mass[i] = full_fit->GetParameter(1);
		Sigma_mass_width[i] = full_fit->GetParameter(2);
		left_fit_bound = Sigma_mass[i] - 3*Sigma_mass_width[i];
		right_fit_bound = Sigma_mass[i] + 3*Sigma_mass_width[i];

		TAxis *axis = inv_mass_in_pT_bin[i]->GetXaxis();
		double n_particles_in_fit_range = inv_mass_in_pT_bin[i]->Integral(axis->FindBin(left_fit_bound), axis->FindBin(right_fit_bound));
		double full_yield = full_fit->Integral(left_fit_bound, right_fit_bound)/bin_width; // normalization factor needed since TH1 and TF1 objects are detached

		printf("Full yield from fit %f\n", full_yield);
		printf("Full yield from histo integral %f\n", n_particles_in_fit_range);
	}


	TFile *results = TFile::Open("results.root", "RECREATE");
	for (int i = 0; i < n_pT_bins; ++i) {
		inv_mass_in_pT_bin[i]->Write();
	}

	TGraphErrors *Sigma_yield_pT_dependency = new TGraphErrors(n_pT_bins, pT_bins_center, Sigma_yield);
	Sigma_yield_pT_dependency->SetTitle("#{Sigma^0} yield #{p_T} dependency");
	Sigma_yield_pT_dependency->SetMarkerStyle(20);
	Sigma_yield_pT_dependency->SetMinimum(0);
	Sigma_yield_pT_dependency->SetMaximum(5000);
	Sigma_yield_pT_dependency->Write();

	TGraphErrors *Sigma_mass_pT_dependency = new TGraphErrors(n_pT_bins, pT_bins_center, Sigma_mass);
	Sigma_mass_pT_dependency->SetTitle("#Sigma^0 mass #p_T dependency");
	Sigma_mass_pT_dependency->SetMarkerStyle(20);
	Sigma_mass_pT_dependency->SetMinimum(1.175);
	Sigma_mass_pT_dependency->SetMaximum(1.225);
	Sigma_mass_pT_dependency->Write();

	TGraphErrors *Sigma_mass_width_pT_dependency = new TGraphErrors(n_pT_bins, pT_bins_center, Sigma_mass_width);
	Sigma_mass_width_pT_dependency->SetTitle("#Sigma^0 mass width #p_T dependency");
	Sigma_mass_width_pT_dependency->SetMarkerStyle(20);
	Sigma_mass_width_pT_dependency->SetMinimum(0);
	Sigma_mass_width_pT_dependency->SetMaximum(0.05);
	Sigma_mass_width_pT_dependency->Write();

	std::cout << "Results written to output file!" << std::endl;
	results->Close();
}
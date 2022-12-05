#include <iostream>
#include <memory>

#include "TMath.h"
#include "TFile.h"
#include "TList.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

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
	const double rapidity_cut = 0.8;
	const int n_pT_bins = 8;

	const double pT_edges[n_pT_bins+1] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 10.0};
	const double pT_bins_center[n_pT_bins] = {2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.5, 8.0};

	double Sigma_yield[n_pT_bins] = {0};
	double Sigma_mass[n_pT_bins] = {0};
	double Sigma_mass_width[n_pT_bins] = {0};
	double significance[n_pT_bins] = {0};
	double efficiency[n_pT_bins] = {0};
	double Sigma_spectrum[n_pT_bins] = {0};


	TFile *file_data = TFile::Open("HISTOS-DT-21nov22.root", "READ");
	if (!file_data) return;
	TDirectoryFile *folder_name = dynamic_cast<TDirectoryFile*>(file_data->Get("MyTask"));
	if (!folder_name) return;
	TList *histos_data = dynamic_cast<TList*>(folder_name->Get("Histos"));
	if (!histos_data) return;

	TH1F *z_vertex = dynamic_cast<TH1F*>(histos_data->FindObject("hZvertex"));
	TH2F *inv_mass_vs_pT_bins = dynamic_cast<TH2F*>(histos_data->FindObject("hLamGv0"));

	TFile *file_MC = TFile::Open("HISTOS-MC-21nov22.root", "READ");
	if (!file_MC) return;
	folder_name = dynamic_cast<TDirectoryFile*>(file_MC->Get("MyTask"));
	if (!folder_name) return;
	TList *histos_MC = dynamic_cast<TList*>(folder_name->Get("Histos"));
	if (!histos_MC) return;

	TH1F *generated_Sigma_spectrum = dynamic_cast<TH1F*>(histos_MC->FindObject("hmc4piSig0Eta"));

	TH1F *inv_mass_in_pT_bin[n_pT_bins]= {0};
	const int n_bins_in_histo = 100;
	const double left_mass_histo_bound = 1.13;
	const double right_mass_histo_bound = 1.23;
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
		inv_mass_in_pT_bin[i]->Fit("background", "IMEQ", "", left_fit_bound, right_fit_bound);

		TF1 *full_fit = new TF1("full_fit", Gaussian_with_background_fit, left_fit_bound, right_fit_bound, polynome_degree + 1 + 3);
		full_fit->SetParLimits(0, 0, 10000);
		full_fit->SetParLimits(1, 1.182, 1.202);
		full_fit->SetParLimits(2, 0.001, 0.1);
		full_fit->SetParameter(3, background->GetParameter(0)*0.8);
		full_fit->SetParameter(4, background->GetParameter(1)*0.8);
		full_fit->SetParameter(5, background->GetParameter(2)*0.8);
		full_fit->SetParameter(6, background->GetParameter(3)*0.8);
		full_fit->SetLineColor(3);

		inv_mass_in_pT_bin[i]->Fit("full_fit", "IMEQ+", "", left_fit_bound, right_fit_bound);

		background->SetParameter(0, full_fit->GetParameter(3));
		background->SetParameter(1, full_fit->GetParameter(4));
		background->SetParameter(2, full_fit->GetParameter(5));
		background->SetParameter(3, full_fit->GetParameter(6));		

		Sigma_yield[i] = 0;
		Sigma_mass[i] = full_fit->GetParameter(1);
		Sigma_mass_width[i] = full_fit->GetParameter(2);
		left_fit_bound = Sigma_mass[i] - 3*Sigma_mass_width[i];
		right_fit_bound = Sigma_mass[i] + 3*Sigma_mass_width[i];

		TAxis *axis = inv_mass_in_pT_bin[i]->GetXaxis();
		double n_particles_in_fit_range = inv_mass_in_pT_bin[i]->Integral(axis->FindBin(left_fit_bound), axis->FindBin(right_fit_bound));
		double full_yield = full_fit->Integral(left_fit_bound, right_fit_bound)/bin_width; // normalization factor needed since TH1 and TF1 objects are detached
		double background_yield = background->Integral(left_fit_bound, right_fit_bound)/bin_width;
		Sigma_yield[i] = full_yield - background_yield;
		significance[i] = Sigma_yield[i]/TMath::Sqrt(full_yield);
		printf("Full yield from fit %f\n", full_yield);
		printf("Full yield from histo integral %f\n", n_particles_in_fit_range);
		printf("Sigma yield is %f\n", Sigma_yield[i]);
		printf("Significance is %f\n", significance[i]);

		int n_pT_bins_generated = generated_Sigma_spectrum->GetNbinsX();

		int bin_index = int(bins_per_1GeV*pT_edges[i] + offset);
		int next_bin_index = int(bins_per_1GeV*pT_edges[i+1] + offset);
		double min_pT = inv_mass_vs_pT_bins->GetYaxis()->GetBinLowEdge(bin_index);
		double max_pT = inv_mass_vs_pT_bins->GetYaxis()->GetBinUpEdge(next_bin_index - 1);

		int n_events_generated = 0;
		for (int j = 0; j < n_pT_bins_generated; ++j) {
			double min_pT_generated = generated_Sigma_spectrum->GetXaxis()->GetBinLowEdge(j);
			double max_pT_generated = generated_Sigma_spectrum->GetXaxis()->GetBinUpEdge(j);

			if (min_pT_generated >= min_pT && max_pT_generated <= max_pT) {
				int n_events_generated_in_bin = generated_Sigma_spectrum->GetBinContent(j);
				n_events_generated += n_events_generated_in_bin;
			}
		}

		printf("Event generated: %d\n", n_events_generated);
		efficiency[i] = Sigma_yield[i]/n_events_generated;
		double pT_bins_center_error = 0.25;
		long int n_inelastic_interactions = z_vertex->GetEntries();
		double base = efficiency[i]*2*rapidity_cut*2*pT_bins_center_error*n_inelastic_interactions*1*2;
		Sigma_spectrum[i] = Sigma_yield[i]/base;
	}

	TFile *results = TFile::Open("results.root", "RECREATE");
	for (int i = 0; i < n_pT_bins; ++i) {
		inv_mass_in_pT_bin[i]->Write();
	}

	TGraphErrors *Sigma_yield_pT_dependency = new TGraphErrors(n_pT_bins, pT_bins_center, Sigma_yield);
	Sigma_yield_pT_dependency->SetTitle("#Sigma^{0} yield p_{T} dependency");
	Sigma_yield_pT_dependency->SetMarkerStyle(20);
	Sigma_yield_pT_dependency->SetMinimum(0);
	Sigma_yield_pT_dependency->SetMaximum(2500);
	Sigma_yield_pT_dependency->Write();

	TGraphErrors *Sigma_mass_pT_dependency = new TGraphErrors(n_pT_bins, pT_bins_center, Sigma_mass);
	Sigma_mass_pT_dependency->SetTitle("#Sigma^{0} mass p_{T} dependency");
	Sigma_mass_pT_dependency->SetMarkerStyle(20);
	Sigma_mass_pT_dependency->SetMinimum(1.185);
	Sigma_mass_pT_dependency->SetMaximum(1.2);
	Sigma_mass_pT_dependency->Write();

	TGraphErrors *Sigma_mass_width_pT_dependency = new TGraphErrors(n_pT_bins, pT_bins_center, Sigma_mass_width);
	Sigma_mass_width_pT_dependency->SetTitle("#Sigma^{0} mass width p_{T} dependency");
	Sigma_mass_width_pT_dependency->SetMarkerStyle(20);
	Sigma_mass_width_pT_dependency->SetMinimum(0);
	Sigma_mass_width_pT_dependency->SetMaximum(0.01);
	Sigma_mass_width_pT_dependency->Write();

	TGraphErrors *significance_pT_dependency = new TGraphErrors(n_pT_bins, pT_bins_center, significance);
	significance_pT_dependency->SetTitle("#Sigma^{0} significance p_{T} dependency");
	significance_pT_dependency->SetMarkerStyle(20);
	significance_pT_dependency->SetMinimum(0);
	significance_pT_dependency->SetMaximum(50);
	significance_pT_dependency->Write();

	TGraphErrors *efficiency_pT_dependency = new TGraphErrors(n_pT_bins, pT_bins_center, efficiency);
	efficiency_pT_dependency->SetTitle("#Sigma^{0} efficiency p_{T} dependency");
	efficiency_pT_dependency->SetMarkerStyle(20);
	//efficiency_pT_dependency->SetMinimum(0);
	//efficiency_pT_dependency->SetMaximum(1);
	efficiency_pT_dependency->Write();

	TGraphErrors *Sigma_spectrum_pT_dependency = new TGraphErrors(n_pT_bins, pT_bins_center, Sigma_spectrum);
	Sigma_spectrum_pT_dependency->SetTitle("#Sigma^{0} spectrum p_{T} dependency");
	Sigma_spectrum_pT_dependency->SetMarkerStyle(20);
	//Sigma_spectrum_pT_dependency->SetMinimum(0);
	//Sigma_spectrum_pT_dependency->SetMaximum(100);
	Sigma_spectrum_pT_dependency->Write();

	TCanvas *dependencies = new TCanvas("dependencies", "Various p_{T} dependencies", 0, 0, 800, 1000);
	gStyle->SetOptStat(111111);
	dependencies->Divide(3,2);
	dependencies->cd(1);
	Sigma_yield_pT_dependency->Draw();
	dependencies->cd(2);
	Sigma_mass_pT_dependency->Draw();	
	dependencies->cd(3);
	Sigma_mass_width_pT_dependency->Draw();	
	dependencies->cd(4);
	significance_pT_dependency->Draw();
	dependencies->cd(5);
	efficiency_pT_dependency->Draw();	
	dependencies->cd(6);
	Sigma_spectrum_pT_dependency->Draw();

	std::cout << "Results written to output file!" << std::endl;
	results->Close();
}
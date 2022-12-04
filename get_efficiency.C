#include <iostream>
#include <memory>

#include "TMath.h"
#include "TH1F.h"
#include "TFile.h"
#include "TList.h"

double background(double *x, double *par) {
	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
}

double Voight_with_background(double *x, double *par) {
	return TMath::Voigt(*x, par[0], par[1]) + background(x, &par[3]);
}

void customize_inv_mass_histogram(TH1F *histo, const int ) {

}

void get_efficiency() {
	const double rapidity_cut = 0.5;
	const int n_pT_bins = 8;

	const double pT_edges[n_pT_bins+1] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 10.0};

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
	int n_bins_in_histo = 100;
	double left_mass_histo_bound = 1.13;
	double right_mass_histo_bound = 1.23;
	const int bins_per_1GeV = 4;
	const int offset = 1;

	for (int i = 0; i < n_pT_bins; ++i) {
		inv_mass_in_pT_bin[i] = new TH1F(Form("inv_mass_in_pT_bin_%d",i), Form("inv_mass_in_pT_bin_%d",i), n_bins_in_histo, left_mass_histo_bound, right_mass_histo_bound);
		
		int bin_index = int(bins_per_1GeV*pT_edges[i] + offset);
		int next_bin_index = int(bins_per_1GeV*pT_edges[i+1] + offset);
		inv_mass_in_pT_bin[i] = reinterpret_cast<TH1F*>(inv_mass_vs_pT_bins->ProjectionX(Form("p_T_%d", i), bin_index, next_bin_index - 1));
	}

	

	TFile *results = TFile::Open("results.root", "RECREATE");
	std::cout << "Results written to output file!" << std::endl;
	results->Close();
}
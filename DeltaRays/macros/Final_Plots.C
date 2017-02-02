#include "TSVDUnfold.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "TLegend.h"
#include "TStyle.h"
#include <iostream>
#include "NewTheoryInfo.h"

int Final_Plots(int kReg) {

	gStyle->SetOptStat(0);

	TFile *file = TFile::Open("results_1_6.root"); //Unscaled_results.root");
	
	TH1D *MC_true = (TH1D*)file->Get("True_Detected_BP");
	TH1D *MC_reco = (TH1D*)file->Get("Reconstructed_BP");
	TH1D *Reco_Data = (TH1D*)file->Get("Detect_sans_back");
	TH2D *MC_response = (TH2D*)file->Get("Adet_BP");

	TCanvas *c1 = new TCanvas("c1", "", 1200, 600);
	c1->Divide(2,1);

	TSVDUnfold *tsvdunf = new TSVDUnfold(Reco_Data, MC_reco, MC_true, MC_response);
	tsvdunf->SetNormalize(kFALSE);
	TH1D *unfresult = tsvdunf->Unfold(kReg);
	
	// Get the computed regularized covariance matrix (always corresponding to total uncertainty passed in constructor) and add uncertainties from finite MC statistics.
	//utaucov->Add( uadetcov );
	TH2D *Xtau = tsvdunf->GetXtau();
	for (int i = 1; i < 61; i++) {
		unfresult->SetBinError(i, TMath::Sqrt(Xtau->GetBinContent(i,i)));
	}

	c1->cd(1)->SetLogy();
	MC_true->SetLineColor(kRed);
	MC_true->Scale(1.99952e-09);
	unfresult->Scale(2.19252e-08);
	unfresult->SetTitle("Unfolding Result");
	unfresult->Draw();
	MC_true->Draw("same");

	c1->cd(2)->SetLogy();
	TH1D *dd = tsvdunf->GetD();
	TH1D *D_vector = new TH1D("D_vector", "", 60, 0, 60);
	for (int row = 1; row <= 60; row++) {
		D_vector->SetBinContent(row, TMath::Abs(dd->GetBinContent(row)));
	}
	D_vector->SetTitle("Abs(D_vector)");
	D_vector->Draw();

	TCanvas *c2 = new TCanvas("c2", "", 600, 600);
	c2->cd()->SetLogy();
	unfresult->Draw();
	TLegend *c2_leg = new TLegend(0.6, 0.6,0.9,0.9);
	c2_leg->AddEntry(unfresult, "Unfolded Reconstructed Data", "l");
	c2_leg->Draw();

	NewTheoryInfo *Be = new NewTheoryInfo("Be", 4, 9.0121831);
	TF1 *Theory = Be->GetTheoryCurve_ElectronEnergySpectrum(5500);

	
	TCanvas *c3 = new TCanvas("c3", "", 600, 600);
	c3->SetLogy();
	TH1D *Efficiency_BP = (TH1D*)file->Get("Efficiency_BP");
	TH1D *Final_Result = (TH1D*)unfresult->Clone();
	Final_Result->SetLineColor(kBlue);
	Final_Result->SetName("Final_Result");
	Final_Result->SetTitle("");
	Final_Result->GetYaxis()->SetTitle("Delta Ray Production Cross Section (cm^{2} g^{-1} MeV^{-1})");
	Final_Result->Divide(Efficiency_BP);
	Final_Result->Draw();
	Theory->Draw("same");

	TLegend *leg = new TLegend(0.6, 0.6, 0.9,0.9);
	leg->AddEntry(Theory, "theory", "l");
	leg->AddEntry(Final_Result, "run3", "l");
	leg->Draw();
	
}

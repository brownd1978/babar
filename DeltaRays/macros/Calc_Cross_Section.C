#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TMatrix.h"
#include "TRandom3.h"
#include "TStyle.h"
#include <iostream>

using namespace std;

void Calc_Cross_Section(TString filename, int option, float mu_len[], int N_mu[]);
void Analyze_Segment(TH1F *Detected_Delta, TH1F *Corrected_Delta, TH1F *All_Delta, TH1F *Cross_Section, TH1F *Efficiency, TH2F *Energy_Smearing, Float_t mu_len, Float_t density, Float_t lower_bound, Float_t upper_bound, Float_t N_mu, Int_t num_bins);
void Invert_Matrix(TH2F *h2, TH2F *Smearing_Correction, TMatrix &Inverted, int size, bool output_hist); 
void Sim_Efficiency(TH1F *Efficiency, TH1F *True_All, TH1F *True_Detected, Int_t num_bins);
void Calc_Acceptance(TString filename, float mu_len[], int N_mu[], float BP_len[]);
void Monte_Carlo_Cross_Section(float mu_len[], int N_mu[], float BP_mu[]);
void ToyMC(TH1F *efficiency, TH2F *energy_smearing, vector<TH1F*> hists, TH1F *True_Delta, TH1F *pull, int num_slices, double slice_size, double step_size);
void Create_Plots();


/* This function calculates the delta ray production cross
 * section from Babar data. It takes both Monte Carlo
 * and real data and can extract a cross section from either.
 * 
 *
 * There are three different data sets that are needed:
 * 
 * 3tracks -> any event with 2 muons and an electron, this 
 * 	is used to calculate the energy spread of the delta
 * 	rays
 * 
 * 2tracks -> any event where there were at least 2 muon tracks,
 * 	this is used for calculated the normalization to scale
 * 	the delta distribution into a physical cross section.
 * 
 * TFilter -> Strictly Monte Carlo data, this is any event where
 * 	there is a true delta ray created. This is used to calculate
 * 	the energy smearing and detector efficiency.
 *
 *
 * There are 3 different functions that are run. These are all run from
 * the main Calc_Cross_Section function:
 * 
 * Calc_Acceptance -> takes 2tracks data and calculates the sum of the 
 * 	path lengths for the muons that traversed the region. It 
 * 	returns a vector with the sum of path lengths, and a vector
 * 	witb the number of muons (5 entries one for each region).
 * 
 * Monte_Carlo_Cross_Section -> Reads the TFilter data to calculate
 * 	the energy smearing and detector efficiency for each region
 * 	of the detector. It also takes the True delta ray distribution
 * 	and scales it by the normalization to get the Monte Carlo 
 * 	cross section with no data corrections. This function saves 
 * 	the histograms to a file called "results.root".
 * 
 * Calc_Cross_Section -> Reads the 3tracks data. Extracts the reconstructed
 * 	delta ray distribution, then passes this to the Analyze_Segment 
 * 	function. This takes the energy smearing and efficiency histograms
 * 	and runs a matrix correction and histogram division. Finally this 
 * 	"corrected" distribution is scaled by the muon path lengths and 
 * 	material density to result in a physical cross section. This function 
 * 	is blind to whether it is getting real data or Monte Carlo
 * 	data at every point except when it saves the histograms to a file.
 * 	If it is Monte Carlo data it is saved to "results0.root"
 * 	If it is real data it is saved to "results1.root"
 * 	If it is TFilter data it is saved to "results2.root"
 * 
 * Calc_Acceptance needs to be run whenever a physical cross section is 
 * wanted, however since Monte_Carlo_Cross_Section saves to a file it 
 * does not need to be run everytime. If there were no changes made to the
 * function, computation time will be greatly reduced by disabling this
 * function. There is also a relatively simple toy monte carlo that can be 
 * run from the Monte_Carlo_Cross_Section function to test the accuracy of 
 * the correction. This is not neccesary to run every time so it is usually 
 * disabled to save on computation time. 
 *
 * The segments are:
 * Beam Pipe                                    (BP)
 * First Silicon Vertex Tracker                 (SVT1)
 * Fourth through fifth Silicon Vertex Trackers (SVT)
 * Carbon Fiber Support Tube                    (ST)
 * Drift Chamber Inner Wall                     (DCH)
 */

// These are the constants that are needed for both the real data
// and for the BaBar Monte Carlo. They are global so that both 
// functions can be run independently or together.

	Int_t const num_files = 102;

	// Densities
	Float_t const density_BP  = 1.548;
	Float_t const density_SVT = 2.192;
	Float_t const density_ST  = 1.804;
	Float_t const density_DCH = 1.552;

	// Histogram Constants
	Int_t   const num_bins_BP      = 60;
	Float_t const lower_bound_BP   = 0;//50;
	Float_t const upper_bound_BP   = 300;

	Int_t   const num_bins_SVT1    = 30;
	Float_t const lower_bound_SVT1 = 60;
	Float_t const upper_bound_SVT1 = 300;

        Int_t   const num_bins_SVT     = 25;
        Float_t const lower_bound_SVT  = 70;
        Float_t const upper_bound_SVT  = 300;

        Int_t   const num_bins_ST      = 25;
        Float_t const lower_bound_ST   = 70;
        Float_t const upper_bound_ST   = 300;

        Int_t   const num_bins_DCH     = 35;
        Float_t const lower_bound_DCH  = 60;
        Float_t const upper_bound_DCH  = 300;

	// Geometric Constants
	
		// Beam Pipe
		Float_t const BP_low_rad    = 2.5; //2.4;
		Float_t const BP_high_rad   = 2.78; //2.9;
		Float_t const BP_z_low      = -7;
		Float_t const BP_z_high     = 9;

		// First Silicon Vertex Tracker
                Float_t const SVT1_low_rad  = 3.0;
                Float_t const SVT1_high_rad = 3.85;
                Float_t const SVT1_low_x    = -0.9;
                Float_t const SVT1_high_x   = 1.6;
                Float_t const SVT1_low_y    = 3.26;
                Float_t const SVT1_high_y   = 3.425;
		Float_t const SVT1_low_U    = -1.25;
		Float_t const SVT1_high_U   = 1.25;
		Float_t const SVT1_z_low    = -17;
		Float_t const SVT1_z_high   = 17;

		// Second Silicon Vertex Tracker
                Float_t const SVT2_low_rad  = 3.85;
                Float_t const SVT2_high_rad = 4.8;
                Float_t const SVT2_low_x    = -1;
                Float_t const SVT2_high_x   = 1.7;
                Float_t const SVT2_low_y    = 3.9;
                Float_t const SVT2_high_y   = 4.05;
		Float_t const SVT2_low_U    = -1.35;
		Float_t const SVT2_high_U   = 1.35;
                Float_t const SVT2_z_low    = -17;
                Float_t const SVT2_z_high   = 17;

		// Third Silicon Vertex Tracker
                Float_t const SVT3_low_rad  = 5.3;
                Float_t const SVT3_high_rad = 7.2;
                Float_t const SVT3_low_x    = -1.2;
                Float_t const SVT3_high_x   = 2.2;
                Float_t const SVT3_low_y    = 5.775;
                Float_t const SVT3_high_y   = 5.925;
		Float_t const SVT3_low_U    = -1.7;
		Float_t const SVT3_high_U   = 1.7;
                Float_t const SVT3_z_low    = -17;
                Float_t const SVT3_z_high   = 17;

		// Fourth and fifth SVT's are really 
		// complicated,
		// For now these are place holder variables
                Float_t const SVT4_low_rad  = 11.9;
                Float_t const SVT4_high_rad = 13.5;
                //Float_t const SVT4_low_x    = -1.2;
                //Float_t const SVT4_high_x   = 2.2;
                //Float_t const SVT4_low_y    = 5.775;
                //Float_t const SVT4_high_y   = 5.925;
		Float_t const SVT4_low_U    = -1.7;
		Float_t const SVT4_high_U   = 1.7;
                Float_t const SVT4_z_low    = -17;
                Float_t const SVT4_z_high   = 17;

                Float_t const SVT5_low_rad  = 13.9;
                Float_t const SVT5_high_rad = 15.55;
                //Float_t const SVT5_low_x    = -1.2;
                //Float_t const SVT5_high_x   = 2.2;
                //Float_t const SVT5_low_y    = 5.775;
                //Float_t const SVT5_high_y   = 5.925;
                //Float_t const SVT5_z_low    = 0;
                //Float_t const SVT5_z_high   = 0;


		// Support Tube
                Float_t const ST_low_rad   = 21.6;
                Float_t const ST_high_rad  = 21.95;
		Float_t const ST_z_low     = -44.3;
		Float_t const ST_z_high    = 74.3;

		// Beryllium Drift Chamber
		Float_t const DCH_low_rad  = 23.7;
		Float_t const DCH_high_rad = 23.95;
		Float_t const DCH_z_low    = -44.3;
		Float_t const DCH_z_high   = 74.3;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

void Calc_Cross_Section() {
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	float MC_mu_len[5], Real_mu_len[5], MC_BP_len[5], Real_BP_len[5];
	int MC_N_mu[5], Real_N_mu[5];

	// Both TFilter and MC data have same normalization so it is only calculated once
	Calc_Acceptance("/data/HD4/babar/mupair/KK2F/mupair_kk2f_split-BTM_2tracks-", MC_mu_len, MC_N_mu, MC_BP_len);

	// Calculates the efficiency and energy smearing and saves it to a file
	// if no changes are made to the code then this does not need to be run
	// every time.
	Monte_Carlo_Cross_Section(MC_mu_len, MC_N_mu, MC_BP_len);
	
	  // Run on TFilter data (looking for discrepencies)
	  // Calc_Cross_Section("/data/HD4/babar/mupair/KK2F/mupair_kk2f_split-BTM_TFilter-", 2, MC_mu_len, MC_N_mu);

	// run analysis on Monte Carlo data (formatted like real data)
	Calc_Cross_Section("/data/HD4/babar/mupair/KK2F/mupair_kk2f_split-BTM_3tracks-", 0, MC_mu_len, MC_N_mu);


	// run the whole analysis on the actual data from BaBar
	Calc_Acceptance("/data/HD4/babar/mupair/run3/AES_run3-BTM_2tracks-", Real_mu_len, Real_N_mu, Real_BP_len);
	Calc_Cross_Section("/data/HD4/babar/mupair/run3/AES_run3-BTM_3tracks-", 1, Real_mu_len, Real_N_mu);	

	Create_Plots();
	
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


// Reads the raw data
void Calc_Cross_Section(TString filename, int option, float mu_len[], int N_mu[]) {

        Float_t step_size_BP = (upper_bound_BP - lower_bound_BP)/num_bins_BP;
        Float_t step_size_SVT1 = (upper_bound_SVT1 - lower_bound_SVT1)/num_bins_SVT1;
        Float_t step_size_SVT = (upper_bound_SVT - lower_bound_SVT)/num_bins_SVT;
        Float_t step_size_ST = (upper_bound_ST - lower_bound_ST)/num_bins_ST;
        Float_t step_size_DCH = (upper_bound_DCH - lower_bound_DCH)/num_bins_DCH;

	// Read efficiency and energy smearing from file:
                TFile *file = TFile::Open("results.root");
                TH1F *Efficiency_BP = (TH1F*)file->Get("Efficiency_BP");
                TH2F *Energy_Smearing_BP = (TH2F*)file->Get("Energy_Smearing_BP");
                Efficiency_BP->SetDirectory(0);
                Energy_Smearing_BP->SetDirectory(0);

        /*        TH1F *Efficiency_SVT1 = (TH1F*)file->Get("Efficiency_SVT1");
                TH2F *Energy_Smearing_SVT1 = (TH2F*)file->Get("Energy_Smearing_SVT1");
                Efficiency_SVT1->SetDirectory(0);
                Energy_Smearing_SVT1->SetDirectory(0);

                TH1F *Efficiency_SVT = (TH1F*)file->Get("Efficiency_SVT");
                TH2F *Energy_Smearing_SVT = (TH2F*)file->Get("Energy_Smearing_SVT");
                Efficiency_SVT->SetDirectory(0);
                Energy_Smearing_SVT->SetDirectory(0);

                TH1F *Efficiency_ST = (TH1F*)file->Get("Efficiency_ST");
                TH2F *Energy_Smearing_ST = (TH2F*)file->Get("Energy_Smearing_ST");
                Efficiency_ST->SetDirectory(0);
                Energy_Smearing_ST->SetDirectory(0);

                TH1F *Efficiency_DCH = (TH1F*)file->Get("Efficiency_DCH");
                TH2F *Energy_Smearing_DCH = (TH2F*)file->Get("Energy_Smearing_DCH");
                Efficiency_DCH->SetDirectory(0);
                Energy_Smearing_DCH->SetDirectory(0);*/

                file->Close();



	TChain *ch = new TChain("ntp667");
	for (int i = 1; i <= num_files; i++) {
		ch->AddFile(TString::Format(filename +"%d" + ".root", i));
	}

	// Initialize the raw data TTree (lots of variable definitions)
        int nevents = ch->GetEntries();
        Int_t ne, nmu, nDelta;
        Int_t eLund[100], muTrkIdx[100], eTrkIdx[100];
        Int_t Deltad1Idx[100], Deltad2Idx[100], Deltad1Lund[100];
        Int_t TRKnDchXY[100], TRKnDchZ[100], TRKnSvtXY[100], TRKnSvtZ[100];
        Float_t pTotalScalar, mucosth[100];
        Float_t eenergy[100], muenergy[100];
        Float_t DeltaDirDot[100], DeltaDOCA[100];
        Float_t DeltaePOCA0[100], DeltaePOCA1[100], DeltaePOCA2[100];
	Float_t eNeutralMag, TRKEMCecal[100];
        Float_t muBP0[100], muBP1[100], muBP2[100], muBPndot[100];
/*      Float_t muDCH0[100],  muDCH1[100],  muDCH2[100],  muDCHndot[100];
        Float_t muST0[100],   muST1[100],   muST2[100],   muSTndot[100];
        Float_t muSVT10[100], muSVT11[100], muSVT12[100], muSVTndot1[100], muSVTU1[100];
        Float_t muSVT20[100], muSVT21[100], muSVT22[100], muSVTndot2[100], muSVTU2[100];
        Float_t muSVT30[100], muSVT31[100], muSVT32[100], muSVTndot3[100], muSVTU3[100];
        Float_t muSVT40[100], muSVT41[100], muSVT42[100], muSVTndot4[100], muSVTU4[100];
        Float_t muSVT50[100], muSVT51[100], muSVT52[100], muSVTndot5[100], muSVTU5[100];
	Float_t ep3[100];*/

        ch->SetBranchStatus("*", 0);
        ch->SetBranchAddress("ne", &ne);
        ch->SetBranchAddress("nmu", &nmu);
        ch->SetBranchAddress("nDelta", &nDelta);
        ch->SetBranchAddress("eLund", eLund);
        ch->SetBranchAddress("muTrkIdx", muTrkIdx);
        ch->SetBranchAddress("eTrkIdx", eTrkIdx);
        ch->SetBranchAddress("Deltad1Idx", Deltad1Idx);
        ch->SetBranchAddress("Deltad2Idx", Deltad2Idx);
        ch->SetBranchAddress("Deltad1Lund", Deltad1Lund);
        ch->SetBranchAddress("TRKnDchXY", TRKnDchXY);
        ch->SetBranchAddress("TRKnDchZ", TRKnDchZ);
        ch->SetBranchAddress("TRKnSvtXY", TRKnSvtXY);
        ch->SetBranchAddress("TRKnSvtZ", TRKnSvtZ);
        ch->SetBranchAddress("eenergy", eenergy);
        ch->SetBranchAddress("muenergy", muenergy);
        ch->SetBranchAddress("DeltaDirDot", DeltaDirDot);
        ch->SetBranchAddress("DeltaDOCA", DeltaDOCA);
        ch->SetBranchAddress("pTotalScalar", &pTotalScalar);
        ch->SetBranchAddress("DeltaePOCA0", DeltaePOCA0);
        ch->SetBranchAddress("DeltaePOCA1", DeltaePOCA1);
        ch->SetBranchAddress("DeltaePOCA2", DeltaePOCA2);
	ch->SetBranchAddress("eNeutralMag", &eNeutralMag);
	ch->SetBranchAddress("TRKEMCecal", TRKEMCecal);
        ch->SetBranchAddress("mucosth", mucosth);
        ch->SetBranchAddress("muBP0", muBP0);
        ch->SetBranchAddress("muBP1", muBP1);
        ch->SetBranchAddress("muBP2", muBP2);
	ch->SetBranchAddress("muBPndot", muBPndot);
/*	ch->SetBranchAddress("muDCH0", muDCH0);
        ch->SetBranchAddress("muDCH1", muDCH1);
        ch->SetBranchAddress("muDCH2", muDCH2);
	ch->SetBranchAddress("muDCHndot", muDCHndot);
	ch->SetBranchAddress("muST0", muST0);
        ch->SetBranchAddress("muST1", muST1);
        ch->SetBranchAddress("muST2", muST2);
	ch->SetBranchAddress("muSTndot", muSTndot);
        ch->SetBranchAddress("muSVT10", muSVT10);
        ch->SetBranchAddress("muSVT11", muSVT11);
        ch->SetBranchAddress("muSVT12", muSVT12);
	ch->SetBranchAddress("muSVTndot1", muSVTndot1);
	ch->SetBranchAddress("muSVTU1", muSVTU1)'
        ch->SetBranchAddress("muSVT20", muSVT20);
        ch->SetBranchAddress("muSVT21", muSVT21);
        ch->SetBranchAddress("muSVT22", muSVT22);
        ch->SetBranchAddress("muSVTndot2", muSVTndot2);
	ch->SetBranchAddress("muSVTU2", muSVTU2);
        ch->SetBranchAddress("muSVT30", muSVT30);
        ch->SetBranchAddress("muSVT31", muSVT31);
        ch->SetBranchAddress("muSVT32", muSVT32);
        ch->SetBranchAddress("muSVTndot3", muSVTndot3);
	ch->SetBranchAddress("muSVTU3", muSVTU3);
        ch->SetBranchAddress("muSVT40", muSVT40);
        ch->SetBranchAddress("muSVT41", muSVT41);
        ch->SetBranchAddress("muSVT42", muSVT42);
        ch->SetBranchAddress("muSVTndot4", muSVTndot4);
	ch->SetBranchAddress("muSVTU4", muSVTU4);
        ch->SetBranchAddress("muSVT50", muSVT50);
        ch->SetBranchAddress("muSVT51", muSVT51);
        ch->SetBranchAddress("muSVT52", muSVT52);
        ch->SetBranchAddress("muSVTndot5", muSVTndot5);
	ch->SetBranchAddress("muSVTU5", muSVTU5);
*/
	// ch->SetBranchAddress("ep3", ep3);
	// of the required Histograms and region specific variables are defined here
	// Efficiency and energy smearing is read from a file

		// Beam Pipe
		int N_mu_BP = 0;
		Float_t mu_len_BP = 0;
		TH1F *Detected_Delta_BP = new TH1F("Detected_Delta_BP", "Raw Detected Deltas;Reco Energy (MeV)", 60, 0, upper_bound_BP);
		TH1F *Detect_sans_back = new TH1F("Detect_sans_back", "Detected Deltas minus Gamma Events;Reco Energy (MeV) ", num_bins_BP, 0, upper_bound_BP);
		TH1F *Corrected_Delta_BP = new TH1F("Corrected_Delta_BP", "Energy Smearing Corrected Deltas;True Energy (MeV)", 60, 0, upper_bound_BP);
		TH1F *All_Delta_BP = new TH1F("All_Delta_BP", "Efficiency Corrected Deltas;True Energy (MeV)", 60, 0, upper_bound_BP);
		TH1F *Cross_Section_BP = new TH1F("Cross_Section_BP", "Delta Ray Production Cross Section;True Energy (MeV);Cross Section cm^{2} g^{-1} MeV^{-1}", 60, 0, upper_bound_BP);
		TH1F *Gamma_BP = new TH1F("Gamma_BP", "Gamma Conversion Positrons;Reco Energy (MeV)", 60, 0, upper_bound_BP);

        	Detected_Delta_BP->Sumw2();
		Detect_sans_back->Sumw2();
		Corrected_Delta_BP->Sumw2();
		All_Delta_BP->Sumw2();
		Cross_Section_BP->Sumw2();
		Gamma_BP->Sumw2();
/*
		// First Silicon Vertex Tracker
        	int N_mu_SVT1 = 0;
        	Float_t mu_len_SVT1 = 0;
       	 	TH1F *Detected_Delta_SVT1 = new TH1F("Detected_Delta_SVT1", "", num_bins_SVT1 + 1, lower_bound_SVT1 - step_size_SVT1, upper_bound_SVT1);
        	TH1F *Corrected_Delta_SVT1 = new TH1F("Corrected_Delta_SVT1", "", num_bins_SVT1, lower_bound_SVT1, upper_bound_SVT1);
        	TH1F *All_Delta_SVT1 = new TH1F("All_Delta_SVT1", "", num_bins_SVT1, lower_bound_SVT1, upper_bound_SVT1);
        	TH1F *Cross_Section_SVT1 = new TH1F("Cross_Section_SVT1", "", num_bins_SVT1, lower_bound_SVT1, upper_bound_SVT1);

        	Detected_Delta_SVT1->Sumw2();
        	Corrected_Delta_SVT1->Sumw2();
        	All_Delta_SVT1->Sumw2();
        	Cross_Section_SVT1->Sumw2();


		// Second through fifth Silicon Vertex Trackers
        	int N_mu_SVT = 0;
        	Float_t mu_len_SVT = 0;
        	TH1F *Detected_Delta_SVT = new TH1F("Detected_Delta_SVT", "", num_bins_SVT + 1, lower_bound_SVT - step_size_SVT, upper_bound_SVT);
        	TH1F *Corrected_Delta_SVT = new TH1F("Corrected_Delta_SVT", "", num_bins_SVT, lower_bound_SVT, upper_bound_SVT);
        	TH1F *All_Delta_SVT = new TH1F("All_Delta_SVT", "", num_bins_SVT, lower_bound_SVT, upper_bound_SVT);
        	TH1F *Cross_Section_SVT = new TH1F("Cross_Section_SVT", "", num_bins_SVT, lower_bound_SVT, upper_bound_SVT);

        	Detected_Delta_SVT->Sumw2();
        	Corrected_Delta_SVT->Sumw2();
        	All_Delta_SVT->Sumw2();
        	Cross_Section_SVT->Sumw2();


		// Carbon Fiber Support Tube
        	int N_mu_ST = 0;
        	Float_t mu_len_ST = 0;
        	TH1F *Detected_Delta_ST = new TH1F("Detected_Delta_ST", "", num_bins_ST + 1, lower_bound_ST - step_size_ST, upper_bound_ST);
        	TH1F *Corrected_Delta_ST = new TH1F("Corrected_Delta_ST", "", num_bins_ST, lower_bound_ST, upper_bound_ST);
        	TH1F *All_Delta_ST = new TH1F("All_Delta_ST", "", num_bins_ST, lower_bound_ST, upper_bound_ST);
        	TH1F *Cross_Section_ST = new TH1F("Cross_Section_ST", "", num_bins_ST, lower_bound_ST, upper_bound_ST);
	
        	Detected_Delta_ST->Sumw2();
        	Corrected_Delta_ST->Sumw2();
        	All_Delta_ST->Sumw2();
        	Cross_Section_ST->Sumw2();


		// Drift Chamber Inner Wall
        	nt N_mu_DCH = 0;
        	Float_t mu_len_DCH = 0;
        	TH1F *Detected_Delta_DCH = new TH1F("Detected_Delta_DCH", "", num_bins_DCH + 1, lower_bound_DCH - step_size_DCH, upper_bound_DCH);
        	TH1F *Corrected_Delta_DCH = new TH1F("Corrected_Delta_DCH", "", num_bins_DCH, lower_bound_DCH, upper_bound_DCH);
        	TH1F *All_Delta_DCH = new TH1F("All_Delta_DCH", "", num_bins_DCH, lower_bound_DCH, upper_bound_DCH);
        	TH1F *Cross_Section_DCH = new TH1F("Cross_Section_DCH", "", num_bins_DCH, lower_bound_DCH, upper_bound_DCH);
	
	        Detected_Delta_DCH->Sumw2();
	        Corrected_Delta_DCH->Sumw2();
	        All_Delta_DCH->Sumw2();
	        Cross_Section_DCH->Sumw2();
*/

	/* Analyze the real data segment by segment. 
 	 * first good muon track is found
 	 * then path length is recorded
 	 * then delta rays from corresponding track are found
 	 * delta ray spectrum is recorded
 	 */

	cout << endl;
	
	if (option == 0) cout << "Reading Monte Carlo Data" << endl;
	if (option == 1) cout << "Reading Real Data" << endl;
	if (option == 2) cout << "Reading TFilter Data" << endl;
	for (int i = 0; i < nevents; i++) {
		ch->GetEntry(i);
		if (i % (nevents/10) == 0) cout << i << "/" << nevents << endl;
		
                if (nmu == 2) {
                        if (TMath::Abs(pTotalScalar - 11.75) <= 0.75 && eNeutralMag < 1.0) {
                                for (int j = 0; j < nmu; j++) {
	                        if (TRKnDchXY[muTrkIdx[j]] >= 15 && TRKnDchZ[muTrkIdx[j]] >= 10) {
					if (muenergy[j] > 3.5 && muenergy[j] < 8 && TRKEMCecal[muTrkIdx[j]] < 0.5) {
						if (true) { //TRKnSvtXY[muTrkIdx[j]] >= 4 && TRKnSvtZ[muTrkIdx[j]] >= 4) {
							float mu_rad = TMath::Sqrt(TMath::Power(muBP0[j], 2) + TMath::Power(muBP1[j], 2));
							// Did muon pass through Beam Pipe
							if (mu_rad > BP_low_rad && mu_rad < BP_high_rad) {
                                                                if (muBP2[j] > BP_z_low && muBP2[j] < BP_z_high && muBPndot[j] > 0.5) {

									// Was the delta ray seen
									for (int l = 0; l < nDelta; l++) {
										if (Deltad2Idx[l] = j) {
											if (DeltaDirDot[l] > 0.98 && DeltaDOCA[l] < 0.1) {
												Float_t Delta_Rad = TMath::Sqrt(TMath::Power(DeltaePOCA0[l], 2) + TMath::Power(DeltaePOCA1[l], 2));
												if (Delta_Rad > BP_low_rad && Delta_Rad < BP_high_rad) {
													if (DeltaePOCA2[l] > BP_z_low && DeltaePOCA2[l] < BP_z_high) {
														if (Deltad1Lund[l] == 11) {
															Detected_Delta_BP->Fill(1000*eenergy[Deltad1Idx[l]]);
															
														}
														if (Deltad1Lund[l] == -11) {
															Gamma_BP->Fill(1000*eenergy[Deltad1Idx[l]]);
														}
													}
												}
											}
										}
									}
								}
							}
/*
							//First Silicon Vertex Tracker
                                                        if (TMath::Abs(TMath::Sqrt(TMath::Power(muSVT10[j], 2) + TMath::Power(muSVT11[j], 2)) - 3.5) < 1) {
                                                                if (muSVT12[j] > SVT1_z_low && muSVT12[j] < SVT1_z_high) {

                                                                        //Was the delta ray seen
                                                                        for (int l = 0; l < nDelta; l++) {
                                                                                if (Deltad2Idx[l] = j) {
                                                                                        if (Deltad1Lund[l] == 11 && DeltaDirDot[l] > 0.98) {
                                                                                                Float_t Delta_Rad = TMath::Sqrt(TMath::Power(DeltaePOCA0[l], 2) + TMath::Power(DeltaePOCA1[l], 2));
                                                                                                if (Delta_Rad > SVT1_low_rad && Delta_Rad < SVT1_high_rad) {

													if((TRKnSvtXY[eTrkIdx[Deltad1Idx[l]]] > 3) && (TRKnSvtZ[eTrkIdx[Deltad1Idx[l]]] > 3)){
        				                                                                        Float_t elec_trans_new = TMath::Sqrt(TMath::Power(DeltaePOCA0[l] + 0.01, 2) + TMath::Power(DeltaePOCA1[l] - 0.013, 2));
                				                                                                Float_t phi = TMath::ATan2(DeltaePOCA1[j] - 0.013, DeltaePOCA0[l] + 0.01) + 2*TMath::Pi() + 0.42705;
                        				                                                        Float_t new_phi = fmodf(TMath::Abs(phi), TMath::Pi()/3.0);
                                		        		                                        Float_t new_x = elec_trans_new*TMath::Cos(new_phi + TMath::Pi()/3.0);
                                        		                		                        Float_t new_y = elec_trans_new*TMath::Sin(new_phi + TMath::Pi()/3.0);

		                                                		                                if ((new_x >= SVT1_low_x) && (new_x <= SVT1_high_x) && (new_y >= SVT1_low_y) && (new_y <= SVT1_high_y)) {
                		               	                                                                        if (DeltaePOCA2[l] > SVT1_z_low && DeltaePOCA2[l] < SVT1_z_high) {
                                		                                                                                Detected_Delta_SVT1->Fill(1000*eenergy[l]);
															}
														}
                                                                                                        }
                                                                                                }
                                                                                        }
                                                                                }
                                                                        }
                                                                }
                                                        }

							//Second through fifth Vertex Trackers
                                                        if (TMath::Abs(TMath::Sqrt(TMath::Power(muSVT20[j], 2) + TMath::Power(muSVT21[j], 2)) - 2.6) < .2) {
                                                                if (muSVT22[j] > SVT2_z_low && muSVT22[j] < SVT2_z_high) {

                                                                        //Was the delta ray seen
                                                                        for (int l = 0; l < nDelta; l++) {
                                                                                if (Deltad2Idx[l] = j) {
                                                                                        if (Deltad1Lund[l] == 11 && DeltaDirDot[l] > 0.98) {
                                                                                                Float_t Delta_Rad = TMath::Sqrt(TMath::Power(DeltaePOCA0[l], 2) + TMath::Power(DeltaePOCA1[l], 2));
                                                                                                if (Delta_Rad > SVT2_low_rad && Delta_Rad < SVT2_high_rad) {
                                                                                                        if((TRKnSvtXY[eTrkIdx[Deltad1Idx[l]]] > 3) && (TRKnSvtZ[eTrkIdx[Deltad1Idx[l]]] > 3)){
                                                                                                                Float_t elec_trans_new = TMath::Sqrt(TMath::Power(DeltaePOCA0[l] + 0.01, 2) + TMath::Power(DeltaePOCA1[l] - 0.013, 2));
                                                                                                                Float_t phi = TMath::ATan2(DeltaePOCA1[j] - 0.013, DeltaePOCA0[l] + 0.01) + 2*TMath::Pi() + 0.42705;
                                                                                                                Float_t new_phi = fmodf(TMath::Abs(phi), TMath::Pi()/3.0);
                                                                                                                Float_t new_x = elec_trans_new*TMath::Cos(new_phi + TMath::Pi()/3.0);
                                                                                                                Float_t new_y = elec_trans_new*TMath::Sin(new_phi + TMath::Pi()/3.0);

                                                                                                                if ((new_x >= SVT2_low_x) && (new_x <= SVT2_high_x) && (new_y >= SVT2_low_y) && (new_y <= SVT2_high_y)) {
                                                                                                                        if (DeltaePOCA2[l] > SVT2_z_low && DeltaePOCA2[l] < SVT2_z_high) {
                                                                                                                                Detected_Delta_SVT->Fill(1000*eenergy[l]);
                                                                                                                        }
                                                                                                                }
                                                                                                        }
                                                                                                }
                                                                                        }
                                                                                }
                                                                        }
                                                                }
                                                        }
                                                        if (TMath::Abs(TMath::Sqrt(TMath::Power(muSVT30[j], 2) + TMath::Power(muSVT31[j], 2)) - 2.6) < .2) {
                                                                if (muSVT32[j] > SVT3_z_low && muSVT32[j] < SVT3_z_high) {

                                                                        //Was the delta ray seen
                                                                        for (int l = 0; l < nDelta; l++) {
                                                                                if (Deltad2Idx[l] = j) {
                                                                                        if (Deltad1Lund[l] == 11 && DeltaDirDot[l] > 0.98) {
                                                                                                Float_t Delta_Rad = TMath::Sqrt(TMath::Power(DeltaePOCA0[l], 2) + TMath::Power(DeltaePOCA1[l], 2));
                                                                                                if (Delta_Rad > SVT3_low_rad && Delta_Rad < SVT3_high_rad) {
                                                                                                        if((TRKnSvtXY[eTrkIdx[Deltad1Idx[l]]] > 3) && (TRKnSvtZ[eTrkIdx[Deltad1Idx[l]]] > 3)){
                                                                                                                Float_t elec_trans_new = TMath::Sqrt(TMath::Power(DeltaePOCA0[l] + 0.01, 2) + TMath::Power(DeltaePOCA1[l] - 0.013, 2));
                                                                                                                Float_t phi = TMath::ATan2(DeltaePOCA1[j] - 0.013, DeltaePOCA0[l] + 0.01) + 2*TMath::Pi() + 0.42705;
                                                                                                                Float_t new_phi = fmodf(TMath::Abs(phi), TMath::Pi()/3.0);
                                                                                                                Float_t new_x = elec_trans_new*TMath::Cos(new_phi + TMath::Pi()/3.0);
                                                                                                                Float_t new_y = elec_trans_new*TMath::Sin(new_phi + TMath::Pi()/3.0);

                                                                                                                if ((new_x >= SVT3_low_x) && (new_x <= SVT3_high_x) && (new_y >= SVT3_low_y) && (new_y <= SVT3_high_y)) {
                                                                                                                        if (DeltaePOCA2[l] > SVT3_z_low && DeltaePOCA2[l] < SVT3_z_high) {
                                                                                                                                Detected_Delta_SVT->Fill(1000*eenergy[l]);
                                                                                                                        }
                                                                                                                }
                                                                                                        }
                                                                                                }
                                                                                        }
                                                                                }
                                                                        }
                                                                }
                                                        } */ /*
                                                        if (TMath::Abs(TMath::Sqrt(TMath::Power(muSVT40[j], 2) + TMath::Power(muSVT41[j], 2)) - 2.6) < .2) {
                                                                if (muSVT42[j] > SVT4_z_low && muSVT42[j] < SVT4_z_high) {

                                                                        //Was the delta ray seen
                                                                        for (int l = 0; l < nDelta; l++) {
                                                                                if (Deltad2Idx[l] = j) {
                                                                                        if (Deltad1Lund[l] == 11 && DeltaDirDot[l] > 0.98) {
                                                                                                Float_t Delta_Rad = TMath::Sqrt(TMath::Power(DeltaePOCA0[l], 2) + TMath::Power(DeltaePOCA1[l], 2));
                                                                                                if (Delta_Rad > SVT4_low_rad && Delta_Rad < SVT4_high_rad) {
                                                                                                        if (DeltaePOCA2[l] > -7 && DeltaePOCA2[l] < 9) {
                                                                                                                Detected_Delta_SVT->Fill(1000*eenergy[l]);
                                                                                                        }
                                                                                                }
                                                                                        }
                                                                                }
									}
								}
							}
                                                        if (TMath::Abs(TMath::Sqrt(TMath::Power(muSVT50[j], 2) + TMath::Power(muSVT51[j], 2)) - 2.6) < .2) {
                                                                if (muSVT52[j] > SVT5_z_low && muSVT52[j] < SVT5_z_high) {

                                                                        //Was the delta ray seen
                                                                        for (int l = 0; l < nDelta; l++) {
                                                                                if (Deltad2Idx[l] = j) {
                                                                                        if (Deltad1Lund[l] == 11 && DeltaDirDot[l] > 0.98) {
                                                                                                Float_t Delta_Rad = TMath::Sqrt(TMath::Power(DeltaePOCA0[l], 2) + TMath::Power(DeltaePOCA1[l], 2));
                                                                                                if (Delta_Rad > SVT5_low_rad && Delta_Rad < SVT5_high_rad) {
                                                                                                        if (DeltaePOCA2[l] > -7 && DeltaePOCA2[l] < 9) {
                                                                                                                Detected_Delta_SVT->Fill(1000*eenergy[l]);
                                                                                                        }
                                                                                                }
                                                                                        }
                                                                                }
                                                                        }
                                                                }
                                                        }
*/
							//Carbon Fiber Support Tube
                                                        /*if (TMath::Abs(TMath::Sqrt(TMath::Power(muST0[j], 2) + TMath::Power(muST1[j], 2)) ) > 1) {
                                                                if (muST2[j] > ST_z_low && muST2[j] < ST_z_high) {

                                                                        //Was the delta ray seen
                                                                        for (int l = 0; l < nDelta; l++) {
                                                                                if (Deltad2Idx[l] = j) {
                                                                                        if (Deltad1Lund[l] == 11 && DeltaDirDot[l] > 0.98) {
                                                                                                Float_t Delta_Rad = TMath::Sqrt(TMath::Power(DeltaePOCA0[l], 2) + TMath::Power(DeltaePOCA1[l], 2));
                                                                                                if (Delta_Rad > ST_low_rad && Delta_Rad < ST_high_rad) {
                                                                                                        if (DeltaePOCA2[l] > ST_z_low && DeltaePOCA2[l] < ST_z_high) {
                                                                                                                Detected_Delta_ST->Fill(1000*eenergy[l]);
                                                                                                        }
                                                                                                }
                                                                                        }
                                                                                *}
                                                                        }
                                                                }
                                                        }

							//Drift Chamber inner wall
							if (TMath::Sqrt(TMath::Power(muDCH0[j], 2) + TMath::Power(muDCH1[j], 2)) > 1) {
                                                                if (muDCH2[j] > DCH_z_low && muDCH2[j] < DCH_z_high) {

                                                                        //Was the delta ray seen
                                                                        for (int l = 0; l < nDelta; l++) {
                                                                                if (Deltad2Idx[l] = j) {
                                                                                        if (Deltad1Lund[l] == 11 && DeltaDirDot[l] > 0.98) {
                                                                                                Float_t Delta_Rad = TMath::Sqrt(TMath::Power(DeltaePOCA0[l], 2) + TMath::Power(DeltaePOCA1[l], 2));
                                                                                                if (Delta_Rad > DCH_low_rad && Delta_Rad < DCH_high_rad) {
                                                                                                        if (DeltaePOCA2[l] > DCH_z_low && DeltaePOCA2[l] < DCH_z_high) {
                                                                                                                Detected_Delta_DCH->Fill(1000*eenergy[l]);
                                                                                                        }
                                                                                                }
                                                                                        }
                                                                                }
                                                                        }
                                                                }
                                                        }*/
						}
					}
				}
			}
		}}
	}


	// Calculate cross section
	
	for (int i = 1; i < 61; i ++) {
		Detect_sans_back->SetBinContent(i, Detected_Delta_BP->GetBinContent(i));
		Detect_sans_back->SetBinError(i, Detected_Delta_BP->GetBinError(i));
	}
	//Detect_sans_back = (TH1F*)Detected_Delta_BP->Clone(); //"Detect_sans_back");
	Detect_sans_back->Add(Gamma_BP, -1);
	//Analyze_Segment(Detect_sans_back, Corrected_Delta_BP, All_Delta_BP, Cross_Section_BP, Efficiency_BP, Energy_Smearing_BP, mu_len[0], density_BP, lower_bound_BP, upper_bound_BP, N_mu[0], num_bins_BP);

        //Analyze_Segment(Detected_Delta_SVT1, Corrected_Delta_SVT1, All_Delta_SVT1, Cross_Section_SVT1, Efficiency_SVT1, Energy_Smearing_SVT1, mu_len_SVT1, density_SVT, lower_bound_SVT1, upper_bound_SVT1, N_mu_SVT1, num_bins_SVT1);

        //Analyze_Segment(Detected_Delta_SVT, Corrected_Delta_SVT, All_Delta_SVT, Cross_Section_SVT, Efficiency_SVT, Energy_Smearing_SVT, mu_len_SVT, density_SVT, lower_bound_SVT, upper_bound_SVT, N_mu_SVT, num_bins_SVT);

        //Analyze_Segment(Detected_Delta_ST, Corrected_Delta_ST, All_Delta_ST, Cross_Section_ST, Efficiency_ST, Energy_Smearing_ST, mu_len_ST, density_ST, lower_bound_ST, upper_bound_ST, N_mu_ST, num_bins_ST);

        //Analyze_Segment(Detected_Delta_DCH, Corrected_Delta_DCH, All_Delta_DCH, Cross_Section_DCH, Efficiency_DCH, Energy_Smearing_DCH, mu_len_DCH, density_DCH, lower_bound_DCH, upper_bound_DCH, N_mu_DCH, num_bins_DCH);


	// Normalize all of the histograms so they have a non-arbitrary y axis
	Float_t norm_BP = 1./(mu_len[0]*density_BP*step_size_BP);
	std::cout << "scale : " << norm_BP << std::endl;
	//Detected_Delta_BP->Scale(norm_BP);
	//Corrected_Delta_BP->Scale(norm_BP);
	//All_Delta_BP->Scale(norm_BP);
	//Gamma_BP->Scale(norm_BP);


	// Save histograms to file, change name based on data type
	if (option == 0) {
                Detected_Delta_BP->SetName("MC_Detected_Delta_BP"); 
        	Detect_sans_back->SetName("MC_Detect_sans_back");
	        Corrected_Delta_BP->SetName("MC_Corrected_Delta_BP");
                All_Delta_BP->SetName("MC_All_Delta_BP");
                Cross_Section_BP->SetName("MC_Cross_Section_BP");
                Gamma_BP->SetName("MC_Gamma_BP");

                Detected_Delta_BP->SetTitle("Reco KK2F");
		Detect_sans_back->SetTitle("Reco KK2F");
		Corrected_Delta_BP->SetTitle("Reco KK2F");
		All_Delta_BP->SetTitle("Reco KK2F");
		Cross_Section_BP->SetTitle("Reco KK2F");
		Gamma_BP->SetTitle("Reco KK2F");
	}

        if (option == 2) {
                Detected_Delta_BP->SetName("TFilter_Detected_Delta_BP");
		Detect_sans_back->SetName("TFilter_Detect_sans_back");
                Corrected_Delta_BP->SetName("TFilter_Corrected_Delta_BP");
                All_Delta_BP->SetName("TFilter_All_Delta_BP");
                Cross_Section_BP->SetName("TFilter_Cross_Section_BP");
                Gamma_BP->SetName("TFilter_Gamma_BP");
		
		Detected_Delta_BP->SetTitle("MC Truth");
		Detect_sans_back->SetTitle("MC Truth");
		Corrected_Delta_BP->SetTitle("MC Truth");
		All_Delta_BP->SetTitle("MC Truth");
		Cross_Section_BP->SetTitle("MC Truth");
		Gamma_BP->SetTitle("MC Truth");
	}

        if (option == 1) {
	        Detected_Delta_BP->SetTitle("Run3");
                Detect_sans_back->SetTitle("Run3");
		Corrected_Delta_BP->SetTitle("Run3");
                All_Delta_BP->SetTitle("Run3");
                Cross_Section_BP->SetTitle("Run3");
                Gamma_BP->SetTitle("Run3");

        }


	TFile *results = new TFile("results.root", "update");
		Detected_Delta_BP->Write();
		Detect_sans_back->Write();
		Corrected_Delta_BP->Write();
		All_Delta_BP->Write();
		Cross_Section_BP->Write();
		Gamma_BP->Write();
	results->Close();

                Detected_Delta_BP->Delete();
                Corrected_Delta_BP->Delete();
                All_Delta_BP->Delete();
                Cross_Section_BP->Delete();
                Gamma_BP->Delete();


}




/* Calculates the cross section from the BaBar Monte Carlo Simulation
 * uses the Monte Carlo to calculate efficiency and energy smearing
 * which is then used from the main function. it also passes the simulated
 * data to Analyze_Segment to return a value for the cross section
 */
void Monte_Carlo_Cross_Section(float mu_len[], int N_mu[], float BP_len[]) {

	TChain *ch = new TChain("ntp667");
        for (int i = 1; i <= num_files; i++) {
		ch->AddFile(TString::Format("/data/HD4/babar/mupair/KK2F/mupair_kk2f_split-BTM_TFilter-%d.root", i));
	}

	//ch->AddFile("/data/HD4/babar/mupair/Acc/mupair_kk2f_2tracks-1.root");
        //for (int i = 1; i <= 102; i++) {
	//        if (i != 17) ch->AddFile(TString::Format("Acc/mupair_kk2f_split-%d.root", i));
        //}
        
        int nevents = ch->GetEntries();
        Int_t ne, nmu, nDelta, mcLen;
        Int_t eLund[100], muTrkIdx[100], eTrkIdx[100], eMCIdx[100], mcCause[100], mothIdx[100];
        Int_t Deltad1Idx[100], Deltad2Idx[100], Deltad1Lund[100], mcLund[100];
        Int_t TRKnDchXY[100], TRKnDchZ[100], TRKnSvtXY[100], TRKnSvtZ[100], muMCIdx[100];
        Float_t pTotalScalar, mucosth[100], mccosth[100];
        Float_t eenergy[100], muenergy[100], mcenergy[100];
        Float_t DeltaDirDot[100], DeltaDOCA[100];
        Float_t DeltaePOCA0[100], DeltaePOCA1[100], DeltaePOCA2[100];
	Float_t eNeutralMag, TRKEMCecal[100];
	Float_t mcGVtxx[100], mcGVtxy[100], mcGVtxz[100];
        Float_t muBP0[100], muBP1[100], muBP2[100], muBPndot[100];
  /*      Float_t muDCH0[100], muDCH1[100], muDCH2[100], muDCHndot[100];
        Float_t muST0[100], muST1[100], muST2[100], muSTndot[100];
        Float_t muSVT10[100], muSVT11[100], muSVT12[100], muSVTndot1[100];
        Float_t muSVT20[100], muSVT21[100], muSVT22[100], muSVTndot2[100];
        Float_t muSVT30[100], muSVT31[100], muSVT32[100], muSVTndot3[100];
        Float_t muSVT40[100], muSVT41[100], muSVT42[100], muSVTndot4[100];
        Float_t muSVT50[100], muSVT51[100], muSVT52[100], muSVTndot5[100];
	Float_t ep3[100];*/
        ch->SetBranchStatus("*", 0);
	ch->SetBranchAddress("ne", &ne);
        ch->SetBranchAddress("nmu", &nmu);
        ch->SetBranchAddress("nDelta", &nDelta);
        ch->SetBranchAddress("eLund", eLund);
        ch->SetBranchAddress("muTrkIdx", muTrkIdx);
        ch->SetBranchAddress("eTrkIdx", eTrkIdx);
        ch->SetBranchAddress("Deltad1Idx", Deltad1Idx);
        ch->SetBranchAddress("Deltad2Idx", Deltad2Idx);
        ch->SetBranchAddress("Deltad1Lund", Deltad1Lund);
        ch->SetBranchAddress("TRKnDchXY", TRKnDchXY);
        ch->SetBranchAddress("TRKnDchZ", TRKnDchZ);
        ch->SetBranchAddress("TRKnSvtXY", TRKnSvtXY);
        ch->SetBranchAddress("TRKnSvtZ", TRKnSvtZ);
        ch->SetBranchAddress("eenergy", eenergy);
        ch->SetBranchAddress("muenergy", muenergy);
        ch->SetBranchAddress("DeltaDirDot", DeltaDirDot);
        ch->SetBranchAddress("DeltaDOCA", DeltaDOCA);
        ch->SetBranchAddress("pTotalScalar", &pTotalScalar);
        ch->SetBranchAddress("DeltaePOCA0", DeltaePOCA0);
        ch->SetBranchAddress("DeltaePOCA1", DeltaePOCA1);
        ch->SetBranchAddress("DeltaePOCA2", DeltaePOCA2);
	ch->SetBranchAddress("eNeutralMag", &eNeutralMag);
	ch->SetBranchAddress("TRKEMCecal", TRKEMCecal);
        ch->SetBranchAddress("mucosth", mucosth);
        ch->SetBranchAddress("muBP0", muBP0);
        ch->SetBranchAddress("muBP1", muBP1);
        ch->SetBranchAddress("muBP2", muBP2);
        ch->SetBranchAddress("muBPndot", muBPndot);
/*      ch->SetBranchAddress("muDCH0", muDCH0);
        ch->SetBranchAddress("muDCH1", muDCH1);
        ch->SetBranchAddress("muDCH2", muDCH2);
        ch->SetBranchAddress("muDCHndot", muDCHndot);
        ch->SetBranchAddress("muST0", muST0);
        ch->SetBranchAddress("muST1", muST1);
        ch->SetBranchAddress("muST2", muST2);
        ch->SetBranchAddress("muSTndot", muSTndot);
        ch->SetBranchAddress("muSVT10", muSVT10);
        ch->SetBranchAddress("muSVT11", muSVT11);
        ch->SetBranchAddress("muSVT12", muSVT12);
        ch->SetBranchAddress("muSVTndot1", muSVTndot1);
        ch->SetBranchAddress("muSVT20", muSVT20);
        ch->SetBranchAddress("muSVT21", muSVT21);
        ch->SetBranchAddress("muSVT22", muSVT22);
        ch->SetBranchAddress("muSVTndot2", muSVTndot2);
        ch->SetBranchAddress("muSVT30", muSVT30);
        ch->SetBranchAddress("muSVT31", muSVT31);
        ch->SetBranchAddress("muSVT32", muSVT32);
        ch->SetBranchAddress("muSVTndot3", muSVTndot3);
        ch->SetBranchAddress("muSVT40", muSVT40);
        ch->SetBranchAddress("muSVT41", muSVT41);
        ch->SetBranchAddress("muSVT42", muSVT42);
        ch->SetBranchAddress("muSVTndot4", muSVTndot4);
        ch->SetBranchAddress("muSVT50", muSVT50);
        ch->SetBranchAddress("muSVT51", muSVT51);
        ch->SetBranchAddress("muSVT52", muSVT52);
        ch->SetBranchAddress("muSVTndot5", muSVTndot5);*/
        ch->SetBranchAddress("mcLen", &mcLen);      
	ch->SetBranchAddress("eMCIdx", eMCIdx);
	ch->SetBranchAddress("mothIdx", mothIdx);
        ch->SetBranchAddress("mcCause", mcCause);
        ch->SetBranchAddress("mcLund", mcLund);
        ch->SetBranchAddress("mcenergy", mcenergy);
        ch->SetBranchAddress("mcGVtxx", mcGVtxx);
        ch->SetBranchAddress("mcGVtxy", mcGVtxy);
        ch->SetBranchAddress("mcGVtxz", mcGVtxz);
        ch->SetBranchAddress("muMCIdx", muMCIdx);
        ch->SetBranchAddress("mccosth", mccosth);

	// Beam Pipe variables and histograms
	Int_t n_BP_mu = 0;
	Float_t mu_len_BP = 0;
        Float_t step_size_BP = (upper_bound_BP - 0)/num_bins_BP;
        TH1F *True_All_BP = new TH1F("True_All_BP", "", num_bins_BP, 0, upper_bound_BP);
        TH1F *True_Detected_BP = new TH1F("True_Detected_BP", "", num_bins_BP, 0, upper_bound_BP);
	TH1F *Reconstructed_BP = new TH1F("Reconstructed_BP", "", num_bins_BP, 0, upper_bound_BP);
	TH1F *Sim_All_BP = new TH1F("Sim_All_BP", "", num_bins_BP, 0, upper_bound_BP);
	TH1F *Sim_Corrected_BP = new TH1F("Sim_Corrected_BP", "", num_bins_BP, 0, upper_bound_BP);
	
	True_All_BP->Sumw2();	
	True_Detected_BP->Sumw2();
	Reconstructed_BP->Sumw2();
	Sim_Corrected_BP->Sumw2();
	Sim_All_BP->Sumw2();

	// First Silicon Vertex Tracker variables and histograms
/*        Int_t n_SVT1_mu = 0;
        Float_t mu_len_SVT1 = 0;
        Float_t step_size_SVT1 = (upper_bound_SVT1 - lower_bound_SVT1)/num_bins_SVT1;
        TH1F *True_All_SVT1 = new TH1F("True_All_SVT1", "", num_bins_SVT1 + 1, lower_bound_SVT1 - step_size_SVT1, upper_bound_SVT1);
        TH1F *True_Detected_SVT1 = new TH1F("True_Detected_SVT1", "", num_bins_SVT1 + 1, lower_bound_SVT1 - step_size_SVT1, upper_bound_SVT1);
        TH1F *Reconstructed_SVT1 = new TH1F("Reconstructed_SVT1", "", num_bins_SVT1 + 1, lower_bound_SVT1 - step_size_SVT1, upper_bound_SVT1);
        TH1F *Sim_All_SVT1 = new TH1F("Sim_All_SVT1", "", num_bins_SVT1, lower_bound_SVT1, upper_bound_SVT1);
        TH1F *Sim_Corrected_SVT1 = new TH1F("Sim_Corrected_SVT1", "", num_bins_SVT1, lower_bound_SVT1, upper_bound_SVT1);
        
	True_Detected_SVT1->Sumw2();
        Reconstructed_SVT1->Sumw2();
        Sim_Corrected_SVT1->Sumw2();
        Sim_All_SVT1->Sumw2();

	// Second through fifth Silicon Vertex trackers
        Int_t n_SVT_mu = 0;
        Float_t mu_len_SVT = 0;
        Float_t step_size_SVT = (upper_bound_SVT - lower_bound_SVT)/num_bins_SVT;
        TH1F *True_All_SVT = new TH1F("True_All_SVT", "", num_bins_SVT + 1, lower_bound_SVT - step_size_SVT, upper_bound_SVT);
        TH1F *True_Detected_SVT = new TH1F("True_Detected_SVT", "", num_bins_SVT + 1, lower_bound_SVT - step_size_SVT, upper_bound_SVT);
        TH1F *Reconstructed_SVT = new TH1F("Reconstructed_SVT", "", num_bins_SVT + 1, lower_bound_SVT - step_size_SVT, upper_bound_SVT);
        TH1F *Sim_All_SVT = new TH1F("Sim_All_SVT", "", num_bins_SVT, lower_bound_SVT, upper_bound_SVT);
        TH1F *Sim_Corrected_SVT = new TH1F("Sim_Corrected_SVT", "", num_bins_SVT, lower_bound_SVT, upper_bound_SVT);
        
	True_Detected_SVT->Sumw2();
        Reconstructed_SVT->Sumw2();
        Sim_Corrected_SVT->Sumw2();
        Sim_All_SVT->Sumw2();

	// Carbon Fiber Support Tube variables and histograms
        Int_t n_ST_mu = 0;
        Float_t mu_len_ST = 0;
        Float_t step_size_ST = (upper_bound_ST - lower_bound_ST)/num_bins_ST;
        TH1F *True_All_ST = new TH1F("True_All_ST", "", num_bins_ST + 1, lower_bound_ST - step_size_ST, upper_bound_ST);
        TH1F *True_Detected_ST = new TH1F("True_Detected_ST", "", num_bins_ST + 1, lower_bound_ST - step_size_ST, upper_bound_ST);
        TH1F *Reconstructed_ST = new TH1F("Reconstructed_ST", "", num_bins_ST + 1, lower_bound_ST - step_size_ST, upper_bound_ST);
        TH1F *Sim_All_ST = new TH1F("Sim_All_ST", "", num_bins_ST, lower_bound_ST, upper_bound_ST);
        TH1F *Sim_Corrected_ST = new TH1F("Sim_Corrected_ST", "", num_bins_ST, lower_bound_ST, upper_bound_ST);
        
	True_Detected_ST->Sumw2();
        Reconstructed_ST->Sumw2();
        Sim_Corrected_ST->Sumw2();
        Sim_All_ST->Sumw2();

	// Drift Chamber histograms
        Int_t n_DCH_mu = 0;
        Float_t mu_len_DCH = 0;
        Float_t step_size_DCH = (upper_bound_DCH - lower_bound_DCH)/num_bins_DCH;
        TH1F *True_All_DCH = new TH1F("True_All_DCH", "", num_bins_DCH + 1, lower_bound_DCH - step_size_DCH, upper_bound_DCH);
        TH1F *True_Detected_DCH = new TH1F("True_Detected_DCH", "", num_bins_DCH + 1, lower_bound_DCH - step_size_DCH, upper_bound_DCH);
        TH1F *Reconstructed_DCH = new TH1F("Reconstructed_DCH", "", num_bins_DCH + 1, lower_bound_DCH - step_size_DCH, upper_bound_DCH);
        TH1F *Sim_All_DCH = new TH1F("Sim_All_DCH", "", num_bins_DCH, lower_bound_DCH, upper_bound_DCH);
        TH1F *Sim_Corrected_DCH = new TH1F("Sim_Corrected_DCH", "", num_bins_DCH, lower_bound_DCH, upper_bound_DCH);
        
	True_Detected_DCH->Sumw2();
        Reconstructed_DCH->Sumw2();
        Sim_Corrected_DCH->Sumw2();
        Sim_All_DCH->Sumw2();
*/
	TH1F *Gamma_ele = new TH1F("Gamma_ele", "", 50, 0, 300);
	TH1F *Gamma_pos = new TH1F("Gamma_pos", "", 50, 0, 300);

        TH1F *Efficiency_BP = new TH1F("Efficiency_BP", "", 60, 0, upper_bound_BP);
        TH2F *Energy_Smearing_BP = new TH2F("Energy_Smearing_BP", "",  60, 0, upper_bound_BP, 60, 0, upper_bound_BP);
	TH2F *Adet_BP = new TH2F("Adet_BP", "", 60, 0, upper_bound_BP, 60, 0, upper_bound_BP);
        TH1F *Sim_Cross_Section_BP = new TH1F("Sim_Cross_Section_BP", "", 60, 0, upper_bound_BP);

        /*TH1F *Efficiency_SVT1 = new TH1F("Efficiency_SVT1", "", num_bins_SVT1, lower_bound_SVT1, upper_bound_SVT1);
        TH2F *Energy_Smearing_SVT1 = new TH2F("Energy_Smearing_SVT1", "",  num_bins_SVT1 + 1, lower_bound_SVT1 - step_size_SVT1, upper_bound_SVT1, num_bins_SVT1 + 1, lower_bound_SVT1 - step_size_SVT1, upper_bound_SVT1);
        TH1F *Sim_Cross_Section_SVT1 = new TH1F("Sim_Cross_Section_SVT1", "", num_bins_SVT1, lower_bound_SVT1, upper_bound_SVT1);

        TH1F *Efficiency_SVT = new TH1F("Efficiency_SVT", "", num_bins_SVT, lower_bound_SVT, upper_bound_SVT);
        TH2F *Energy_Smearing_SVT = new TH2F("Energy_Smearing_SVT", "",  num_bins_SVT + 1, lower_bound_SVT - step_size_SVT, upper_bound_SVT, num_bins_SVT + 1, lower_bound_SVT - step_size_SVT, upper_bound_SVT);
        TH1F *Sim_Cross_Section_SVT = new TH1F("Sim_Cross_Section_SVT", "", num_bins_SVT, lower_bound_SVT, upper_bound_SVT);

        TH1F *Efficiency_ST = new TH1F("Efficiency_ST", "", num_bins_ST, lower_bound_ST, upper_bound_ST);
        TH2F *Energy_Smearing_ST = new TH2F("Energy_Smearing_ST", "",  num_bins_ST + 1, lower_bound_ST - step_size_ST, upper_bound_ST, num_bins_ST + 1, lower_bound_ST - step_size_ST, upper_bound_ST);
        TH1F *Sim_Cross_Section_ST = new TH1F("Sim_Cross_Section_ST", "", num_bins_ST, lower_bound_ST, upper_bound_ST);

        TH1F *Efficiency_DCH = new TH1F("Efficiency_DCH", "", num_bins_DCH, lower_bound_DCH, upper_bound_DCH);
        TH2F *Energy_Smearing_DCH = new TH2F("Energy_Smearing_DCH", "",  num_bins_DCH + 1, lower_bound_DCH - step_size_DCH, upper_bound_DCH, num_bins_DCH + 1, lower_bound_DCH - step_size_DCH, upper_bound_DCH);
        TH1F *Sim_Cross_Section_DCH = new TH1F("Sim_Cross_Section_DCH", "", num_bins_DCH, lower_bound_DCH, upper_bound_DCH);
*/
	TH2F *Mu_En_Dist = new TH2F("Mu_En_Dist", "", 60, 0, 300, 5, 3000, 8000);
	TH2F *Mu_En_Cross_Temp = new TH2F("Mu_En_Cross_Temp", "", 60, 0, 300, 5, 3000, 8000);
        TH2F *Mu_En_Cross = new TH2F("Mu_En_Cross", "", 60, 0, 300, 5, 3000, 8000);

	Mu_En_Dist->Sumw2();
	Mu_En_Cross->Sumw2();

	cout << endl;
	cout << "Reading TFilter Monte Carlo Data" << endl;

	//TH2F *POCA_3 = new TH2F("POCA_3", "", 400, -25, 25, 400, -25, 25);
	int num_slices = 40;
	float norm = TMath::Sqrt(2);
	float slice_size = norm*(upper_bound_BP - lower_bound_BP - step_size_BP)/num_slices;
	vector<TH1F*> hists;
	for (int i = 0; i < num_slices; i++) {
		TH1F *hist = new TH1F(TString::Format("hist%d", i), "", 200, -50, 50);
		hists.push_back(hist);
	}


	// Actually fill the histograms.
        for (int i = 0; i < nevents; i++) {
                ch->GetEntry(i);
                if (i % (nevents/10) == 0) cout << i << "/" << nevents << endl;
                
                if (nmu == 2 && TMath::Abs(pTotalScalar - 11.75) <= 0.75 && eNeutralMag < 1.0) {

                        // Look at individual muon tracks to see where they hit
			for (int j = 0; j < nmu; j++) {
                        	if (muenergy[j] > 3.5 && muenergy[j] < 8 && TRKEMCecal[muTrkIdx[j]] < 0.5) {
                                	if (true) { //TRKnSvtXY[muTrkIdx[j]] >= 4 && TRKnSvtZ[muTrkIdx[j]] >= 4) {
		                        if (TRKnDchXY[muTrkIdx[j]] >= 15 && TRKnDchZ[muTrkIdx[j]] >= 10) {
                                        	// Select muon track through Beam Pipe
                                         	float mu_rad = TMath::Sqrt(TMath::Power(muBP0[j],2) + TMath::Power(muBP1[j], 2));
                                                if (mu_rad > BP_low_rad && mu_rad < BP_high_rad) {
                                                	if (muBP2[j] > BP_z_low && muBP2[j] < BP_z_high && muBPndot[j] > 0.5) {

                                                                // Was a delta ray produced?
                                                                for (int k = 0; k < mcLen; k++) {
                                                                	if (mothIdx[k] == muMCIdx[j] && mcCause[k] == 202) {
                                                                        	float mc_Rad = TMath::Sqrt(TMath::Power(mcGVtxx[k], 2) + TMath::Power(mcGVtxy[k], 2));
                                                                                if (mc_Rad > BP_low_rad && mc_Rad < BP_high_rad) {
                                                                	                if (mcGVtxz[k] > BP_z_low && mcGVtxz[k] < BP_z_high) {
												True_All_BP->Fill(1000*mcenergy[k]);
												//Mu_En_Cross_Temp->Fill(1000*mcenergy[k], 1000*muenergy[j]);
												Mu_En_Cross_Temp->Fill(1000*mcenergy[k], 1000*mcenergy[muMCIdx[j]]);

												// was delta ray seen?
 												for (int l = 0; l < nDelta; l++) {
											        	if (Deltad2Idx[l] == j) {
											       	        	if (eMCIdx[Deltad1Idx[l]] == k) {
												                        if (Deltad1Lund[l] == 11 && DeltaDirDot[l] > 0.98 && DeltaDOCA[l] < 0.1) {
												                                Float_t Delta_Rad = TMath::Sqrt(TMath::Power(DeltaePOCA0[l], 2) + TMath::Power(DeltaePOCA1[l], 2));
												                                if (Delta_Rad > BP_low_rad && Delta_Rad < BP_high_rad) {
												                                        if (DeltaePOCA2[l] > BP_z_low && DeltaePOCA2[l] < BP_z_high) {
												                                        	True_Detected_BP->Fill(1000*mcenergy[k]);
																		Reconstructed_BP->Fill(1000*eenergy[Deltad1Idx[l]]);
																		//Energy_Smearing_BP->Fill(1000*mcenergy[k], 1000*eenergy[Deltad1Idx[l]]);
																		Mu_En_Dist->Fill(1000*mcenergy[k], 1000*muenergy[j]);
																		//POCA_3->Fill(DeltaePOCA0[l], DeltaePOCA1[l]);
																		//ToyMC stuff
																		//Float_t True_prime = 1000*mcenergy[k] + 1000*eenergy[Deltad1Idx[l]];
																		//Float_t Recon_prime = -1000*mcenergy[k] + 1.1*1000*eenergy[Deltad1Idx[l]] + 5;
																		//if (TMath::Abs(Recon_prime) <= 50) {
																			Energy_Smearing_BP->Fill(1000*mcenergy[k], 1000*eenergy[Deltad1Idx[l]]);

						// TSVDUnfold requires axes flipped
						Adet_BP->Fill(1000*eenergy[Deltad1Idx[l]], 1000*mcenergy[k]);
																		//}
																		//for (int slice = 0; slice < num_slices; slice++) {
																		//	if (True_prime >=(norm*slice*slice_size) && True_prime < (norm*(slice + 1)*slice_size)) {
																		//		hists.at(slice)->Fill(Recon_prime);
																		//	}
																		//}
                                         												}
                                 												}
                         												}
                 												}
         												}
 												}
                                                                                        }
										}
                                                                                }
                                                                        
										if (mcCause[k] == 201) {
											if (mcLund[k] == 11) {
												Gamma_ele->Fill(1000*mcenergy[k]);
											}
											if (mcLund[k] == -11) {
												Gamma_pos->Fill(1000*mcenergy[k]);
											}
									}
									}
                                                                }
                                                        }
						}
/*
						//First Silicon Vertex Tracker
						if (TMath::Sqrt(TMath::Power(muSVT10[j], 2) + TMath::Power(muSVT11[j], 2)) > 1) {
							if (muSVT12[j] > SVT1_z_low && muSVT12[j] < SVT1_z_high) {

								for (int k = 0; k < mcLen; k++) {
									if (mothIdx[k] == muMCIdx[j]) {
										Float_t mc_Rad = TMath::Sqrt(TMath::Power(mcGVtxx[k], 2) + TMath::Power(mcGVtxy[k], 2));
										if (mc_Rad > SVT1_low_rad && mc_Rad < SVT1_high_rad) {
											if (mcGVtxz[k] > SVT1_z_low && mcGVtxz[k] < SVT1_z_high) {
												True_All_SVT1->Fill(1000*mcenergy[k]);

												for (int l = 0; l < nDelta; l++) {
													if (Deltad2Idx[l] == j) {
														if (eMCIdx[Deltad1Idx[l]] == k) {
															if (Deltad1Lund[l] == 11 && DeltaDirDot[l] > 0.98) {
                                                                                                                                Float_t Delta_Rad = TMath::Sqrt(TMath::Power(DeltaePOCA0[l], 2) + TMath::Power(DeltaePOCA1[l], 2));
                                                                                                                                if (Delta_Rad > SVT1_low_rad && Delta_Rad < SVT1_high_rad) {
                                                                                                                                        if (DeltaePOCA2[l] > SVT1_z_low && DeltaePOCA2[l] < SVT1_z_high) {

                                                                                                                                                if((TRKnSvtXY[eTrkIdx[Deltad1Idx[l]]] > 3) && (TRKnSvtZ[eTrkIdx[Deltad1Idx[l]]] > 3)){
                                                                                                                                                        Float_t elec_trans_new = TMath::Sqrt(TMath::Power(DeltaePOCA0[l] + 0.01, 2) + TMath::Power(DeltaePOCA1[l] - 0.013, 2));
                                                                                                                                                        Float_t phi = TMath::ATan2(DeltaePOCA1[j] - 0.013, DeltaePOCA0[l] + 0.01) + 2*TMath::Pi() + 0.42705;
                                                                                                                                                        Float_t new_phi = fmodf(TMath::Abs(phi), TMath::Pi()/3.0);
                                                                                                                                                        Float_t new_x = elec_trans_new*TMath::Cos(new_phi + TMath::Pi()/3.0);
                                                                                                                                                        Float_t new_y = elec_trans_new*TMath::Sin(new_phi + TMath::Pi()/3.0);

                                                                                                                                                        if ((new_x >= SVT1_low_x) && (new_x <= SVT1_high_x) && (new_y >= SVT1_low_y) && (new_y <= SVT1_high_y)) {
                                                                                                                                                                True_Detected_SVT1->Fill(1000*mcenergy[k]);
                                                                                                                                               	 		Reconstructed_SVT1->Fill(1000*eenergy[l]);
                                                                                                                                                                Energy_Smearing_SVT1->Fill(1000*mcenergy[k], 1000*eenergy[l]);
                                                                                                                                                POCA_3->Fill(DeltaePOCA0[l], DeltaePOCA1[l]);
                                                                                                                                                        }
                                                                                                                                                }
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}

						//Second through fifth SVT's
                                                if (TMath::Sqrt(TMath::Power(muSVT20[j], 2) + TMath::Power(muSVT21[j], 2)) > 1) {
                                                        if (muSVT22[j] > SVT2_z_low && muSVT22[j] < SVT2_z_high) {

                                                                for (int k = 0; k < mcLen; k++) {
                                                                        if (mothIdx[k] == muMCIdx[j]) {
                                                                                Float_t mc_Rad = TMath::Sqrt(TMath::Power(mcGVtxx[k], 2) + TMath::Power(mcGVtxy[k], 2));
                                                                                if (mc_Rad > SVT2_low_rad && mc_Rad < SVT2_high_rad) {
                                                                                        if (mcGVtxz[k] > SVT2_z_low && mcGVtxz[k] < SVT2_z_high) {
                                                                                                True_All_SVT->Fill(1000*mcenergy[k]);

                                                                                                for (int l = 0; l < nDelta; l++) {
                                                                                                        if (Deltad2Idx[l] == j) {
                                                                                                                if (eMCIdx[Deltad1Idx[l]] == k) {
                                                                                                                        if (Deltad1Lund[l] == 11 && DeltaDirDot[l] > 0.98) {
                                                                                                                                Float_t Delta_Rad = TMath::Sqrt(TMath::Power(DeltaePOCA0[l], 2) + TMath::Power(DeltaePOCA1[l], 2));
                                                                                                                                if (Delta_Rad > SVT2_low_rad && Delta_Rad < SVT2_high_rad) {
                                                                                                                                        if (DeltaePOCA2[l] > SVT2_z_low && DeltaePOCA2[l] < SVT2_z_high) {

                                                                                                                                                if((TRKnSvtXY[eTrkIdx[Deltad1Idx[l]]] > 3) && (TRKnSvtZ[eTrkIdx[Deltad1Idx[l]]] > 3)){
                                                                                                                                                        Float_t elec_trans_new = TMath::Sqrt(TMath::Power(DeltaePOCA0[l] + 0.01, 2) + TMath::Power(DeltaePOCA1[l] - 0.013, 2));
                                                                                                                                                        Float_t phi = TMath::ATan2(DeltaePOCA1[j] - 0.013, DeltaePOCA0[l] + 0.01) + 2*TMath::Pi() + 0.42705;
                                                                                                                                                        Float_t new_phi = fmodf(TMath::Abs(phi), TMath::Pi()/3.0);
                                                                                                                                                        Float_t new_x = elec_trans_new*TMath::Cos(new_phi + TMath::Pi()/3.0);
                                                                                                                                                        Float_t new_y = elec_trans_new*TMath::Sin(new_phi + TMath::Pi()/3.0);

                                                                                                                                                        if ((new_x >= SVT2_low_x) && (new_x <= SVT2_high_x) && (new_y >= SVT2_low_y) && (new_y <= SVT2_high_y)) {
                                                                                                                                                                True_Detected_SVT->Fill(1000*mcenergy[k]);
                                                                                                                                                		Reconstructed_SVT->Fill(1000*eenergy[l]);
                                                                                                                                                                Energy_Smearing_SVT->Fill(1000*mcenergy[k], 1000*eenergy[l]);
                                                                                                                                                                                                                                                                                                        POCA_3->Fill(DeltaePOCA0[l], DeltaePOCA1[l]);
																			}
                                                                                                                                                }
                                                                                                                                        }
                                                                                                                                }
                                                                                                                        }
                                                                                                                }
                                                                                                        }
                                                                                                }
                                                                                        }
                                                                                }
                                                                        }
                                                                }
                                                        }
                                                }
                                                if (TMath::Sqrt(TMath::Power(muSVT30[j], 2) + TMath::Power(muSVT31[j], 2)) > 1) {
                                                        if (muSVT32[j] > SVT3_z_low && muSVT32[j] < SVT3_z_high) {

                                                                for (int k = 0; k < mcLen; k++) {
                                                                        if (mothIdx[k] == muMCIdx[j]) {
                                                                                Float_t mc_Rad = TMath::Sqrt(TMath::Power(mcGVtxx[k], 2) + TMath::Power(mcGVtxy[k], 2));
                                                                                if (mc_Rad > SVT3_low_rad && mc_Rad < SVT3_high_rad) {
                                                                                        if (mcGVtxz[k] > SVT3_z_low && mcGVtxz[k] < SVT3_z_high) {
                                                                                                True_All_SVT->Fill(1000*mcenergy[k]);

                                                                                                for (int l = 0; l < nDelta; l++) {
                                                                                                        if (Deltad2Idx[l] == j) {
                                                                                                                if (eMCIdx[Deltad1Idx[l]] == k) {
                                                                                                                        if (Deltad1Lund[l] == 11 && DeltaDirDot[l] > 0.98) {
                                                                                                                                Float_t Delta_Rad = TMath::Sqrt(TMath::Power(DeltaePOCA0[l], 2) + TMath::Power(DeltaePOCA1[l], 2));
                                                                                                                                if (Delta_Rad > SVT3_low_rad && Delta_Rad < SVT3_high_rad) {
                                                                                                                                        if (DeltaePOCA2[l] > SVT3_z_low && DeltaePOCA2[l] < SVT3_z_high) {

                                                                                                                                                if((TRKnSvtXY[eTrkIdx[Deltad1Idx[l]]] > 3) && (TRKnSvtZ[eTrkIdx[Deltad1Idx[l]]] > 3)){
                                                                                                                                                        Float_t elec_trans_new = TMath::Sqrt(TMath::Power(DeltaePOCA0[l] + 0.01, 2) + TMath::Power(DeltaePOCA1[l] - 0.013, 2));
                                                                                                                                                        Float_t phi = TMath::ATan2(DeltaePOCA1[j] - 0.013, DeltaePOCA0[l] + 0.01) + 2*TMath::Pi() + 0.42705;
                                                                                                                                                        Float_t new_phi = fmodf(TMath::Abs(phi), TMath::Pi()/3.0);
                                                                                                                                                        Float_t new_x = elec_trans_new*TMath::Cos(new_phi + TMath::Pi()/3.0);
                                                                                                                                                        Float_t new_y = elec_trans_new*TMath::Sin(new_phi + TMath::Pi()/3.0);

                                                                                                                                                        if ((new_x >= SVT3_low_x) && (new_x <= SVT3_high_x) && (new_y >= SVT3_low_y) && (new_y <= SVT3_high_y)) {
                                                                                                                                                                True_Detected_SVT->Fill(1000*mcenergy[k]);
                                                                                                                                                		Reconstructed_SVT->Fill(1000*eenergy[l]);
                                                                                                                                                                Energy_Smearing_SVT->Fill(1000*mcenergy[k], 1000*eenergy[l]);
                                                                                                                                                                                                POCA_3->Fill(DeltaePOCA0[l], DeltaePOCA1[l]); 
						                                                                                                       }
                                                                                                                                                }
                                                                                                                                        }
                                                                                                                                }
                                                                                                                        }
                                                                                                                }
                                                                                                        }
                                                                                                }
                                                                                        }
                                                                                }
                                                                        }
                                                                }
                                                        }
                                                }

						//Carbon Fiber Support Tube
                                                if (TMath::Abs(TMath::Sqrt(TMath::Power(muST0[j], 2) + TMath::Power(muST1[j], 2)) - 21) < 1) {
                                                        if (muST2[j] > ST_z_low && muST2[j] < ST_z_high) {

                                                                for (int k = 0; k < mcLen; k++) {
                                                                        if (mothIdx[k] == muMCIdx[j]) {
                                                                                Float_t mc_Rad = TMath::Sqrt(TMath::Power(mcGVtxx[k], 2) + TMath::Power(mcGVtxy[k], 2));
                                                                                if (mc_Rad > ST_low_rad && mc_Rad < ST_high_rad) {
                                                                                        if (mcGVtxz[k] > ST_z_low && mcGVtxz[k] < ST_z_high) {
                                                                                                True_All_ST->Fill(1000*mcenergy[k]);
 
                                                                                                for (int l = 0; l < nDelta; l++) {
                                                                                                        if (Deltad2Idx[l] == j) {
                                                                                                                if (eMCIdx[Deltad1Idx[l]] == k) {
                                                                                                                        if (Deltad1Lund[l] == 11 && DeltaDirDot[l] > 0.98) {
                                                                                                                                Float_t Delta_Rad = TMath::Sqrt(TMath::Power(DeltaePOCA0[l], 2) + TMath::Power(DeltaePOCA1[l], 2));
                                                                                                                                if (Delta_Rad > ST_low_rad && Delta_Rad < ST_high_rad) {
                                                                                                                                        if (DeltaePOCA2[l] > ST_z_low && DeltaePOCA2[l] < ST_z_high) {
                                                                                                                                                True_Detected_ST->Fill(1000*mcenergy[k]);
                                                                                                                                                Reconstructed_ST->Fill(1000*eenergy[l]);
																		Energy_Smearing_ST->Fill(1000*mcenergy[k], 1000*eenergy[l]);
                                                                                                                                                POCA_3->Fill(DeltaePOCA0[l], DeltaePOCA1[l]);
                                                                                                                                        }
                                                                                                                                }
                                                                                                                        }
                                                                                                                }
                                                                                                        }
                                                                                                }

                                                                                        }
                                                                                }
                                                                        }
                                                                }
                                                        }
                                                }

                                                //muon path through Drift Chamber
                                                if (TMath::Sqrt(TMath::Power(muDCH0[j], 2) + TMath::Power(muDCH1[j], 2)) > 1) {
                                              	 	if (muDCH2[j] > DCH_z_low && muDCH2[j] < DCH_z_high) {

                                                                for (int k = 0; k < mcLen; k++) {
                                                                	if (mothIdx[k] == muMCIdx[j]) {
                                                                        	Float_t mc_Rad = TMath::Sqrt(TMath::Power(mcGVtxx[k], 2) + TMath::Power(mcGVtxy[k], 2));
                                                                                if (mc_Rad > DCH_low_rad && mc_Rad < DCH_high_rad) {
                                                                                	if (mcGVtxz[k] > DCH_z_low && mcGVtxz[k] < DCH_z_high) {
												True_All_DCH->Fill(1000*mcenergy[k]);

												for (int l = 0; l < nDelta; l++) {
													if (Deltad2Idx[l] == j) {
														if (eMCIdx[Deltad1Idx[l]] == k) {
															if (Deltad1Lund[l] == 11 && DeltaDirDot[l] > 0.98) {
																Float_t Delta_Rad = TMath::Sqrt(TMath::Power(DeltaePOCA0[l], 2) + TMath::Power(DeltaePOCA1[l], 2));
																if (Delta_Rad > DCH_low_rad && Delta_Rad < DCH_high_rad) {
																	if (DeltaePOCA2[l] > DCH_z_low && DeltaePOCA2[l] < DCH_z_high) {
                                                                                        							True_Detected_DCH->Fill(1000*mcenergy[k]);
                                                                                                                                                Reconstructed_DCH->Fill(1000*eenergy[l]);
																		Energy_Smearing_DCH->Fill(1000*mcenergy[k], 1000*eenergy[l]);
                                                                                                                                                POCA_3->Fill(DeltaePOCA0[l], DeltaePOCA1[l]);
																	}
																}
															}
														}
													}
												}
                                                                                        }
                                                                                }
                                                                        }
                                                                }
                                                        }
						}*/
					}
				}
			}
		}
	}

	

	//TH1F *pull = new TH1F("pull", "", 100, -10, 10);
	//TMatrix mat(52, 52);
	//TH2F *Smearing_Correction_BP = new TH2F("Smearing_Correction_BP", "Inverted Correction of Energy Smearing", 52, 40, 300, 52, 40, 300);
	//Invert_Matrix(Energy_Smearing_BP, Smearing_Correction_BP, mat, 52, true);
	Sim_Efficiency(Efficiency_BP, True_All_BP, True_Detected_BP, num_bins_BP);

	std::cout << "scale : " << 1./(mu_len[0]*density_BP*step_size_BP) << std::endl;

	//True_All_BP->Scale(1./(mu_len[0]*density_BP*step_size_BP));
	//True_Detected_BP->Scale(1./(mu_len[0]*density_BP*step_size_BP));
        
	//ToyMC(Efficiency_BP, Energy_Smearing_BP, hists, True_All_BP, pull, num_slices, slice_size, step_size_BP);
	//Analyze_Segment(Reconstructed_BP, Sim_Corrected_BP, Sim_All_BP, Sim_Cross_Section_BP, Efficiency_BP, Energy_Smearing_BP, mu_len[0], density_BP, lower_bound_BP, upper_bound_BP, N_mu[0], num_bins_BP);
/*        
        Sim_Efficiency(Efficiency_SVT1, True_All_SVT1, True_Detected_SVT1, num_bins_SVT1 + 1);
        Analyze_Segment(Reconstructed_SVT1, Sim_Corrected_SVT1, Sim_All_SVT1, Sim_Cross_Section_SVT1, Efficiency_SVT1, Energy_Smearing_SVT1, mu_len_SVT1, density_SVT, lower_bound_SVT1, upper_bound_SVT1, n_SVT1_mu, num_bins_SVT1);

        Sim_Efficiency(Efficiency_SVT, True_All_SVT, True_Detected_SVT, num_bins_SVT + 1);
        Analyze_Segment(Reconstructed_SVT, Sim_Corrected_SVT, Sim_All_SVT, Sim_Cross_Section_SVT, Efficiency_SVT, Energy_Smearing_SVT, mu_len_SVT, density_SVT, lower_bound_SVT, upper_bound_SVT, n_SVT_mu, num_bins_SVT);

        Sim_Efficiency(Efficiency_ST, True_All_ST, True_Detected_ST, num_bins_ST + 1);
        Analyze_Segment(Reconstructed_ST, Sim_Corrected_ST, Sim_All_ST, Sim_Cross_Section_ST, Efficiency_ST, Energy_Smearing_ST, mu_len_ST, density_ST, lower_bound_ST, upper_bound_ST, n_ST_mu, num_bins_ST);

	Sim_Efficiency(Efficiency_DCH, True_All_DCH, True_Detected_DCH, num_bins_DCH + 1);
	Analyze_Segment(Reconstructed_DCH, Sim_Corrected_DCH, Sim_All_DCH, Sim_Cross_Section_DCH, Efficiency_DCH, Energy_Smearing_DCH, mu_len_DCH, density_DCH, lower_bound_DCH, upper_bound_DCH, n_DCH_mu, num_bins_DCH);
*/        


	// This calculates the differential cross section for delta ray energy and muon energy
	// each bin of muon energies needs to be normalized by the respective number of muons
	// with that energy.
	TH1F *Muon_Cross[5];
	for (int col = 1; col <= 5; col++) {
		Muon_Cross[col - 1] = new TH1F(TString::Format("Muon_Cross_%d", col), TString::Format("Muon Cross Section %d GeV; Energy (MeV); Production Cross Section (cm^{2}g^{-1}MeV^{-1})", col + 2), 60, 0, 300);
		for (int row = 1; row <= 60; row++) {
			//Mu_En_Dist->SetBinContent(row, col, Mu_En_Dist->GetBinContent(row, col)/(BP_len[row - 1]*density_BP*5));
                        //Mu_En_Dist->SetBinError(row, col, Mu_En_Dist->GetBinError(row, col)/(BP_len[row - 1]*density_BP*5));

			Mu_En_Cross->SetBinContent(row, col, Mu_En_Cross_Temp->GetBinContent(row, col)/(BP_len[col - 1]*density_BP*5.0*1000));
                        Mu_En_Cross->SetBinError(row, col, Mu_En_Cross_Temp->GetBinError(row, col)/(BP_len[col - 1]*density_BP*5.0*1000));

			Muon_Cross[col - 1]->SetBinContent(row, Mu_En_Cross->GetBinContent(row, col));
			Muon_Cross[col - 1]->SetBinError(row, Mu_En_Cross->GetBinError(row, col));
			//cout << "At " << 1000*(col + 2) << "GeV   " << BP_len[col - 1] << endl;
		}
	}

	// Save the needed histograms to a file
	// if (num_files == 102) {
		TFile f("results.root", "recreate");
	        Efficiency_BP->Write();
        	Energy_Smearing_BP->Write();
                True_Detected_BP->Write();
                True_All_BP->Write();
		Reconstructed_BP->Write();
		Sim_Corrected_BP->Write();
		Sim_All_BP->Write();
		Sim_Cross_Section_BP->Write();
        	//Smearing_Correction_BP->Write();
		Mu_En_Cross->Write();
		Mu_En_Dist->Write();
		for (int i = 0; i < 5; i++) {
			Muon_Cross[i]->Write();
		}
		Adet_BP->Write();
/*

		Efficiency_SVT1->Write();
	        Energy_Smearing_SVT1->Write();
	        Efficiency_SVT->Write();
	        Energy_Smearing_SVT->Write();
	        Efficiency_ST->Write();
	        Energy_Smearing_ST->Write();
	        Efficiency_DCH->Write();
	        Energy_Smearing_DCH->Write();*/
        	f.Close();
	/*}
	else {
                TFile f("results.root", "update");
                True_Detected_BP->Write();
                True_All_BP->Write();
                Reconstructed_BP->Write();
                Sim_Corrected_BP->Write();
                Sim_All_BP->Write();
                Sim_Cross_Section_BP->Write();

                f.Close();
	}
*/

}


// are prescaled by filescale, so the final numbers will be multiplied by this number to get the 
// correct amount of muons. This code is identical for real and Monte Carlo data.
void Calc_Acceptance(TString filename, float mu_len[], int N_mu[], float BP_len[]) {

	TChain *chain = new TChain("ntp667");
	for (int i = 1; i <= num_files; i++) {
		chain->AddFile(TString::Format(filename + "%d" + ".root", i));
	}

	Int_t nevents = chain->GetEntries();
	Int_t nmu, muTrkIdx[100], TRKnSvtXY[100], TRKnSvtZ[100];
	Int_t TRKnDchXY[100], TRKnDchZ[100];
	Float_t eNeutralMag, TRKEMCecal[100];
	Float_t pTotalScalar, muenergy[100], muBPndot[100];
	Float_t muBP0[100], muBP1[100], muBP2[100];
        Float_t muDCH0[100], muDCH1[100], muDCH2[100], muDCHndot[100];
        Float_t muST0[100], muST1[100], muST2[100], muSTndot[100];
        Float_t muSVT10[100], muSVT11[100], muSVT12[100], muSVTndot1[100];
        Float_t muSVT20[100], muSVT21[100], muSVT22[100], muSVTndot2[100];
        Float_t muSVT30[100], muSVT31[100], muSVT32[100], muSVTndot3[100];
        Float_t muSVT40[100], muSVT41[100], muSVT42[100], muSVTndot4[100];
        Float_t muSVT50[100], muSVT51[100], muSVT52[100], muSVTndot5[100];

        chain->SetBranchStatus("*", 0);
	chain->SetBranchAddress("nmu", &nmu);
	chain->SetBranchAddress("muTrkIdx", muTrkIdx);
	chain->SetBranchAddress("eNeutralMag", &eNeutralMag);
	chain->SetBranchAddress("TRKEMCecal", TRKEMCecal);
        chain->SetBranchAddress("TRKnSvtXY", TRKnSvtXY);
        chain->SetBranchAddress("TRKnSvtZ", TRKnSvtZ);
        chain->SetBranchAddress("TRKnDchXY", TRKnDchXY);
        chain->SetBranchAddress("TRKnDchZ", TRKnDchZ);
        chain->SetBranchAddress("pTotalScalar", &pTotalScalar);
        chain->SetBranchAddress("muenergy", muenergy);
        chain->SetBranchAddress("muBPndot", muBPndot);
        chain->SetBranchAddress("muBP0", muBP0);
        chain->SetBranchAddress("muBP1", muBP1);
        chain->SetBranchAddress("muBP2", muBP2);
/*        chain->SetBranchAddress("muDCH0", muDCH0);
        chain->SetBranchAddress("muDCH1", muDCH1);
        chain->SetBranchAddress("muDCH2", muDCH2);
        chain->SetBranchAddress("muDCHndot", muDCHndot);
        chain->SetBranchAddress("muST0", muST0);
        chain->SetBranchAddress("muST1", muST1);
        chain->SetBranchAddress("muST2", muST2);
        chain->SetBranchAddress("muSTndot", muSTndot);
        chain->SetBranchAddress("muSVT10", muSVT10);
        chain->SetBranchAddress("muSVT11", muSVT11);
        chain->SetBranchAddress("muSVT12", muSVT12);
        chain->SetBranchAddress("muSVTndot1", muSVTndot1);
        chain->SetBranchAddress("muSVT20", muSVT20);
        chain->SetBranchAddress("muSVT21", muSVT21);
        chain->SetBranchAddress("muSVT22", muSVT22);
        chain->SetBranchAddress("muSVTndot2", muSVTndot2);
        chain->SetBranchAddress("muSVT30", muSVT30);
        chain->SetBranchAddress("muSVT31", muSVT31);
        chain->SetBranchAddress("muSVT32", muSVT32);
        chain->SetBranchAddress("muSVTndot3", muSVTndot3);
        chain->SetBranchAddress("muSVT40", muSVT40);
        chain->SetBranchAddress("muSVT41", muSVT41);
        chain->SetBranchAddress("muSVT42", muSVT42);
        chain->SetBranchAddress("muSVTndot4", muSVTndot4);
        chain->SetBranchAddress("muSVT50", muSVT50);
        chain->SetBranchAddress("muSVT51", muSVT51);
        chain->SetBranchAddress("muSVT52", muSVT52);
        chain->SetBranchAddress("muSVTndot5", muSVTndot5);
*/
	cout << endl;
	cout << "Calculating Acceptance" << endl;

	Float_t mu_len_BP = 0;
	Int_t N_mu_BP = 0;
        Float_t mu_len_SVT1 = 0;
        Int_t N_mu_SVT1 = 0;
        Float_t mu_len_SVT = 0;
        Int_t N_mu_SVT = 0;
        Float_t mu_len_ST = 0;
        Int_t N_mu_ST = 0;
        Float_t mu_len_DCH = 0;
        Int_t N_mu_DCH = 0;

	for (int i = 0; i < 5; i++) {
		BP_len[i] = 0;
	}
	int num_mu = 0;

        for (int i = 0; i < nevents; i++) {
                chain->GetEntry(i);
                if (i % (nevents/10) == 0) cout << i << "/" << nevents << endl;

                if (nmu == 2 && TMath::Abs(pTotalScalar - 11.75) < 0.75 && eNeutralMag < 1.0) {

                        // Look at individual muon tracks to see where they hit
                        for (int j = 0; j < nmu; j++) {
                                if (muenergy[j] > 3.5 && muenergy[j] < 8 && TRKEMCecal[muTrkIdx[j]] < 0.5) {
				if (TRKnDchXY[muTrkIdx[j]] >= 15 && TRKnDchZ[muTrkIdx[j]] >= 10) {
                                        if (true) { //TRKnSvtXY[muTrkIdx[j]] >= 4 && TRKnSvtZ[muTrkIdx[j]] >= 4) {						
					float mu_rad = TMath::Sqrt(TMath::Power(muBP0[j], 2) + TMath::Power(muBP1[j], 2));

                                                // Select muon track through Beam Pipe
                                                if (mu_rad > BP_low_rad && mu_rad < BP_high_rad) {
                                                        if (muBP2[j] > BP_z_low && muBP2[j] < BP_z_high && muBPndot[j] > 0.5) {
								N_mu_BP++;
								mu_len_BP = mu_len_BP + 0.284/muBPndot[j];

								for (int l = 0; l < 5; l++) {
									if (muenergy[j] > (3 + l) && muenergy[j] <= 4 + l) {
										
										BP_len[l] = BP_len[l] + 0.284/muBPndot[j];
										num_mu++;

									}
								}
							}
						}
	/*					//First Silicon Vertex Tracker
                                                if (TMath::Sqrt(TMath::Power(muSVT10[j], 2) + TMath::Power(muSVT11[j], 2)) > 1) {
                                                        if (muSVT12[j] > SVT1_z_low && muSVT12[j] < SVT1_z_high) {
                                                                N_mu_SVT1++;
                                                                mu_len_SVT1 = mu_len_SVT1 + 0.04/muSVTndot1[j];
							}
						}
						//Second through Fifth SVT
                                                if (TMath::Sqrt(TMath::Power(muSVT20[j], 2) + TMath::Power(muSVT21[j], 2)) > 1) {
                                                        if (muSVT22[j] > SVT2_z_low && muSVT22[j] < SVT2_z_high) {
                                                                N_mu_SVT++;
                                                                mu_len_SVT = mu_len_SVT + 0.04/muSVTndot2[j];
							}
						}
                                                if (TMath::Sqrt(TMath::Power(muSVT30[j], 2) + TMath::Power(muSVT31[j], 2)) > 1) {
                                                        if (muSVT32[j] > SVT3_z_low && muSVT32[j] < SVT3_z_high) {
                                                                N_mu_SVT++;
                                                                mu_len_SVT = mu_len_SVT + 0.04/muSVTndot3[j];

							}
						}
                                                if (TMath::Sqrt(TMath::Power(muSVT30[j], 2) + TMath::Power(muSVT31[j], 2)) > 1) {
                                                        if (muSVT32[j] > SVT3_z_low && muSVT32[j] < SVT3_z_high) {
                                                                N_mu_SVT++;
                                                                mu_len_SVT = mu_len_SVT + 0.04/muSVTndot3[j];
							}
						}
						//Support Tube
                                                if (TMath::Abs(TMath::Sqrt(TMath::Power(muST0[j], 2) + TMath::Power(muST1[j], 2)) - 21) < 1) {
                                                        if (muST2[j] > ST_z_low && muST2[j] < ST_z_high) {
                                                                mu_len_ST = mu_len_ST + 0.208/muSTndot[j];
                                                                N_mu_ST++;
							}
						}
						//Support Tube
                                                if (TMath::Sqrt(TMath::Power(muDCH0[j], 2) + TMath::Power(muDCH1[j], 2)) > 1) {
                                                        if (muDCH2[j] > DCH_z_low && muDCH2[j] < DCH_z_high) {
                                                                mu_len_DCH = mu_len_DCH + 0.11/muDCHndot[j];
                                                                N_mu_DCH++;
							}
						}*/
					}
				}
		}	}
		}
	}

	N_mu[0]   = 100*N_mu_BP;
	mu_len[0] = 100*mu_len_BP;

        N_mu[1]   = 100*N_mu_SVT1;
        mu_len[1] = 100*mu_len_SVT1;

        N_mu[2]   = 100*N_mu_SVT;
        mu_len[2] = 100*mu_len_SVT;

        N_mu[3]   = 100*N_mu_ST;
        mu_len[3] = 100*mu_len_ST;

        N_mu[4]   = 100*N_mu_DCH;
        mu_len[4] = 100*mu_len_DCH;

	for (int i = 0; i < 5; i++) {
		BP_len[i] = 100*BP_len[i];
	}
}

// Simply divides two histograms and calculates error
void Sim_Efficiency(TH1F *Efficiency, TH1F *True_All, TH1F *True_Detected, Int_t num_bins) {

	// Assume True_All has no error
        for (int i = 1; i <= num_bins; i++) {
		if (True_All->GetBinContent(i) > 0) {
	                Efficiency->SetBinContent(i, True_Detected->GetBinContent(i)/True_All->GetBinContent(i));
			Efficiency->SetBinError(i, TMath::Sqrt(True_Detected->GetBinContent(i))/True_All->GetBinContent(i));
		}
        }
}



// Corrects for efficiency and energy smearing, then scales the result to be the cross section
void Analyze_Segment(TH1F *Detected_Delta, TH1F *Corrected_Delta, TH1F *All_Delta, TH1F *Cross_Section, TH1F *Efficiency, TH2F *Energy_Smearing, Float_t mu_len, Float_t density, Float_t lower_bound, Float_t upper_bound, Float_t N_mu, Int_t num_bins) {

	double step_size = 5;
	TMatrix Transfer(52, 52);
	Invert_Matrix(Energy_Smearing, NULL, Transfer, 52, false); 
	
	// Normalization for region of detector
	Float_t norm = 1/(mu_len*density*step_size);
	Float_t norm_err = norm*TMath::Sqrt(1/N_mu); //TMath::Sqrt(97*N_mu)*(0.005/0.38);


/*
        // This is a new method I am testing using an "approximate unsmearing" technique
        TH1F *True_Dist = (TH1F*)Energy_Smearing->ProjectionX();
        TH1F *Meas_Dist = (TH1F*)Energy_Smearing->ProjectionY();
        True_Dist->Sumw2();
        Meas_Dist->Sumw2();

        True_Dist->Divide(Meas_Dist);

        // Divide the detected distribution by this function
        for (int i = 1; i < 61; i++) {
                Corrected_Delta->SetBinContent(i, Detected_Delta->GetBinContent(i));
                Corrected_Delta->SetBinError(i, Detected_Delta->GetBinError(i));
        }

        Corrected_Delta->Multiply(True_Dist);

        for (int j = 1; j < 61; j++) {
                All_Delta->SetBinContent(j, Corrected_Delta->GetBinContent(j));
                All_Delta->SetBinError(j, Corrected_Delta->GetBinError(j));
        }
        All_Delta->Divide(Efficiency);

        for (int j = 1; j < 61; j++) {
                Cross_Section->SetBinContent(j, All_Delta->GetBinContent(j)*norm);
                Cross_Section->SetBinError(j, All_Delta->GetBinError(j)*norm);
        }
*/



	// Matrix correction for energy smearing including the correction
	// and the error associated from the correction
	int energy_cut = 8;
        for (int row = 0; row < 60; row++) {
                Double_t weight = 0;
                Double_t Err_weight = 0;
		if (row > energy_cut) {
                	for (int col = 0; col < 60; col++) {
				if (col > energy_cut) {
	                        	weight = weight + Transfer[row - energy_cut][col - energy_cut]*Detected_Delta->GetBinContent(col + 1);
        	                	Err_weight = Err_weight + TMath::Power(Transfer[row - energy_cut][col - energy_cut], 2)*Detected_Delta->GetBinContent(col + 1);
				}
                	}
		}
		
                Corrected_Delta->SetBinContent(row + 1, weight);
		Corrected_Delta->SetBinError(row + 1, TMath::Sqrt(Err_weight));
	}
	
	for (int i = 1; i < 61; i++) {
		All_Delta->SetBinContent(i, Corrected_Delta->GetBinContent(i));
		All_Delta->SetBinError(i, Corrected_Delta->GetBinError(i));
	}
	All_Delta->Divide(Efficiency);

	
	for (int i = 1; i < 61; i++) {
		Cross_Section->SetBinContent(i, All_Delta->GetBinContent(i));
		Cross_Section->SetBinError(i, All_Delta->GetBinError(i));
	}
	Cross_Section->Scale(norm);
/*

	// Perform the Efficiency correction, and then scale the distribution by
	// norm to get the calculated cross section (and error).
	for (int i = 1; i <= 60; i++) {

			// Using matrix
			Float_t All_err = TMath::Sqrt(TMath::Power(Corrected_Delta->GetBinError(i)/Efficiency->GetBinContent(i + 10), 2) + TMath::Power(Efficiency->GetBinError(i + 10)*Corrected_Delta->GetBinContent(i)/TMath::Power(Efficiency->GetBinContent(i + 10),2),2));

			All_Delta->SetBinContent(i, Corrected_Delta->GetBinContent(i)/Efficiency->GetBinContent(i + 10));
			All_Delta->SetBinError(i, All_err);


			// not using matrix
			All_Delta->SetBinContent(i, Detected_Delta->GetBinContent(i)/Efficiency->GetBinContent(i + 1));
			Float_t All_err = TMath::Sqrt(TMath::Power(Detected_Delta->GetBinError(i + 1)/Efficiency->GetBinContent(i + 1), 2) + TMath::Power(Efficiency->GetBinError(i + 1)*Detected_Delta->GetBinContent(i + 1)/TMath::Power(Efficiency->GetBinContent(i + 1), 2), 2));
			All_Delta->SetBinError(i, All_err);
			

		if (Efficiency->GetBinContent(i) > 0 && Detected_Delta->GetBinContent(i) > 0) { 

			All_Delta->SetBinContent(i, Detected_Delta->GetBinContent(i)/Efficiency->GetBinContent(i));
			Float_t All_err = (Detected_Delta->GetBinContent(i)/Efficiency->GetBinContent(i))*TMath::Sqrt(TMath::Power(Detected_Delta->GetBinError(i)/Detected_Delta->GetBinContent(i), 2) + TMath::Power(Efficiency->GetBinError(i)/Efficiency->GetBinContent(i), 2));
			//All_Delta->SetBinError(i, All_err);

			Float_t Cross_Err = TMath::Sqrt(TMath::Power(All_err*norm, 2) + TMath::Power(All_Delta->GetBinContent(i)*norm_err, 2));
			Cross_Section->SetBinContent(i, All_Delta->GetBinContent(i)*norm);
			Cross_Section->SetBinError(i, Cross_Err);
*/		//}
	//	else {
	//		All_Delta->SetBinContent(i, 0);
	//		Cross_Section->SetBinContent(i, 0);
	//	}
	//}
}


// Similar to Analyze_Segment, however all inputs are 2d histograms, and the muon length number is now an array
void Analyze_BP(TH2F *Detected_Delta, TH2F *Corrected_Delta, TH2F *All_Delta, TH2F *Cross_Section, TH1F *Efficiency, TH2F *Energy_Smearing, Float_t mu_len[], Float_t density, Float_t lower_bound) {

	TMatrix Transfer(52, 52);
	Invert_Matrix(Energy_Smearing, NULL, Transfer, 52, false);
	
	for (int i = 0; i < 5; i++) {
	        Float_t norm = 1/(mu_len[i]*density*5000);
		for (int row = 0; row < 52; row++){
	                Double_t weight = 0;
        	        Double_t Err_weight = 0;
                	for (int col = 0; col < 52; col++) {
                        	weight = weight + Transfer[row][col]*Detected_Delta->GetBinContent(col + 9, i + 1);
	                        Err_weight = Err_weight + TMath::Power(Transfer[row][col], 2)*Detected_Delta->GetBinContent(col + 9, i + 1);
        	        }

                	Corrected_Delta->SetBinContent(row + 9, i + 1, weight);
               		Corrected_Delta->SetBinError(row + 9, i + 1, TMath::Sqrt(Err_weight));
		}
	

        	for (int j = 1; j < 61; j++) {
                	All_Delta->SetBinContent(j, i + 1, Corrected_Delta->GetBinContent(j));
	                All_Delta->SetBinError(j, i + 1, Corrected_Delta->GetBinError(j));
       		}
	        All_Delta->Divide(Efficiency);


        	for (int j = 1; j < 61; j++) {
                	Cross_Section->SetBinContent(j, i + 1, All_Delta->GetBinContent(j)*norm);
	                Cross_Section->SetBinError(j, i + 1, All_Delta->GetBinError(j)*norm);
        	}
	        //Cross_Section->Scale(norm);
	}

}


// Takes a TH2F and the address of an appropriately sized square matrix
// along with the size of the matrix, and computes the inverse matrix.
// perform_invert is included to test how TH2's get mapped to matrices
void Invert_Matrix(TH2F *h2, TH2F *Smearing_Correction, TMatrix &Inverted, int size, bool output_hist) {

        TMatrix Transfer(size, size);
        for (int col = 0; col < size; col++) {
                Double_t entries = 0;
                for (int row = 0; row < size; row++) { //for (int row = col - 4; row <= col + 4; row++) {
			//if (row > -1 && row < size) {
                        	Transfer[row][col] = h2->GetBinContent(col + 8, row + 8);
                        	entries = entries + h2->GetBinContent(col + 8, row + 8);
			//}
                }
                
		for (int row = 0; row < size; row++) {
                        Transfer[row][col] = (Transfer[row][col])/entries;
                }
        }
        Inverted = Transfer.Invert();

	if (output_hist == true) {
		for (int row = 0; row < size; row++) {
			for (int col = 0; col < size; col++) {
				Smearing_Correction->SetBinContent(col + 1, row + 1, Inverted[row][col]);
			}
		}
	}
}

// This is a quick monte carlo that tests the accuracy of the reconstruction
// process. It randomly samples a a random number of points from the true
// distribution then runs the correction on an ensemble of 1000 monte carlos.
void ToyMC(TH1F *efficiency, TH2F *energy_smearing, vector<TH1F*> hists, TH1F *True_Delta, TH1F *pull, int num_slices, double slice_size, double step_size) {

        TRandom3 *rando = new TRandom3(0);

	Double_t norm = TMath::Sqrt(2);
	Double_t num_points = True_Delta->Integral(1, num_bins_BP + 1);
	Double_t seen;
	TF1 *Dist = new TF1("Dist", TString::Format("TMath::Poisson(x, %g)", num_points), 0.8*num_points, 1.2*num_points);
	
	//TMatrix Transfer(num_bins_BP + 1, num_bins_BP + 1);
	//Invert_Matrix(energy_smearing, NULL, Transfer, num_bins_BP + 1, false);

                TH1F *Orig = new TH1F("Orig", "", num_bins_BP, 0, upper_bound_BP);
                TH1F *detected = new TH1F("detected", "", num_bins_BP, 0, upper_bound_BP);
                TH1F *All_Delta = new TH1F("All_Delta", "", num_bins_BP, 0, upper_bound_BP);
                TH1F *Residuals = new TH1F("Residuals", "", num_bins_BP, 0, upper_bound_BP);

        TH1F *True_Dist = (TH1F*)energy_smearing->ProjectionX();
        TH1F *Meas_Dist = (TH1F*)energy_smearing->ProjectionY();
        TH1F *Corrected_Delta = new TH1F("Corrected_Delta", "", 60, 0, 300);
        True_Dist->Sumw2();
        Meas_Dist->Sumw2();
        Corrected_Delta->Sumw2();

        True_Dist->Divide(Meas_Dist);

	Orig->Sumw2();
	detected->Sumw2();
	All_Delta->Sumw2();
	Residuals->Sumw2();

	for (int i = 0; i < 1000; i++) {
		if (i%100 == 0) cout << i << "/" << 1000 << endl;
		int entries = (int)(Dist->GetRandom());
		Orig->Reset();
		detected->Reset();
		All_Delta->Reset();
		Residuals->Reset();
		Corrected_Delta->Reset();
		for (int j = 0; j < entries; j++) {
			Float_t energy = True_Delta->GetRandom();
			Orig->Fill(energy);
			for (int k = 0; k < num_bins_BP; k++) {
	                        if (energy >= k*step_size && energy < (k + 1)*step_size) {
					seen = rando->Rndm();
					if (seen <= efficiency->GetBinContent(k + 1)) {
						for (int l = 0; l < num_slices; l++) {
							if (energy >= l*slice_size/norm && energy < (l + 1)*(slice_size)/norm) {
								detected->Fill(energy + (hists[l]->GetRandom())/norm);
							}
						}
					}
				}
			}
		}

        // This is a new method I am testing using an "approximate unsmearing" technique

        True_Dist->Divide(Meas_Dist);

        // Divide the detected distribution by this function
        for (int j = 1; j < 61; j++) {
                Corrected_Delta->SetBinContent(i, detected->GetBinContent(i));
                Corrected_Delta->SetBinError(i, detected->GetBinError(i));
        }

        Corrected_Delta->Multiply(True_Dist);

        for (int j = 1; j < 61; j++) {
                All_Delta->SetBinContent(j, Corrected_Delta->GetBinContent(j));
                All_Delta->SetBinError(j, Corrected_Delta->GetBinError(j));
        }
        All_Delta->Divide(efficiency);

        for (int row = 11; row < 61; row++) {

/*        	for (int row = 1; row < num_bins_BP + 1; row++) {
                	Double_t weight = 0;
	                Double_t Err_weight = 0;
        	        for (int col = 0; col < num_bins_BP + 1; col++) {
                	        weight = weight + Transfer[row][col]*detected->GetBinContent(col + 1);
                        	Err_weight = Err_weight + TMath::Power(Transfer[row][col], 2)*detected->GetBinContent(col + 1);
                	}

	                All_Delta->SetBinContent(row, weight/efficiency->GetBinContent(row + 1));
			Float_t All_err = TMath::Sqrt((Err_weight/(weight*weight)) + TMath::Power(efficiency->GetBinError(row + 1)/efficiency->GetBinContent(row + 1), 2));
			All_Delta->SetBinError(row, All_Delta->GetBinContent(row)*All_err);
*/
			Residuals->SetBinContent(row, All_Delta->GetBinContent(row) - Orig->GetBinContent(row + 1));
			Residuals->SetBinError(row, TMath::Sqrt(TMath::Power(All_Delta->GetBinError(row), 2) + TMath::Power(Orig->GetBinError(row + 1), 2)));
		
			pull->Fill(Residuals->GetBinContent(row)/Residuals->GetBinError(row));
		}
        }
	
	pull->Draw();
}


// Takes all of the created histograms and plots them
//   MC Truth ->  Green
//   MC Reco  ->  Blue
//   Real     ->  Red
void Create_Plots() {

	TFile *file = new TFile("results.root", "update");
	gStyle->SetOptStat(0);
	

	TH1F *True_KK2F_BP = (TH1F*)file->Get("True_All_BP");
	TH1F *Reco_KK2F_BP = (TH1F*)file->Get("MC_Cross_Section_BP");
	TH1F *Run3_BP      = (TH1F*)file->Get("Cross_Section_BP");

	// Beam Pipe Cross Sections
	TCanvas *c1 = new TCanvas("c1", "Cross Section BP", 700, 700);
		True_KK2F_BP->GetYaxis()->SetTitle("Cross Section (cm^{2} g^{-1} MeV^{-1})");
		True_KK2F_BP->GetYaxis()->SetTitleOffset(1.5);
		True_KK2F_BP->GetXaxis()->SetTitle("Energy (MeV)");
		True_KK2F_BP->SetLineColor(kGreen);
		True_KK2F_BP->Draw();
		Reco_KK2F_BP->SetLineColor(kBlue);
		Reco_KK2F_BP->Draw("same");
		Run3_BP->SetLineColor(kRed);
		Run3_BP->Draw("same");
	c1->BuildLegend();

	
	TH1F *MC_Detect_BP = (TH1F*)file->Get("MC_Detected_Delta_BP");
	TH1F *Run3_Detect_BP = (TH1F*)file->Get("Detected_Delta_BP");
	
	// Raw Distribution
	TCanvas *c2 = new TCanvas("c2", "Scaled raw distribution", 700, 700);
		Run3_Detect_BP->GetYaxis()->SetTitle("cm^{2} g^{-1} MeV^{-1}");
		Run3_Detect_BP->GetYaxis()->SetTitleOffset(1.4);
		Run3_Detect_BP->GetXaxis()->SetTitle("Reconstructed Energy (MeV)");
		Run3_Detect_BP->SetLineColor(kRed);
		Run3_Detect_BP->Draw();
		MC_Detect_BP->Draw("same");
	c2->BuildLegend();
	

	TH1F *Efficiency_BP = (TH1F*)file->Get("Efficiency_BP");

	// Efficiency
	TCanvas *c3 = new TCanvas("c3", "Delta Ray Detection Efficiency", 700, 700);
		Efficiency_BP->GetYaxis()->SetTitle("Efficiency");
		Efficiency_BP->GetYaxis()->SetTitleOffset(1.4);
		Efficiency_BP->GetXaxis()->SetTitle("Energy (MeV)");
		Efficiency_BP->Draw("h");
	

	TH1F* Energy_Smearing_BP = (TH1F*)file->Get("Energy_Smearing_BP");

	// Energy Smearing
	TCanvas *c4 = new TCanvas("c4", "Delta Ray Energy Smearing", 700, 700);
		Energy_Smearing_BP->GetYaxis()->SetTitle("Reconstructed Energy (MeV)");
		Energy_Smearing_BP->GetYaxis()->SetTitleOffset(1.3);
		Energy_Smearing_BP->GetXaxis()->SetTitle("True Energy (MeV)");
		Energy_Smearing_BP->Draw("colz");


	TH1F *Gamma_BP         = (TH1F*)file->Get("Gamma_BP");
	TH1F *MC_Gamma_BP      = (TH1F*)file->Get("MC_Gamma_BP");
	
	// Detected Pair Production
	TCanvas *c5 = new TCanvas("c5", "Gamma Conversion", 700, 700);
		Gamma_BP->GetYaxis()->SetTitle("cm^{2} g^{-1} MeV^{-1}");
		Gamma_BP->GetYaxis()->SetTitleOffset(1.4);
		Gamma_BP->GetXaxis()->SetTitle("Reconstructed Energy (MeV)");
		Gamma_BP->SetLineColor(kRed);
		Gamma_BP->Draw();
		MC_Gamma_BP->SetLineColor(kBlue);
		MC_Gamma_BP->Draw("same");
	c5->BuildLegend();

/*
	TH1F *Smearing_Correction_BP = (TH1F*)file->Get("Smearing_Correction_BP");

	TCanvas *c6 = new TCanvas("c6", "Smearing Correction", 700, 700);
		Smearing_Correction_BP->Draw("colz");
*/
/*
	TH1F *Residuals = new TH1F("Residuals", "", 60, 0, 300);

	// Plot the residuals between Reco KK2F and MC Truth
	for (int i = 1; i < 61; i++) {
		Residuals->SetBinContent(i, Reco_KK2F_BP->GetBinContent(i) - True_KK2F_BP->GetBinContent(i));
		Residuals->SetBinError(i, Reco_KK2F_BP->GetBinError(i));
	}
	TCanvas *c7 = new TCanvas("c7", "Residuals from reconstruction", 700, 700);
		Residuals->GetXaxis()->SetTitle("Energy (MeV)");
		Residuals->Draw();		
	Residuals->Write();*/
}

int BP_Offset() {

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);

	TChain *ch = new TChain("ntp667");
	for (int i = 0; i < 102; i++) {
		ch->AddFile(TString::Format("/data/HD4/babar/mupair/run3/AES_run3-BTM_3tracks-%d.root", i + 1));
	}
	
	TProfile* bprad = new TProfile("bprad","BP radius vs azimumth",50,-3.142,3.142,2.5,2.8,"i");
	ch->Project("bprad","sqrt(muBP0^2+muBP1^2):atan2(muBP1,muBP0)");

	TF1* sine = new TF1("sine","[0] + [1]*sin(x+[2])",3);
	sine->SetParameters(2.65,0.0,0.0);
	bprad->Fit(sine);

	TProfile* prad = new TProfile("prad","POCA radius vs azimumth;azimuth;radius (cm)",50,-3.142,3.142,2.5,2.8,"i");
	ntp667->Project("prad","sqrt(DeltaePOCA0^2+DeltaePOCA1^2):atan2(DeltaePOCA1,DeltaePOCA0)");
	sine->SetParameters(2.65, 0.0, 0.0);
	prad->Fit(sine);

	TCanvas *c1 = new TCanvas("c1", "", 600, 600);
	bprad->SetLineColor(kBlack);
	prad->Draw();
	//bprad->Draw();//"same");
}

void plot3region(){
	TFile* f[3];
	TGraphErrors* gr_theta13_Anorm[3];
	TGraph* gr_CL_sigma1[3];
	TGraph* gr_CL_sigma2[3];
	TCanvas* c = new TCanvas("c","c",800,600);
	c->cd();
	gPad->SetGrid();
	for(int r = 0;r<3;r++){
		f[r] = new TFile(TString::Format("rootFile2/region%i.info.root",r),"READ");
		gr_theta13_Anorm[r] = (TGraphErrors*)f[r]->Get("gr_2sin2theta13_Anorm");
		gr_CL_sigma1[r] = (TGraph*)f[r]->Get("CLsigma1");
		gr_CL_sigma2[r] = (TGraph*)f[r]->Get("CLsigma2");

		if(r==0){
			gr_theta13_Anorm[r]->SetMaximum(1.02);
			gr_theta13_Anorm[r]->SetMinimum(0.92);
			gr_theta13_Anorm[r]->GetXaxis()->SetLimits(0.0,0.12);
			gr_theta13_Anorm[r]->Draw("ap");
			gr_theta13_Anorm[r]->SetTitle("Confidence level contour sin^{2}2#theta_{13} versus A_{norm};sin^{2}2#theta_{13};A_{norm}");
			gr_theta13_Anorm[r]->GetXaxis()->CenterTitle();
			gr_theta13_Anorm[r]->GetYaxis()->CenterTitle();
		}
		else
			gr_theta13_Anorm[r]->Draw("pSAME");
		gr_CL_sigma1[r]->SetFillColorAlpha(2,0.5);
		gr_CL_sigma2[r]->SetFillColorAlpha(8,0.5);
		gr_CL_sigma2[r]->Draw("pfSAME");
		gr_CL_sigma1[r]->Draw("pfSAME");
		gr_theta13_Anorm[r]->Draw("pSAME");
	}
    	TGraphErrors* grGd = new TGraphErrors();
    	grGd->SetMarkerColor(kRed);
    	grGd->SetMarkerStyle(21);
    	grGd->SetMarkerSize(2.0);
    	grGd->SetPoint(0,0.0848,0.9538);
    	grGd->SetPointError(0,0.0044,0.0237);
    	grGd->SetLineWidth(2);
    	grGd->SetLineColor(kRed);
    	TGraphErrors* grH = new TGraphErrors();
    	grH->SetMarkerColor(4);
    	grH->SetMarkerStyle(21);
    	grH->SetMarkerSize(2.0);
    	grH->SetPoint(0,0.071,0.962);
    	grH->SetPointError(0,0.01,0.02);
    	grH->SetLineWidth(2);
    	grH->SetLineColor(4);
	grGd->Draw("pSAME");
	grH->Draw("pSAME");
	TLegend* lg = new TLegend(0.1,0.1,0.3,0.3);
	lg->AddEntry(grGd,"PRD nGd","p");
	lg->AddEntry(grH,"PRD nH","p");
	lg->AddEntry(gr_CL_sigma1[0],"1#sigma CL (68%)","f");
	lg->AddEntry(gr_CL_sigma2[0],"2#sigma CL (95%)","f");
	lg->Draw("SAME");
}

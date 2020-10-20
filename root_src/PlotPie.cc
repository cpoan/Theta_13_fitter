void PlotPie(){
    ifstream inTxt;
    inTxt.open("unifiedErrors.txt");
    int N = 4;
    double error;
    double errors[4] = {0,0,0,0};
    int color[4] = {2,3,4,5};
    for(int idx = 0;idx<15;idx++){
        inTxt>>error;
        if(idx==0)
            errors[0]+=error;
        else if(idx==11)
            errors[2]+=error;
        else
            errors[3]+=error;
    }
    double PureIBDErr = 0.123*errors[0];
    double AccErr = 0.875*errors[0];
    double OtherErr = 0.002*errors[0]; 
    errors[0] = PureIBDErr;
    errors[1] = AccErr;
    errors[3] +=OtherErr;

    const char* entryName[4] = {
        "Stat(Signal)",
        "#splitline{Stat(Backgrounds)}{Accidental pairs}",
        "#splitline{Uncorrelated Syst of}{detection efficiency}",
        "Others Syst"
    };
    TCanvas* c = new TCanvas("c","c",2000,1800);
    c->cd();
    TPie* p = new TPie("p","p",N,errors,color);
    p->SetTitle("Unified analysis: sin^{2}#theta_{13} error budgets");
    p->SetLabelFormat("#splitline{%txt}{(%perc)}");
    p->SetRadius(.32);
    p->SetTextSize(0.03);
    p->SetLabels(entryName);
    p->Draw(">");
    c->SaveAs("unifiedErrorPie.png");
}


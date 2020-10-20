void transfernH(){
    ofstream output;
    output.open("transfer.txt");
    cout<<fixed<<setprecision(3);
    output<<fixed<<setprecision(3);
    double eps_multi[8] = {
        0.9844,
        0.9845,
        0.9846,
        0.9844,
        0.9841,
        0.9839,
        0.9845,
        0.9846
    };
    double acc[8][2] ={
        {59.31,0.11},
        {58.34,0.11},
        {54.54,0.09},
        {55.07,0.04},
        {55.35,0.04},
        {59.27,0.04},
        {56.73,0.04},
        {52.71,0.09}
    };
    double lihe[3][2] = {
        {2.36,1.02},
        {1.73,0.75},
        {0.19,0.09}
    };
    double fast[3][2] = {
        {2.11,0.18},
        {1.81,0.17},
        {0.16,0.03}
    };
    double Amc6[8][2] = {
        {0.07,0.04},
        {0.07,0.04},
        {0.07,0.03},
        {0.03,0.02},
        {0.03,0.02},
        {0.03,0.02},
        {0.02,0.01},
        {0.07,0.03}
    };
    double Amc8[8][2] = {
        {0.07,0.04},
        {0.07,0.04},
        {0.07,0.03},
        {0.03,0.02},
        {0.03,0.02},
        {0.03,0.02},
        {0.02,0.01},
        {0.07,0.03}
    };
    for(int d = 0;d<8;d++){
        output<<acc[d][0]*eps_multi[d]<<"\t"<<acc[d][1]*eps_multi[d]<<"\t";
    }
    output<<"\n";
    for(int site = 0;site<3;site++){
        double multi = 0;
        if(site == 0){
            multi = eps_multi[0]+eps_multi[1];
            multi /= 2.;
            output<<lihe[site][0]*multi<<"\t"<<lihe[site][1]*multi<<"\t";
        }else if(site == 1){
            multi = eps_multi[2]+eps_multi[7];
            multi /= 2.;
            output<<lihe[site][0]*multi<<"\t"<<lihe[site][1]*multi<<"\t";
        }else{
            multi = eps_multi[3]+eps_multi[4]+eps_multi[5]+eps_multi[6];
            multi /= 4.;
            output<<lihe[site][0]*multi<<"\t"<<lihe[site][1]*multi<<"\t";
        }
    }
    output<<"\n";
    for(int site = 0;site<3;site++){
        double multi = 0;
        if(site == 0){
            multi = eps_multi[0]+eps_multi[1];
            multi /= 2.;
            output<<fast[site][0]*multi<<"\t"<<fast[site][1]*multi<<"\t";
        }else if(site == 1){
            multi = eps_multi[2]+eps_multi[7];
            multi /= 2.;
            output<<fast[site][0]*multi<<"\t"<<fast[site][1]*multi<<"\t";
        }else{
            multi = eps_multi[3]+eps_multi[4]+eps_multi[5]+eps_multi[6];
            multi /= 4.;
            output<<fast[site][0]*multi<<"\t"<<fast[site][1]*multi<<"\t";
        }
    }
    output<<"\n";
    for(int d = 0;d<8;d++){
        if(d!=6&&d!=7){
            output<<(217./1230*Amc6[d][0]+1013./1230*Amc8[d][0])*eps_multi[d]<<"\t"<<sqrt(pow(217./1230*Amc6[d][1],2)+pow(1013./1230*Amc8[d][1],2))<<"\t";
       }else{
            output<<(Amc8[d][0]*eps_multi[d])<<"\t"<<(Amc8[d][1]*eps_multi[d])<<"\t";
        }
    }
    output<<"\n";
}

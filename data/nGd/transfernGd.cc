void transfernGd(){
    ofstream output;
    output.open("transfer.txt");
    cout<<fixed<<setprecision(3);
    output<<fixed<<setprecision(3);
    double eps_multi[8] = {
        0.9744,
        0.9747,
        0.9757,
        0.9759,
        0.9758,
        0.9756,
        0.9758,
        0.9757
    };
    double acc[8][2] ={
        {8.46,0.09},
        {8.46,0.09},
        {6.29,0.06},
        {1.27,0.01},
        {1.19,0.01},
        {1.20,0.01},
        {0.98,0.01},
        {6.18,0.06}
    };
    double lihe[3][2] = {
        {2.46,1.06},
        {1.72,0.77},
        {0.15,0.06}
    };
    double fast[3][2] = {
        {0.79,0.10},
        {0.57,0.07},
        {0.05,0.01}
    };
    double Amc6[8][2] = {
        {0.27,0.12},
        {0.25,0.11},
        {0.28,0.13},
        {0.22,0.10},
        {0.21,0.10},
        {0.21,0.10},
        {0.00,0.00},
        {0.00,0.00}
    };
    double Amc8[8][2] = {
        {0.15,0.07},
        {0.16,0.07},
        {0.13,0.06},
        {0.04,0.02},
        {0.03,0.02},
        {0.03,0.02},
        {0.05,0.02},
        {0.15,0.07}
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

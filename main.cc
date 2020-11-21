#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <vector>
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPie.h"
#include "TPieSlice.h"
#include "TString.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/Interpolator.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "Fit/LogLikelihoodFCN.h"
#include "Fit/PoissonLikelihoodFCN.h"
#include "HFitInterface.h"
#include "TFitter.h"
#include "TMinuit.h"
#include "TLine.h"
#include "TMatrixD.h"
#include "TArrayD.h"
using namespace std;
TFile* outf;
int norm = 1;
int debug = 0;
int debug_average_power = 0;
int showfinal = 0;
int draw = 1;
bool firstRun = true;
const int Selection = 2;//choose nGd or nH or Unified
string selectionName[3] = {
    "nGd",
    "nH",
    "region0"
};
///////////////////////////////////////////Parameters of program needs
const int N_detectors = 8;
const int N_reactors = 6;
const int N_nuclears = 4;
const size_t energy_bins = 10;
const double energy_min = 2;
const double energy_max = 8;
const double energy_step = (energy_max - energy_min) / energy_bins;
///////////////////////////////////////////Layout data section
const double dist_map[N_detectors][N_reactors] = {
	{362.38, 371.76, 903.47, 817.16, 1353.62, 1265.32},
	{357.94, 368.41, 903.35, 816.90, 1354.23, 1265.89},
	{1332.48, 1358.15, 467.57, 489.58, 557.58, 499.21},
	{1919.63, 1894.34, 1533.18, 1533.63, 1551.38, 1524.94},
	{1917.52, 1891.98, 1534.92, 1535.03, 1554.77, 1528.05},
	{1925.26, 1899.86, 1538.93, 1539.47, 1556.34, 1530.08},
	{1923.15, 1897.51, 1540.67, 1540.87, 1559.72, 1533.18},
	{1337.43, 1362.88, 472.97, 495.35, 558.71, 501.07}
};
/////////////////////////////////////////// Reactor data section
void parseAveragePower(const char* data);
const double SystThermalPower = 0.005;
const char* week_average_beda_data = "./data/WeeklyAvg_P17B_by_Beda.txt";
double ThermalPower6AD[N_reactors] = {//GW
    2.082,
    2.874,
    2.516,
    2.554,
    2.825,
    1.976
};
double ThermalPower8AD[N_reactors] = {//GW
    2.514,
    2.447,
    2.566,
    2.519,
    2.519,
    2.550
};
double ThermalPower7AD[N_reactors] = {
    0,
    0,
    0,
    0,
    0,
    0
};
double ThermalPower_nH_period2[N_reactors] = {
    0,0,0,0,0,0
};
/////////////////////////////////////////////
const double SystActinideFraction = 0.05;
const double Fraction[N_nuclears] = {
    0.566,
    0.0755,
    0.30075,
    0.05775
};
const double FracCorr[N_nuclears][N_nuclears] = {
    { 1.00, -0.22, -0.53, -0.18},
    {-0.22,  1.00,  0.18,  0.26},
    {-0.53,  0.18,  1.00,  0.49},
    {-0.18,  0.26,  0.49,  1.00}
};
TMatrixD* corrMatrix;
double* invertCorrArray;
const double SystReleaseEnergy = 0.002;
const double releaseEnergyPerFission = 205.95;//MeV
const char *HuberMiller_data[N_nuclears] = {
	"./data/U235.dat",
	"./data/U238.dat",
	"./data/Pu239.dat",
	"./data/Pu241.dat"
};
const char* LBNL_Huber_data = "./data/LBNL_Huber_v0.txt";
const double SystPopu = 0.0015;
const double SystSNF = 0.0038;
//////////////////////////////////////////////////IBD candidates data section
const char* candidates_data[3] = {
    "./data/nGd/candidates.txt",
    "./data/nH/candidates.txt",
    "./data/unified/region0/candidates.txt"
};
void parseCandidates(const char* data);
double R_candidates[8][2];
//////////////////////////////////////////////////Background data section
const char* backgrounds_data[3] = {
    "./data/nGd/backgrounds.txt",
    "./data/nH/backgrounds.txt",
    "./data/unified/region0/backgrounds.final.txt"
};
void parseBackgrounds(const char* data);
double R_Accidentals[8][2];
double Syst_Accidentals[8];
double R_LiHe[3][2];
double Syst_LiHe[3];
double R_FastN[8][2];
double Syst_FastN[3];
double R_AmC[8][2];
double Syst_AmC;
//////////////////////////////////////////////////Efficiency section 
const char* efficiencies_data[3] = {
    "./data/nGd/efficiencies.txt",
    "./data/nH/efficiencies.txt",
    "./data/unified/region0/efficiencies.final.txt"
};
void parseEfficiencies(const char* data);
double eps_multi[8];
double Syst_eps_multi[8];
double eps_det,Syst_eps_det_corr,Syst_eps_det_uncorr;
const double effCorrection[3][6][8]={
    {
        {1.00084,1.00070,0.99941,0.99837,0.99698,0.99815,0.99675,1.00067},
        {1.00032,1.00028,0.99994,0.99639,0.99756,0.99842,0.99805,0.99998},
        {1.00144,1.00066,1.00013,0.99982,1.00047,0.99921,0.99973,1.00121},
        {1.00190,1.00039,1.00095,0.99912,1.00055,0.99937,1.00120,1.00087},
        {1.00038,1.00033,1.00177,0.99775,0.99937,0.99972,0.99812,1.00117},
        {1.00180,1.00145,1.00103,0.99797,1.00120,0.99972,0.99887,1.00115}
    },
    {
        {1.00052,1.00053,1.00093,0.99731,0.99711,0.99753,0.99684,1.00086},
        {1.00077,1.00092,1.00061,0.99727,0.99794,0.99783,0.99768,1.00074},
        {1.00210,1.00208,1.00109,0.99964,0.99968,1.00016,0.99970,1.00103},
        {1.00181,1.00216,1.00122,0.99939,0.99978,0.99982,0.99907,1.00101},
        {1.00028,1.00108,1.00162,0.99934,0.99935,0.99978,0.99933,1.00138},
        {1.00138,1.00149,1.00127,0.99989,1.00010,1.00014,0.99968,1.00098}
    },
    {
        {1.00199,1.00152,1.00130,0.99646,0.99930,0.99641,0.00599,1.00398},
        {1.00192,1.00122,1.00111,1.00105,0.99597,0.99435,0.99492,1.00222},
        {1.00443,1.00296,1.00211,1.00351,1.00211,1.00397,0.99899,1.00381},
        {1.00672,1.00472,1.00037,0.99779,1.00201,1.00269,0.99318,1.00103},
        {1.00573,1.00830,0.99951,1.00175,1.00319,1.00571,1.00055,1.00116},
        {1.00264,1.00643,1.00017,1.00236,1.00052,1.00055,1.00456,1.00019}
    }
};
double effCorrectionUnified[6][8];
//////////////////////////////////////////////////Targer Protons 
const char* targetMass_data = "./data/targetMass.txt";
void parseTargetMass(const char* data);
double Mass_GdLS[8];
double Syst_Mass_GdLS[8];
double Mass_LS[8];
double Syst_Mass_LS;
double Mass_Acrylic[8];
double Syst_Mass_Acrylic;
const double protonDensity_GdLS = 7.169e25;
const double Syst_protonDensity_GdLS = 0.0047;
const double protonDensity_LS = 7.116e25;
const double Syst_protonDensity_LS = 0.0060;
const double protonDensity_Acrylic = 4.780e25;
//////////////////////////////////////////////////Survival Probability 
const bool normal_order = true;

const double sin2_theta12[2] = {0.307, 0.013};
const double delta_m2_21[2] = {7.53e-5, 0.18e-5}; 
const double delta_m2_32[2] = {2.43e-3, 0.07e-3}; 
double d_ji(double d_m2_ji, double l, double e_nu);
double survival_prob(double sin_theta12, double sin2_2theta13, double d_21, double d_ee, double l, double e_nu);
//////////////////////////////////////////////////Function section 
double fcn(const double*);
double singleChi2(double x1,double x2,double sigma);
double corrSingleChi2(double x1,double mu1,double x2,double mu2,double sigma,double corr);
//void parseInput(double ibd_rate[8][2], double bkg_rate[8][2], double obs_rate[8][2]);
void parseCrossSection(ROOT::Math::Interpolator *&f, double &e_min, double &e_max);
void parseReactorSpectrum(const char *data);

ROOT::Math::Interpolator* nuclear_spectrum[N_nuclears];
ROOT::Math::Interpolator *cross_section;
double _cross_section[energy_bins];
double cs_min,cs_max;
double nullExpected[8];
double Expected[8];
double Observed[8][2];
double RateError[8];
//////////////////////////////////////////////////Fitting configuration 
int N_pars = 0;
void initialize_minimizer(ROOT::Math::Minimizer* mini);
void profile_minimizer(ROOT::Math::Minimizer* mini);
double array_sin2_2theta13[15][2];
//////////////////////////////////////////////////Draw working 
const double* final_pars;
const double* final_errors;
void DrawOscillationCurve(ROOT::Math::Minimizer* mini);
void DrawContour(ROOT::Math::Minimizer* mini);
void scanChi2(ROOT::Math::Minimizer*,unsigned int,unsigned int&,double*,double*,double,double);
double eff_prob(double,double,double);
double myEff_Prob(double* x,double* p){
    return eff_prob(p[0],x[0],p[1]);
}
void finalFCN(const double*);
double effSum(double L, int det);
double myEffSum(double* x,double* p){
    return effSum(x[0],p[0]);
}
double normalSum(int det);
double myNormalSum(double* x,double* p){
    return normalSum(p[0])+x[0]-x[0];
}
double oscSum(double L);
double myOsc(double* x,double* p){
    return oscSum(x[0]);
}

int main(){
	outf = new TFile("./region0.info.root","RECREATE");
    ///////////////////////////////////////////Correlation matrix
    corrMatrix = new TMatrixD(4,4);
    TArrayD data(16);
    for(int n1 = 0;n1<N_nuclears;n1++)
        for(int n2 = 0;n2<N_nuclears;n2++){
            data[n1*4+n2] = FracCorr[n1][n2];
        }
    corrMatrix->SetMatrixArray(data.GetArray());
    invertCorrArray = corrMatrix->Invert().GetMatrixArray();
    if(debug){
        //corrMatrix->Print();
        //corrMatrix->Invert().Print();
        for(int n1 = 0;n1<N_nuclears;n1++){
            for(int n2 = 0;n2<N_nuclears;n2++){
                cout<<" "<<invertCorrArray[n1*4+n2];
            }
            cout<<"\n";
        }
    }
    ///////////////////////////////////////////CrossSection
	parseCrossSection(cross_section, cs_min, cs_max);
	for(size_t bin = 0; bin < energy_bins; ++bin){
		double _e = ( 0.5 + bin ) * energy_step + energy_min;
		_cross_section[bin] = cross_section->Eval(_e);
	}
    //for(int n = 0;n<N_nuclears;n++)
    //    parseReactorSpectrum(HuberMiller_data[n],nuclear_spectrum[n]);
    parseAveragePower(week_average_beda_data);
    parseReactorSpectrum(LBNL_Huber_data);
    parseCandidates(candidates_data[Selection]);
    parseBackgrounds(backgrounds_data[Selection]);
    parseEfficiencies(efficiencies_data[Selection]);
    parseTargetMass(targetMass_data);
    double par[150];
    for(int i = 0;i<150;i++){
        if(i==101)
            par[i] = sin2_theta12[0];
        else if(i==102)
            par[i] = delta_m2_21[0];
        else if(i==103)
            par[i] = delta_m2_32[0];
        else if(i==104)
            par[i] = 0.0;
        else if(i==105)
            par[i] = 0.95; 
        else
            par[i] = 0.0;
    }
    cout<<"Chi2 = "<<fcn(par)<<"\n";
	ROOT::Math::Functor f_tmp(&fcn, N_pars);
	ROOT::Math::Minimizer *mini = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Minimize");
    if(!debug&&!debug_average_power){
        mini->SetFunction(f_tmp);
        initialize_minimizer(mini);
        mini->Minimize();
        mini->Hesse();
        final_pars = mini->X();
        final_errors = mini->Errors();
        if(draw){
            DrawOscillationCurve(mini);
            DrawContour(mini);
        }
        profile_minimizer(mini);
    }
    //TF1* test = new TF1("",myEff_Prob,0,2000,2);
    //test->SetParameters(2.43e-3,4);
    //TCanvas* ctest = new TCanvas("ctest","ctest",1600,900);
    //ctest->cd();
    //test->Draw();
    //ctest->SaveAs("./Figures/nGd/testEffProb.png");
    outf->Close();
};

double fcn(const double* par){
    double chi2 = 0;
    //////////////////////////////////////////////////////Reactor Shape proccesing
    const int powerPullIndexBegin = 0;
    const int fracPullIndexBegin = 7;
    const int populationPullIndexBegin=31;
    const int SNFPullIndexBegin=37;

    double fraction[N_reactors][N_nuclears];
    int fracCount = fracPullIndexBegin;
    for(int r = 0;r<N_reactors;r++){
        double sum = 0;
        for(int n = 0;n<N_nuclears;n++){
            fraction[r][n] = Fraction[n]*(1.+par[fracCount++]);        
            sum+=fraction[r][n];
        }
        for(int n = 0;n<N_nuclears;n++){
            fraction[r][n]/=sum;
            if(debug)
                cout<<"Frac"<<r<<"_"<<n<<" "<<fraction[r][n]<<"\n";
        }
    }

    double power[N_reactors][N_detectors];
    double NumberOfFission[N_reactors][N_detectors];
    double ReactorShape[N_reactors][N_detectors][energy_bins];
    for(int r = 0;r<6;r++){
        for(int d = 0;d<8;d++){
            if(Selection==0){
                if(d!=6&&d!=7){
                    power[r][d] = (ThermalPower6AD[r]*0.1756+ThermalPower8AD[r]*(1.-0.1756))*(1.+par[r+powerPullIndexBegin]);
                }else{
                    power[r][d] = ThermalPower8AD[r]*(1.+par[r+powerPullIndexBegin]);
                }
            }else if(Selection==2){
                if(d==0)
                    power[r][d] = (ThermalPower6AD[r]*217.+ThermalPower8AD[r]*1524.)*(1.+par[r+powerPullIndexBegin])/(1524.+217.);
                else if(d!=6&&d!=7){
                    power[r][d] = 0;
                    power[r][d]+=ThermalPower6AD[r]*217.;
                    power[r][d]+=ThermalPower8AD[r]*1524.;
                    power[r][d]+=ThermalPower7AD[r]*217.;
                    power[r][d]*=1/(1958.)*(1.+par[r+powerPullIndexBegin]);
                }else
                    power[r][d] = (ThermalPower7AD[r]*217.+ThermalPower8AD[r]*1524.)*(1.+par[r+powerPullIndexBegin])/(1524.+217.);
            }else if(Selection==1){
                if(d!=6&&d!=7){
                    power[r][d] = (ThermalPower6AD[r]*217.+ThermalPower_nH_period2[r]*404.)*(1.+par[r+powerPullIndexBegin])/(621);
                }else{
                    power[r][d] = ThermalPower_nH_period2[r]*(1.+par[r+powerPullIndexBegin]);
                }
            }
            if(debug)
                cout<<"Power Testing r"<<r+1<<" "<<d+1<<"="<<power[r][d]<<"\n";
            NumberOfFission[r][d] = power[r][d]*1.e3*1.e19/(releaseEnergyPerFission*1.60217*(1.+par[6]));
            if(debug)
                cout<<"NumburOfFission "<<r<<" "<<NumberOfFission[r][d]<<"\n";
            for(int bin = 0;bin<energy_bins;bin++){
		        double _e = ( 0.5 + bin ) * energy_step + energy_min;
                ReactorShape[r][d][bin] = 0;
                for(int n = 0;n<N_nuclears;n++){
                    ReactorShape[r][d][bin] += NumberOfFission[r][d] * fraction[r][n] * nuclear_spectrum[n]->Eval(_e) * energy_step;
                }
                ReactorShape[r][d][bin] *= (1.+par[r+populationPullIndexBegin]+par[r+SNFPullIndexBegin]);
            }
        }
    }
    //////////////////////////////////////////////////////Background calculating
    int AccPullIndexBegin = 43;
    int LiHePullIndexBegin = 51;
    int FastNPullIndexBegin = 54;
    int AmCPullIndexBegin = 57;

    double Bkg[8][2];
    for(int d = 0;d<8;d++){
        int site = 0;
        if(d==0||d==1)
            site = 0;
        else if(d==2||d==7)
            site = 1;
        else
            site = 2;
        Bkg[d][0] = 0;
        Bkg[d][1] = 0;
        Bkg[d][0]+=R_Accidentals[d][0]*(1+par[d+AccPullIndexBegin]);
        Bkg[d][0]+=R_LiHe[site][0]*(1+par[site+LiHePullIndexBegin]);
        Bkg[d][0]+=R_FastN[d][0]*(1+par[site+FastNPullIndexBegin]);
        Bkg[d][0]+=R_AmC[d][0]*(1+par[AmCPullIndexBegin]);
        
        Bkg[d][1]+=pow(R_Accidentals[d][1],2);
        Bkg[d][1]+=pow(R_LiHe[site][1],2);
        Bkg[d][1]+=pow(R_FastN[d][1],2);
        Bkg[d][1]+=pow(R_AmC[d][1],2);
        Bkg[d][1] = sqrt(Bkg[d][1]);
        if(debug)
            cout<<"Total backgrounds rate at AD"<<d+1<<" = "<<Bkg[d][0]<<"\n";
    }
    //////////////////////////////////////////////////////Efficiency calculating
    int epsCorPullIndexBegin = 58;
    int epsUncorrPullIndexBegin = 59;
    int epsMultiPullIndexBegin = 67;

    double eps[8];
    double eps_nH_SamePart[3] = {
        0.1570*0.9496*0.9467*0.9855*0.7504,
        0.9624*0.9123*0.6739*0.8456*0.7504,
        0.8459*0.1581*0.4667*0.8764*0.7504
    };
    double eps_nH[3][8];
    for(int d = 0;d<8;d++){
        eps[d] = (eps_det+par[epsCorPullIndexBegin]+par[d+epsUncorrPullIndexBegin])*(eps_multi[d]+par[d+epsMultiPullIndexBegin]);
        if(debug)
            cout<<"Total efficiency of AD"<<d+1<<" = "<<eps[d]<<"\n";
        if(Selection==1){
            for(int v = 0;v<3;v++)
                eps_nH[v][d] = eps_nH_SamePart[v]*(1.+par[epsCorPullIndexBegin])*(1.+par[d+epsUncorrPullIndexBegin])*(eps_multi[d]+par[d+epsMultiPullIndexBegin]);
        }
    }
    //////////////////////////////////////////////////////Target Protons Calculations
    int GdLSMassPullIndexBegin = 75;
    int LSMassPullIndexBegin = 83;
    int AcrylicMassPullIndexBegin = 91;
    int GdLSdensityPullIndexBegin = 99;
    int LSdensityPullIndexBegin = 100;
    
    double np[8];
    double np_nH[3][8];
    for(int d = 0;d<8;d++){
        np[d]   = (Mass_GdLS[d]+par[d+GdLSMassPullIndexBegin])*(protonDensity_GdLS*(1.+par[GdLSdensityPullIndexBegin]));
        if(Selection==0)
            continue;
        else if(Selection==2||Selection==1){
            np[d]   +=(Mass_LS[d]+par[d+LSMassPullIndexBegin])*(protonDensity_LS)*(1.+par[LSdensityPullIndexBegin]);
            np[d]   +=(Mass_Acrylic[d]*(1+par[d+AcrylicMassPullIndexBegin]))*protonDensity_Acrylic;
            np_nH[0][d] = (Mass_GdLS[d]+par[d+GdLSMassPullIndexBegin])*(protonDensity_GdLS*(1.+par[GdLSdensityPullIndexBegin]));
            np_nH[1][d] =(Mass_LS[d]+par[d+LSMassPullIndexBegin])*(protonDensity_LS)*(1.+par[LSdensityPullIndexBegin]);
            np_nH[2][d] =(Mass_Acrylic[d]*(1+par[d+AcrylicMassPullIndexBegin]))*protonDensity_Acrylic;
            for(int r = 0;r<6;r++){
                effCorrectionUnified[r][d]  =(effCorrection[0][r][d]*np_nH[0][d])/(np_nH[0][d]+np_nH[1][d]+np_nH[2][d]); 
                effCorrectionUnified[r][d]  +=(effCorrection[1][r][d]*np_nH[1][d])/(np_nH[0][d]+np_nH[1][d]+np_nH[2][d]); 
                effCorrectionUnified[r][d]  +=(effCorrection[2][r][d]*np_nH[2][d])/(np_nH[0][d]+np_nH[1][d]+np_nH[2][d]); 
            }
        }
    }
    if(debug)
        for(int d = 0;d<8;d++)
            cout<<"This is targetmass times proton density of AD"<<d+1<<" "<<np[d]<<"\n";
    //////////////////////////////////////////////////////Deal with oscillation parameter
    int sin2_theta12PullIndexBegin = 101;
    int deltaM2_21PullIndexBegin = 102;
    int deltaM2_32PullIndexBegin = 103;
    int sin2_2theta13Index = 104;
    double _sin2_theta12 = par[sin2_theta12PullIndexBegin];
    double _sin_theta12 = sqrt(_sin2_theta12);
    double _delta_m2_21 = par[deltaM2_21PullIndexBegin];
    double _delta_m2_32 = par[deltaM2_32PullIndexBegin];
    double _sin2_2theta13 = par[sin2_2theta13Index];
	double _delta_m2_ee = normal_order? _delta_m2_32 + 5.2e-5: _delta_m2_32 - 5.2e-5;
    //////////////////////////////////////////////////////Expected Rate calculating
    for(int d = 0;d<8;d++){
        Expected[d] = 0;
        for(int bin = 0;bin<energy_bins;bin++){
				double _e = ( 0.5 + bin ) * energy_step + energy_min;
            for(int r = 0;r<N_reactors;r++){
				double _prob = survival_prob(_sin_theta12, _sin2_2theta13, _delta_m2_21, _delta_m2_ee, dist_map[d][r], _e);
                if(Selection==0||Selection==2)
                    Expected[d]+=par[105]*ReactorShape[r][d][bin]*np[d]*_cross_section[bin]*1.e-4/(4*3.14*dist_map[d][r]*dist_map[d][r])*86400*eps[d]*_prob;//0.80 is efficiency
                else if(Selection==1){
                    for(int v = 0;v<3;v++)
                        Expected[d]+=par[105]*ReactorShape[r][d][bin]*_cross_section[bin]*1.e-4/(4*3.14*dist_map[d][r]*dist_map[d][r])*86400*np_nH[v][d]*eps_nH[v][d]*_prob*effCorrection[v][r][d];//0.80 is efficiency
                }
            }
        }
        if(debug)
            cout<<"Expected Daily Rate of AD"<<d+1<<" = "<<Expected[d]<<"\n";
    }
    for(int d = 0;d<8;d++){
        RateError[d] = sqrt(pow(Bkg[d][1],2)+pow(R_candidates[d][1],2));
        if(showfinal){
            cout<<"Effective survival probability at AD"<<d+1<<" = "<<(R_candidates[d][0]-Bkg[d][0])/Expected[d]<<"\n";
            cout<<"R_IBD = "<<R_candidates[d][0]-Bkg[d][0]<<", Expected = "<< Expected[d]<<"\n";
            //cout<<"Rate Error of AD"<<d+1<<" = "<<RateError[d]<<"\n";
        }
    }
    if(firstRun){
        firstRun = false;
        for(int ad = 0;ad<8;ad++)
            cout<<"Real IBD daily rates AD"<<ad+1<<" = "<<R_candidates[ad][0]-Bkg[ad][0]<<" +- "<<RateError[ad]<<"\n";
    }
    //////////////////////////////////////////////////////Chi2 calculating
    int count = 0;
    //Expected rate
    for(int d = 0;d<8;d++){
        chi2+=singleChi2(R_candidates[d][0]-Bkg[d][0],Expected[d],RateError[d]);
        //count++;
    }
    //Power pull term
    //int powerPullIndexBegin = 0;
    //int fracPullIndexBegin = 7;
    //int populationPullIndexBegin=31;
    //int SNFPullIndexBegin=37;
    for(int r = 0;r<6;r++){
        chi2+=singleChi2(par[r],0,SystThermalPower);
        count++;
    }
    
    //frac
    int fraccount = 0;
    for(int r = 0;r<6;r++){
        for(int n1 = 0;n1<4;n1++){
            for(int n2 = 0;n2<4;n2++){
                chi2+=corrSingleChi2(par[fracPullIndexBegin+r*4+n1],0,par[fracPullIndexBegin+r*4+n2],0,SystActinideFraction,invertCorrArray[n1*4+n2]);
                if(debug)
                    cout<<fracPullIndexBegin+r*4+n1<<" "<<fracPullIndexBegin+r*4+n2<<" "<<FracCorr[n1][n2]<<"\n";
            }
            count++;
        }
    }
    //popu
    //SNF
    for(int r = 0;r<6;r++){
        chi2+=singleChi2(par[r+populationPullIndexBegin],0,SystPopu);
        chi2+=singleChi2(par[r+SNFPullIndexBegin],0,SystSNF);
        count+=2;
    }
    //energyRelease
    chi2+=singleChi2(par[6],0,SystReleaseEnergy);
    count++;

    //int AccPullIndexBegin = 43;
    //int LiHePullIndexBegin = 51;
    //int FastNPullIndexBegin = 54;
    //int AmCPullIndexBegin = 57;
    //background
    for(int d = 0;d<8;d++){
        chi2+=singleChi2(par[d+AccPullIndexBegin],0,Syst_Accidentals[d]);
        count++;
    }
    for(int site = 0;site<3;site++){
        chi2+=singleChi2(par[site+LiHePullIndexBegin],0,Syst_LiHe[site]);
        chi2+=singleChi2(par[site+FastNPullIndexBegin],0,Syst_FastN[site]);
        count+=2;
    }
    chi2+=singleChi2(par[AmCPullIndexBegin],0,Syst_AmC);
    count++;
    //chi2+=singleChi2(par[AmCPullIndexBegin],0,Syst_AmC);
    //int epsCorPullIndexBegin = 58;
    //int epsUncorrPullIndexBegin = 59;
    //int epsMultiPullIndexBegin = 67;
    //efficiency
    chi2+=singleChi2(par[epsCorPullIndexBegin],0,Syst_eps_det_corr);
    count++;
    for(int d = 0;d<8;d++){
        chi2+=singleChi2(par[d+epsUncorrPullIndexBegin],0,Syst_eps_det_uncorr);
        chi2+=singleChi2(par[d+epsMultiPullIndexBegin],0,Syst_eps_multi[d]);
        count+=2;
    }
    //int GdLSMassPullIndexBegin = 75;
    //int LSMassPullIndexBegin = 83;
    //int AcrylicMassPullIndexBegin = 91;
    //int GdLSdensityPullIndexBegin = 99;
    //int LSdensityPullIndexBegin = 100;
    //Target mass and density
    for(int d = 0;d<8;d++){
        chi2+=singleChi2(par[d+GdLSMassPullIndexBegin],0,Syst_Mass_GdLS[d]);
        chi2+=singleChi2(par[d+LSMassPullIndexBegin],0,Syst_Mass_LS);
        chi2+=singleChi2(par[d+AcrylicMassPullIndexBegin],0,Syst_Mass_Acrylic);
        count+=3;
    }
    chi2+=singleChi2(par[GdLSdensityPullIndexBegin],0,Syst_protonDensity_GdLS);
    chi2+=singleChi2(par[LSdensityPullIndexBegin],0,Syst_protonDensity_LS);
    count+=2;
    //int sin2_theta12PullIndexBegin = 101;
    //int deltaM2_21PullIndexBegin = 102;
    //int deltaM2_32PullIndexBegin = 103;
    //int sin2_2theta13Index = 104;
    //oscillation parameters
    chi2+=singleChi2(par[sin2_theta12PullIndexBegin],sin2_theta12[0],sin2_theta12[1]);
    chi2+=singleChi2(par[deltaM2_21PullIndexBegin],delta_m2_21[0],delta_m2_21[1]);
    chi2+=singleChi2(par[deltaM2_32PullIndexBegin],delta_m2_32[0],delta_m2_32[1]);
    count+=3;
    count++;
    N_pars = count;
    if(norm)
        N_pars++;
    if(debug)
        cout<<"Total parameters = "<<N_pars<<"\n";

    return chi2;
}

void parseReactorSpectrum(const char *data){
	cout << "PARSING REACTOR SPECTRUM: " << data << endl;
	ifstream ifile(data);
	string str;
	vector<double> x, y0, y1, y2, y3;
	while(getline(ifile, str)){
		if(str[0] == '#') 
			continue;

		istringstream tmp(str);
		
		double _x, _y0, _y1, _y2, _y3;
		tmp >> _x >> _y0 >> _y1 >> _y2 >> _y3;
        if(debug)
            cout << _x << " " << _y0 << " " << _y1 << " " << _y2 << " " << _y3 <<"\n";
	//	cout << _x << " " << _y << endl;			
		x.push_back(_x);
		y0.push_back(_y0);
		y1.push_back(_y1);
		y2.push_back(_y2);
		y3.push_back(_y3);
	}
	nuclear_spectrum[0] = new ROOT::Math::Interpolator(x, y0, ROOT::Math::Interpolation::kLINEAR);
	nuclear_spectrum[1] = new ROOT::Math::Interpolator(x, y1, ROOT::Math::Interpolation::kLINEAR);
	nuclear_spectrum[2] = new ROOT::Math::Interpolator(x, y2, ROOT::Math::Interpolation::kLINEAR);
	nuclear_spectrum[3] = new ROOT::Math::Interpolator(x, y3, ROOT::Math::Interpolation::kLINEAR);
}

void parseCrossSection(ROOT::Math::Interpolator *&f, double &e_min, double &e_max){
	char *cross_section_data = "./data/cross_section";
	cout << "PARSING CROSS SECTION: " << cross_section_data << endl;
	ifstream ifile(cross_section_data);
	double _e, _cs;
	vector<double> x, y;
	while( ifile >> _e >> _cs ){
		x.push_back(_e);
		y.push_back(_cs);
		e_min = min(_e, e_min);
		e_max = max(_e, e_max);
	}
	f = new ROOT::Math::Interpolator(x, y, ROOT::Math::Interpolation::kLINEAR);
}

void parseCandidates(const char* data){
    cout<<"PARSING CANDIDATES\n"; 
    ifstream ifile(data);

    double eps_mu[8];
    double livetime[8];
    double N_candidates[8];

    for(int d = 0;d<8;d++)
        ifile>>N_candidates[d];
    for(int d = 0;d<8;d++)
        ifile>>livetime[d];
    for(int d = 0;d<8;d++)
        ifile>>eps_mu[d];
    for(int d = 0;d<8;d++){
        R_candidates[d][0] = N_candidates[d]/livetime[d]/eps_mu[d];
        R_candidates[d][1] = sqrt(N_candidates[d])/livetime[d]/eps_mu[d];
    }

    if(debug){
        for(int d = 0;d<8;d++)
            cout<<"IBD candidates of AD"<<d<<" "<<R_candidates[d][0]<<" "<<R_candidates[d][1]<<"\n";
    }

    cout<<"PARSING CANDIDATES OVER\n"; 
}
void parseBackgrounds(const char* data){
    cout<<"PARSING BACKGROUNDS\n";
    ifstream ifile(data);

    for(int d = 0;d<8;d++)
        ifile>>R_Accidentals[d][0]>>R_Accidentals[d][1];
    for(int d = 0;d<8;d++)
        ifile>>Syst_Accidentals[d];
    for(int site = 0;site<3;site++)
        ifile>>R_LiHe[site][0]>>R_LiHe[site][1];
    for(int site = 0;site<3;site++)
        ifile>>Syst_LiHe[site];
    for(int d = 0;d<8;d++)
        ifile>>R_FastN[d][0]>>R_FastN[d][1];
    for(int site = 0;site<3;site++)
        ifile>>Syst_FastN[site];
    for(int d = 0;d<8;d++)
        ifile>>R_AmC[d][0]>>R_AmC[d][1];
    ifile>>Syst_AmC;

    if(debug){
        for(int d = 0;d<8;d++)
            cout<<" "<<R_Accidentals[d][0]<<" "<<R_Accidentals[d][1];
        cout<<"\n";
        for(int d = 0;d<8;d++)
            cout<<" "<<Syst_Accidentals[d];
        cout<<"\n";
        for(int site = 0;site<3;site++)
            cout<<" "<<R_LiHe[site][0]<<" "<<R_LiHe[site][1];
        cout<<"\n";
        for(int site = 0;site<3;site++)
            cout<<" "<<Syst_LiHe[site];
        cout<<"\n";
        for(int d = 0;d<8;d++)
            cout<<" "<<R_FastN[d][0]<<" "<<R_FastN[d][1];
        cout<<"\n";
        for(int site = 0;site<3;site++)
            cout<<" "<<Syst_FastN[site];
        cout<<"\n";
        for(int d = 0;d<8;d++)
            cout<<" "<<R_AmC[d][0]<<" "<<R_AmC[d][1];
        cout<<"\n";
        cout<<Syst_AmC<<"\n";
    }

    cout<<"PARSING BACKGROUNDS OVER\n";
}

void parseEfficiencies(const char* data){
    cout<<"PARSING EFFICIENCIES\n";
    ifstream ifile(data);

    for(int d = 0;d<8;d++)
        ifile>>eps_multi[d];
    for(int d = 0;d<8;d++)
        ifile>>Syst_eps_multi[d];
    ifile>>eps_det>>Syst_eps_det_corr>>Syst_eps_det_uncorr;

    if(debug){
        for(int d = 0;d<8;d++)
            cout<<eps_multi[d]<<" ";
        cout<<"\n";
        for(int d = 0;d<8;d++)
            cout<<Syst_eps_multi[d]<<" ";
        cout<<"\n";
        cout<<eps_det<<Syst_eps_det_corr<<Syst_eps_det_uncorr<<"\n";
    }
    cout<<"PARSING EFFICIENCIES OVER\n";
}

void parseTargetMass(const char* data){
    cout<<"PARSING TARGET PROTONS\n";
    ifstream ifile(data);
    
    for(int d = 0;d<8;d++)
        ifile>>Mass_GdLS[d];
    for(int d = 0;d<8;d++)
        ifile>>Syst_Mass_GdLS[d];
    for(int d = 0;d<8;d++)
        ifile>>Mass_LS[d];
    ifile>>Syst_Mass_LS;
    for(int d = 0;d<8;d++)
        ifile>>Mass_Acrylic[d];
    ifile>>Syst_Mass_Acrylic;

    if(debug){
        for(int d = 0;d<8;d++)
            cout<<Mass_GdLS[d]<<" ";
        cout<<"\n";
        for(int d = 0;d<8;d++)
            cout<<Syst_Mass_GdLS[d]<<" ";
        cout<<"\n";
        for(int d = 0;d<8;d++)
            cout<<Mass_LS[d]<<" ";
        cout<<"\n";
        cout<<Syst_Mass_LS<<"\n";
        for(int d = 0;d<8;d++)
            cout<<Mass_Acrylic[d]<<" ";
        cout<<"\n";
        cout<<Syst_Mass_Acrylic<<"\n";
    }
    cout<<"PARSING TARGET PROTONS OVER\n";
}

double d_ji(double d_m2_ji, double l, double e_nu){
	return 1.267 * d_m2_ji * l / e_nu;
}

double survival_prob(double sin_theta12, double sin2_2theta13, double d_21, double d_ee, double l, double e_nu){
	double _theta13 = asin(sqrt(sin2_2theta13)) / 2;
    //if(isnan(_theta13)){
    //    _theta13=0;
    //}
	double _cos_theta13 = cos(_theta13);
	double _cos4_theta13 = _cos_theta13 * _cos_theta13 * _cos_theta13 * _cos_theta13;
	double _theta12 = asin(sin_theta12); 
	double _sin2_2theta12 = sin(2 * _theta12) * sin(2 * _theta12);
	double _sin2_delta_21 = sin(d_ji(d_21, l, e_nu)) * sin(d_ji(d_21, l, e_nu));
	double _sin2_delta_ee = sin(d_ji(d_ee, l, e_nu)) * sin(d_ji(d_ee, l, e_nu));
	//cout << _cos4_theta13 * _sin2_2theta12 * _sin2_delta_21 << " " << sin2_2theta13 * _sin2_delta_ee << " " << _cos4_theta13 * _sin2_2theta12 * _sin2_delta_21 / (sin2_2theta13 * _sin2_delta_ee) << endl;
	return 1 - _cos4_theta13 * _sin2_2theta12 * _sin2_delta_21 - sin2_2theta13 * _sin2_delta_ee;
}

double singleChi2(double x1,double x2,double sigma){
    return pow(x1-x2,2)/pow(sigma,2);
}

double corrSingleChi2(double x1,double mu1,double x2,double mu2,double sigma,double corr){
    return (x1-mu1)*(x2-mu2)/pow(sigma,2)*corr;
}

void initialize_minimizer(ROOT::Math::Minimizer *mini){
	const char *init_str = "Initializing parameter %2d %s %10.2e %10.2e\n";
    double step_ratio = 1.0e-3;

    mini->SetErrorDef(1.0);
    mini->SetTolerance(0.001);
    mini->SetPrintLevel(2);
    mini->SetStrategy(1);
    mini->SetMaxFunctionCalls(1000000);
    mini->SetMaxIterations(1000000);
    
    int idx = 0;
    string parname;
    for(int r = 0;r<6;r++){
        parname = TString::Format("Pull_Power_%i",r);
       mini->SetVariable(idx++,parname,0,SystThermalPower*step_ratio);  
    }
    mini->SetVariable(idx++,"Pull_ReleaseEnergyPerFission",0,SystReleaseEnergy*step_ratio);
    for(int r = 0;r<6;r++)
        for(int d = 0;d<4;d++){
            parname = TString::Format("Pull_Fraction_Reactor_Nuclear_%i_%i",r,d);
            mini->SetVariable(idx++,parname,0,SystActinideFraction*step_ratio);
        }
    for(int r = 0;r<6;r++){
        parname = TString::Format("Pull_population_ractor_%i",r);
        mini->SetVariable(idx++,parname,0,SystPopu*step_ratio);
    }
    for(int r = 0;r<6;r++){
        parname = TString::Format("Pull_SNF_%i",r);
        mini->SetVariable(idx++,parname,0,SystSNF*step_ratio);
    }
    for(int d = 0;d<8;d++){
        parname = TString::Format("Pull_Acc_%i",d);
        mini->SetVariable(idx++,parname,0,Syst_Accidentals[d]*step_ratio);
    }
    for(int site = 0;site<3;site++){
        parname = TString::Format("Pull_LiHe_%i",site);
        mini->SetVariable(idx++,parname,0,Syst_LiHe[site]*step_ratio);
    }
    for(int site = 0;site<3;site++){
        parname = TString::Format("Pull_FastN_%i",site);
        mini->SetVariable(idx++,parname,0,Syst_FastN[site]*step_ratio);
    }
    mini->SetVariable(idx++,"Pull_AmC",0,Syst_AmC*step_ratio);
    mini->SetVariable(idx++,"Pull_eps_corr",0,Syst_eps_det_corr*step_ratio);
    for(int d = 0;d<8;d++){
        parname = TString::Format("Pull_eps_uncorr_%i",d);
        mini->SetVariable(idx++,parname,0,Syst_eps_det_uncorr*step_ratio);
    }
    for(int d = 0;d<8;d++){
        parname = TString::Format("Pull_eps_multi_%i",d);
        mini->SetVariable(idx++,parname,0,Syst_eps_multi[d]*step_ratio);
    }
    for(int d = 0;d<8;d++){
        parname = TString::Format("Pull_TargetMass_GdLS_%i",d);
        mini->SetVariable(idx++,parname,0,Syst_Mass_GdLS[d]*step_ratio);
    }
    for(int d = 0;d<8;d++){
        parname = TString::Format("Pull_TargetMass_LS_%i",d);
        mini->SetVariable(idx++,parname,0,Syst_Mass_LS*step_ratio);
    }
    for(int d = 0;d<8;d++){
        parname = TString::Format("Pull_TargetMass_Acrylic_%i",d);
        mini->SetVariable(idx++,parname,0,Syst_Mass_Acrylic*step_ratio);
    }
    mini->SetVariable(idx++,"Pull_protonDensity_GdLS",0.,Syst_protonDensity_GdLS*step_ratio);
    mini->SetVariable(idx++,"Pull_protonDensity_LS",0.,Syst_protonDensity_LS*step_ratio);
    mini->SetVariable(idx++,"Pull_sin2_theta12",sin2_theta12[0],sin2_theta12[1]*step_ratio);
    mini->SetVariable(idx++,"Pull_delta_m2_21",delta_m2_21[0],delta_m2_21[1]*step_ratio);
    mini->SetVariable(idx++,"Pull_delta_m2_32",delta_m2_32[0],delta_m2_32[1]*step_ratio);
    mini->SetLowerLimitedVariable(idx++,"sin2_2theta13",0.086,0.0001,0);
    if(norm)
        mini->SetVariable(idx++,"NORM",0.95,0.0001);
    cout<<"Number of fit parameters = "<<idx<<"\n";

    if(debug){
        for(int i = 0;i<N_pars;i++)
            printf(init_str,i,mini->VariableName(i).c_str(),mini->X()[i], mini->Errors()[i]);
        cout<<"Number of free parameter in combined fit :"<< mini->NFree()<<"\n";
    }
}
void profile_minimizer(ROOT::Math::Minimizer *mini){
    int set = 0;
	const char *init_str = "Initializing parameter %2d %s %10.2e %10.2e\n";
    for(int idx = 0;idx<104;idx++)
        mini->FixVariable(idx);
    mini->Minimize();mini->Hesse();
    array_sin2_2theta13[set][0] = mini->X()[104];
    array_sin2_2theta13[set++][1] = mini->Errors()[104];
    for(int idx = 0;idx<104;idx++){
        //mini->SetVariable(idx,mini->VariableName(idx).c_str(),final_pars[idx],final_pars[idx]*step_ratio);
        mini->ReleaseVariable(idx);
        if(idx == 5){
            mini->Minimize();mini->Hesse();
            array_sin2_2theta13[set][0] = mini->X()[104];
            array_sin2_2theta13[set++][1] = mini->Errors()[104];
        }else if(idx == 6){
            mini->Minimize();mini->Hesse();
            array_sin2_2theta13[set][0] = mini->X()[104];
            array_sin2_2theta13[set++][1] = mini->Errors()[104];
        }else if(idx == 30){
            mini->Minimize();mini->Hesse();
            array_sin2_2theta13[set][0] = mini->X()[104];
            array_sin2_2theta13[set++][1] = mini->Errors()[104];
        }else if(idx == 36){
            mini->Minimize();mini->Hesse();
            array_sin2_2theta13[set][0] = mini->X()[104];
            array_sin2_2theta13[set++][1] = mini->Errors()[104];
        }else if(idx == 42){
            mini->Minimize();mini->Hesse();
            array_sin2_2theta13[set][0] = mini->X()[104];
            array_sin2_2theta13[set++][1] = mini->Errors()[104];
        }else if(idx == 50){
            mini->Minimize();mini->Hesse();
            array_sin2_2theta13[set][0] = mini->X()[104];
            array_sin2_2theta13[set++][1] = mini->Errors()[104];
        }else if(idx == 53){
            mini->Minimize();mini->Hesse();
            array_sin2_2theta13[set][0] = mini->X()[104];
            array_sin2_2theta13[set++][1] = mini->Errors()[104];
        }else if(idx == 56){
            mini->Minimize();mini->Hesse();
            array_sin2_2theta13[set][0] = mini->X()[104];
            array_sin2_2theta13[set++][1] = mini->Errors()[104];
        }else if(idx == 57){
            mini->Minimize();mini->Hesse();
            array_sin2_2theta13[set][0] = mini->X()[104];
            array_sin2_2theta13[set++][1] = mini->Errors()[104];
        }else if(idx == 58){
            mini->Minimize();mini->Hesse();
            array_sin2_2theta13[set][0] = mini->X()[104];
            array_sin2_2theta13[set++][1] = mini->Errors()[104];
        }else if(idx == 66){
            mini->Minimize();mini->Hesse();
            array_sin2_2theta13[set][0] = mini->X()[104];
            array_sin2_2theta13[set++][1] = mini->Errors()[104];
        }else if(idx == 74){
            mini->Minimize();mini->Hesse();
            array_sin2_2theta13[set][0] = mini->X()[104];
            array_sin2_2theta13[set++][1] = mini->Errors()[104];
        }else if(idx == 100){
            mini->Minimize();mini->Hesse();
            array_sin2_2theta13[set][0] = mini->X()[104];
            array_sin2_2theta13[set++][1] = mini->Errors()[104];
        }else if(idx == 103){
            mini->Minimize();mini->Hesse();
            array_sin2_2theta13[set][0] = mini->X()[104];
            array_sin2_2theta13[set++][1] = mini->Errors()[104];
        }
    }
    cout<<"Final Results!!\n";
    double errorBudget[set];
    double errorSum = 0;
    for(int s = 0;s<set;s++){
        cout<<"SET "<<s+1<<" sin2_2theta13 =  "<<array_sin2_2theta13[s][0]<<" +- "<<array_sin2_2theta13[s][1]<<"\n";
        if(s==0)
            errorBudget[s] = pow(array_sin2_2theta13[s][1],2);
        else
            errorBudget[s] =abs(pow(array_sin2_2theta13[s][1],2)-pow(array_sin2_2theta13[s-1][1],2));
        if(isnan(errorBudget[s]))
            errorBudget[s]=0.00001;
        errorSum+=errorBudget[s];
    }
    cout<<fixed<<setprecision(2);
    for(int s = 0;s<set;s++){
        cout<<"SET "<<s+1<< ",errorBudget(%)  = "<< errorBudget[s]/errorSum*100.<<"\n";
        errorBudget[s] = errorBudget[s]/errorSum*100.;
    }
    double sin22theta13 = array_sin2_2theta13[0][0];
    double totalError =  mini->Errors()[104];
    double statError = array_sin2_2theta13[0][1];
    double systError = sqrt(totalError*totalError-statError*statError);
    cout<<setprecision(3);
    cout<<"Final results of Sin2 2theta13 = "<<sin22theta13<<" +- "<<statError<<"(stat)"<<" +- "<<systError<<"(syst)\n";
    cout<<"Final results of Anorm(reactor neutrino anomaly) = " << mini->X()[105] <<" +- "<< mini->Errors()[105]<<"\n";

    ofstream fout;
    fout.open(TString::Format("./root_src/%sErrors.txt",selectionName[Selection].c_str()));
    for(int idx = 0;idx<set;idx++)
        fout<<errorBudget[idx]<<"\n";
    fout.close();

}

void parseAveragePower(const char* data){
    cout<<"PARSING AVERAGE POWER\n";
    ifstream ifile(data);
    double totalPower[3][6] = {
        {0,0,0,0,0,0},
        {0,0,0,0,0,0},
        {0,0,0,0,0,0}
    };
    int stageWeek[3] = {0,0,0};
    int week,core,start_sec,end_sec;
    double frac_power,nono,f0,f1,f2,f3;
    int end_6AD = 1343404800;
    int start_8AD = 1350576000;
    int end_8AD = 1482163200;
    int start_7AD = 1485360000;
    //For nH
    double totalPower_nH[6] = {0,0,0,0,0,0};
    int stageWeek_nH = 0;
    int end_nH_period2 =start_8AD+86400*(621-217);
    while(!ifile.eof()){
        ifile>>week>>core>>start_sec>>end_sec>>frac_power>>nono>>f0>>f1>>f2>>f3;
        if(debug_average_power){
            cout<<" "<<week;
            cout<<" "<<core;
            cout<<" "<<start_sec;
            cout<<" "<<end_sec;
            cout<<" "<<frac_power;
            cout<<" "<<nono;
            cout<<" "<<f0;
            cout<<" "<<f1;
            cout<<" "<<f2;
            cout<<" "<<f3;
            cout<<"\n";
        }
        if(start_sec<end_6AD){
            stageWeek[0]++;
            totalPower[0][core-1]+=frac_power;
        }else if(start_8AD<start_sec&&start_sec<end_8AD){
            stageWeek[1]++;
            totalPower[1][core-1]+=frac_power;
            if(start_sec<end_nH_period2){
                stageWeek_nH++;
                totalPower_nH[core-1]+=frac_power;
            }
        }else if(start_7AD<start_sec){
            stageWeek[2]++;
            totalPower[2][core-1]+=frac_power;
        }
    }
    cout<<"6AD average power\n";
    for(int r = 0;r<6;r++){
        cout<<" "<<totalPower[0][r]/stageWeek[0]*2.9*6;
    }
    cout<<"\n8AD average power\n";
    for(int r = 0;r<6;r++){
        cout<<" "<<totalPower[1][r]/stageWeek[1]*2.9*6;
    }
    cout<<"\n7AD average power\n";
    for(int r = 0;r<6;r++){
        cout<<" "<<totalPower[2][r]/stageWeek[2]*2.9*6;
        ThermalPower7AD[r] = totalPower[2][r]/stageWeek[2]*2.9*6;
    }
    cout<<"\nnH period 2 average power\n";
    for(int r= 0;r<6;r++){
        cout<<" "<<totalPower_nH[r]/stageWeek_nH*2.9*6;
        ThermalPower_nH_period2[r] = totalPower_nH[r]/stageWeek_nH*2.9*6;
    }
    //cout<<start_7AD<<"\n";
    //cout<<"\n";
    //cout<<stageWeek[2]<<" "<<totalPower[2][0]<<"\n";
    cout<<"\nPARSING AVERAGE POWER END\n";
}

void DrawOscillationCurve(ROOT::Math::Minimizer* mini){
    double Leff[8] = {
        510,
        516,
        522,
        1620,
        1630,
        1640,
        1650,
        528
    };
    finalFCN(mini->X());
    TF1* fsum[8];
    TF1* feff[8];
    TLegend* lg;
    lg = new TLegend(0.7,0.4,0.9,0.6);

    TCanvas* c1 = new TCanvas("c","c",1.05*1000,1.7*1000);
    c1->Divide(2,4);
    for(int d = 0;d<8;d++){
        fsum[d] = new TF1("",myNormalSum,0.0,3000,1);
        fsum[d]->SetParameter(0,d);
        feff[d] = new TF1(TString::Format("AD%i:F_{normal} and F_{effective};L_{eff};F",d+1),myEffSum,0.0,3000,1);
        feff[d]->SetParameter(0,d);
        c1->cd(d+1);
        gPad->SetGrid();
        fsum[d]->SetLineColor(kBlue);
        feff[d]->SetLineWidth(2);
        fsum[d]->SetLineWidth(2);
        feff[d]->GetXaxis()->CenterTitle(kTRUE);
        feff[d]->GetYaxis()->CenterTitle(kTRUE);
        feff[d]->GetXaxis()->SetTitleSize(0.04);
        feff[d]->GetYaxis()->SetTitleSize(0.04);
        feff[d]->Draw("");
        fsum[d]->Draw("SAME");
        if(d==0){
            lg->AddEntry(fsum[d],"F_{normal}","l");
            lg->AddEntry(feff[d],"F_{effective}","l");
        }
        lg->Draw("SAME");
    }
    c1->SaveAs(TString::Format("./Figures/%s/effectiveCurve.png",selectionName[Selection].c_str()));

    TCanvas* c2 = new TCanvas("c2","c2",1600,900);
    TGraphErrors* data[3];
    data[0] = new TGraphErrors(2);
    data[1] = new TGraphErrors(2);
    data[2] = new TGraphErrors(4);
    data[0]->SetPoint(1,510,Observed[0][0]/nullExpected[0]);
    data[0]->SetPointError(1,0,Observed[0][1]/nullExpected[0]);
    data[0]->SetPoint(2,516,Observed[1][0]/nullExpected[1]);
    data[0]->SetPointError(2,0,Observed[1][1]/nullExpected[1]);
    data[1]->SetPoint(1,545,Observed[2][0]/nullExpected[2]);
    data[1]->SetPointError(1,0,Observed[2][1]/nullExpected[2]);
    data[1]->SetPoint(2,551,Observed[7][0]/nullExpected[7]);
    data[1]->SetPointError(2,0,Observed[7][1]/nullExpected[7]);
    for(int d = 3;d<7;d++){
        data[2]->SetPoint(d-2,Leff[d],Observed[d][0]/nullExpected[d]);
        data[2]->SetPointError(d-2,0,Observed[d][1]/nullExpected[d]);
    }
    c2->cd();
    gPad->SetGrid();
    data[2]->SetMaximum(1.01);
    data[2]->SetMinimum(0.916);
    data[2]->GetXaxis()->SetLimits(0,2000);
    for(int site = 0;site<3;site++){
        data[site]->SetLineWidth(2);
        data[site]->SetMarkerStyle(8);
        data[site]->SetMarkerSize(1.5);
    }
    data[0]->SetLineColor(kBlue);
    data[1]->SetLineColor(kGreen);
    data[2]->SetLineColor(kRed);
    data[2]->SetTitle(TString::Format("Oscillation results of %s analysis;            Effective baseline L_{eff}(meter);R^{obs}/R^{pred}_{no_osc}",selectionName[Selection].c_str()));
    data[2]->GetXaxis()->CenterTitle(kTRUE);
    data[2]->GetYaxis()->CenterTitle(kTRUE);
    //data[0]->SetMarkerColor(kBlue);
    //data[1]->SetMarkerColor(kGreen);
    //data[2]->SetMarkerColor(kRed);
    data[2]->Draw("ap");
    data[1]->Draw("pSAME");
    data[0]->Draw("pSAME");
    data[2]->Draw("pSAME");

    outf->cd();
    data[0]->Write("gr_EH1");
    data[1]->Write("gr_EH2");
    data[2]->Write("gr_EH3");

    TF1* fnull = new TF1("","1.0",0,3000);
    TF1* fosc = new TF1("",myOsc,0,3000,0);
    
    fnull->SetLineColor(kBlue);
    fnull->SetLineWidth(2);
    fosc->SetLineWidth(2);
    fnull->Draw("SAME");
    fosc->Draw("SAME");
    fosc->Write("tf1_osc_curve");

    TLegend* lgOsc = new TLegend(0.12,0.12,0.4,0.4); 
    lgOsc->AddEntry(fosc,"Best fit oscillations","l");
    lgOsc->AddEntry(fnull,"No oscillations","l");
    lgOsc->AddEntry(data[0],"EH1","pl");
    lgOsc->AddEntry(data[1],"EH2","pl");
    lgOsc->AddEntry(data[2],"EH3","pl");
    lgOsc->Draw("SAME");
    c2->SaveAs(TString::Format("./Figures/%s/oscillationCurve.png",selectionName[Selection].c_str()));
    lgOsc->Write("lg_Osc");


}

void DrawContour(ROOT::Math::Minimizer* mini){
    double mean = final_pars[104];
    double sigma = final_errors[104];
    double mean2 = final_pars[105];
    double sigma2 = final_errors[105];
    TGraphErrors* grFinalPoint = new TGraphErrors();
    grFinalPoint->SetMarkerStyle(22);
    grFinalPoint->SetMarkerSize(2.0);
    grFinalPoint->SetPoint(0,mean,mean2);
    grFinalPoint->SetPointError(0,sigma,sigma2);
    grFinalPoint->SetLineWidth(2);
    grFinalPoint->SetLineColor(1);
    //////////////////////////////////////////
    TGraphErrors* grGd = new TGraphErrors();
    grGd->SetMarkerColor(7);
    grGd->SetMarkerStyle(21);
    grGd->SetMarkerSize(2.0);
    grGd->SetPoint(0,0.0848,0.9538);
    grGd->SetPointError(0,0.0044,0.0237);
    grGd->SetLineWidth(2);
    grGd->SetLineColor(7);
    //////////////////////////////////////////
    TGraphErrors* grH = new TGraphErrors();
    grH->SetMarkerColor(4);
    grH->SetMarkerStyle(21);
    grH->SetMarkerSize(2.0);
    grH->SetPoint(0,0.071,0.962);
    grH->SetPointError(0,0.01,0.02);
    grH->SetLineWidth(2);
    grH->SetLineColor(4);
    //////////////////////////////////////////
    unsigned int npoints = 100;
    double x1[3][npoints];
    double x2[3][npoints];
    TGraph* gr2d[3];
    mini->Contour(104,105,npoints,x1[0],x2[0]);
    mini->SetErrorDef(4.0);
    mini->Contour(104,105,npoints,x1[1],x2[1]);
    mini->SetErrorDef(9.0);
    mini->Contour(104,105,npoints,x1[2],x2[2]);
    gr2d[0] = new TGraph(npoints,x1[0],x2[0]);
    gr2d[1] = new TGraph(npoints,x1[1],x2[1]);
    gr2d[2] = new TGraph(npoints,x1[2],x2[2]);

    gr2d[2]->SetTitle(TString::Format("Confidence level contour of %s analysis;\\sin^{2}2\\theta_{13};A_{norm}",selectionName[Selection].c_str()));
    gr2d[2]->GetXaxis()->CenterTitle(kTRUE);
    gr2d[2]->GetYaxis()->CenterTitle(kTRUE);
    TLegend* lg = new TLegend(0.7,0.7,0.9,0.9);
    TCanvas* c1 = new TCanvas("c1","c1",1600,1000);
    c1->cd();
    TPad* pad1 = new TPad("pad1","pad1",0,0,0.7,0.7);
    TPad* pad2 = new TPad("pad2","pad2",0,0.7,0.7,1.0);
    TPad* pad3 = new TPad("pad3","pad3",0.7,0,1.0,0.7);
    pad1->cd();
    pad1->SetGrid();
    gr2d[0]->SetFillColorAlpha(2,0.5);
    gr2d[1]->SetFillColorAlpha(8,0.5);
    gr2d[2]->SetFillColorAlpha(4,0.5);
    gr2d[2]->GetXaxis()->SetLimits(mean-5*sigma,mean+5*sigma);
    gr2d[2]->SetMaximum(mean2+5*sigma2);
    gr2d[2]->SetMinimum(mean2-5*sigma2);
    gr2d[2]->Draw("alf");
    gr2d[1]->Draw("lfSAME");
    gr2d[0]->Draw("lfSAME");
    grFinalPoint->Draw("plSAME");

    outf->cd();
    grFinalPoint->Write("gr_2sin2theta13_Anorm");
    gr2d[0]->Write("CLsigma1");
    gr2d[1]->Write("CLsigma2");
    gr2d[2]->Write("CLsigma3");
    //lg->AddEntry(grFinalPoint,TString::Format("(%s)Best fit",selectionName[Selection].c_str()),"pl");
    lg->AddEntry(grFinalPoint,"Best fit","pl");
    lg->AddEntry(gr2d[0],"1 #sigma (68\% CL)","f");
    lg->AddEntry(gr2d[1],"2 #sigma (95\% CL)","f");
    lg->AddEntry(gr2d[2],"3 #sigma (99.7\% CL)","f");
    lg->Write("lg_CL");
    if(Selection==2||Selection==1){
        grGd->Draw("plSAME");
        lg->AddEntry(grGd,"previous nGd","pl");
        if(Selection==2){
            grH->Draw("plSAME");
            lg->AddEntry(grH,"previous nH","pl");
        }
    }
    lg->Draw("SAME");
    c1->cd();
    pad1->Draw();
    
    mini->SetErrorDef(1.0);
    double s2_2th13[npoints];
    double dChi[npoints];
    double A[npoints];
    double dChi2[npoints];
    //mini->Scan(104,npoints,s2_2th13,dChi,0.07,0.1);
    //mini->Scan(105,npoints,A,dChi2,0.88,1.03);
    scanChi2(mini,104,npoints,s2_2th13,dChi,mean-3*sigma,mean+3*sigma);
    scanChi2(mini,105,npoints,A,dChi2,mean2-3*sigma2,mean2+3*sigma2);
    //for(int p = 0;p<npoints;p++){
    //    dChi[p]/=6.;
    //    dChi2[p]/=6.;
    //}
    TGraph* gr1d[2];
    gr1d[0] = new TGraph(npoints,s2_2th13,dChi);
    gr1d[1] = new TGraph(npoints,dChi2,A);
    pad2->cd();
    pad2->SetGrid();
    gr1d[0]->SetTitle(";;\\Delta\\chi^{2}");
    gr1d[0]->GetYaxis()->CenterTitle(kTRUE);
    gr1d[0]->GetYaxis()->SetTitleSize(0.08);
    gr1d[0]->GetYaxis()->SetLabelSize(0.06);
    gr1d[0]->GetXaxis()->SetLabelOffset(999);
    gr1d[0]->GetXaxis()->SetLimits(mean-5*sigma,mean+5*sigma);
    gr1d[0]->SetMaximum(10);
    gr1d[0]->SetMinimum(0);
    gr1d[0]->SetLineWidth(2);
    gr1d[0]->Draw("aplY+");
    TLine* line_gr1d[2][3][3];
    for(int dev = 1;dev<=3;dev++){
        line_gr1d[0][dev-1][0] = new TLine(mini->X()[104]-dev*mini->Errors()[104],0,mini->X()[104]-dev*mini->Errors()[104],dev*dev);
        line_gr1d[0][dev-1][1] = new TLine(mini->X()[104]+dev*mini->Errors()[104],0,mini->X()[104]+dev*mini->Errors()[104],dev*dev);
        line_gr1d[0][dev-1][2] = new TLine(mini->X()[104]-dev*mini->Errors()[104],dev*dev,mini->X()[104]+dev*mini->Errors()[104],dev*dev);
        for(int line = 0;line<3;line++){
            line_gr1d[0][dev-1][line]->SetLineStyle(9);
            line_gr1d[0][dev-1][line]->SetLineColor(kBlue);
            line_gr1d[0][dev-1][line]->SetLineWidth(2);
            line_gr1d[0][dev-1][line]->Draw("SAME");
        }
    }
    c1->cd();
    pad2->Draw("SAME");
    pad3->cd();
    pad3->SetGrid();
    gr1d[1]->SetTitle(";\\Delta\\chi^{2};");
    gr1d[1]->GetXaxis()->CenterTitle(kTRUE);
    gr1d[1]->GetXaxis()->SetTitleSize(0.06);
    gr1d[1]->GetXaxis()->SetTitleOffset(0.7);
    gr1d[1]->GetYaxis()->SetLabelOffset(999);
    gr1d[1]->GetXaxis()->SetLimits(0,10);
    gr1d[1]->SetMaximum(mean2+5*sigma2);
    gr1d[1]->SetMinimum(mean2-5*sigma2);
    gr1d[1]->SetLineWidth(2);
    gr1d[1]->Draw("aplX+");
    for(int dev = 1;dev<=3;dev++){
        line_gr1d[1][dev-1][0] = new TLine(0,mini->X()[105]-dev*mini->Errors()[105],dev*dev,mini->X()[105]-dev*mini->Errors()[105]);
        line_gr1d[1][dev-1][1] = new TLine(0,mini->X()[105]+dev*mini->Errors()[105],dev*dev,mini->X()[105]+dev*mini->Errors()[105]);
        line_gr1d[1][dev-1][2] = new TLine(dev*dev,mini->X()[105]-dev*mini->Errors()[105],dev*dev,mini->X()[105]+dev*mini->Errors()[105]);
        for(int line = 0;line<3;line++){
            line_gr1d[1][dev-1][line]->SetLineStyle(9);
            line_gr1d[1][dev-1][line]->SetLineColor(kBlue);
            line_gr1d[1][dev-1][line]->SetLineWidth(2);
            line_gr1d[1][dev-1][line]->Draw("SAME");
        }
    }
    c1->cd();
    pad3->Draw("SAME");
    c1->SaveAs(TString::Format("./Figures/%s/contour.png",selectionName[Selection].c_str()));

    for(int i = 0;i<2;i++)
        for(int j = 0;j<3;j++)
            for(int k = 0;k<3;k++)
                delete line_gr1d[i][j][k];
    for(int i = 0;i<2;i++)
        delete gr1d[i];
    delete pad3;
    delete pad2;
    delete pad1;
    delete c1;
    for(int i = 0;i<3;i++)
        delete gr2d[i];
    delete grGd;
    delete grH;
    delete grFinalPoint;
};

void scanChi2(ROOT::Math::Minimizer *mini,unsigned int ivar,unsigned int& nstep,double* x,double* y,double xmin,double xmax){
    mini->FixVariable(ivar);
    mini->SetPrintLevel(1);

    double chi2_min = 0;

    if(Selection==0)
        chi2_min = 4.51;
    else if(Selection==1)
        chi2_min = 6.56;
    else if(Selection==2)
        chi2_min = 5.25;

    for(int s = 0;s<nstep;s++){
        x[s] = (xmax-xmin)/(nstep-1)*s+xmin;
        mini->SetVariableValue(ivar,x[s]);
        mini->SetStrategy(0);
        mini->Minimize();
        mini->SetStrategy(1);
        mini->Minimize();
        y[s] = mini->MinValue() - chi2_min;
    }
    mini->SetPrintLevel(1);
    for(int s = 0;s<nstep;s++)
        cout<<"Var "<<ivar<<", STEP "<<s<<", Chi2 = "<<y[s]<<"\n";
    mini->ReleaseVariable(ivar);
    mini->Minimize();
}

void finalFCN(const double* par){
    double chi2 = 0;
    //////////////////////////////////////////////////////Reactor Shape proccesing
    const int powerPullIndexBegin = 0;
    const int fracPullIndexBegin = 7;
    const int populationPullIndexBegin=31;
    const int SNFPullIndexBegin=37;

    double fraction[N_reactors][N_nuclears];
    int fracCount = fracPullIndexBegin;
    for(int r = 0;r<N_reactors;r++){
        double sum = 0;
        for(int n = 0;n<N_nuclears;n++){
            fraction[r][n] = Fraction[n]*(1.+par[fracCount++]);        
            sum+=fraction[r][n];
        }
        for(int n = 0;n<N_nuclears;n++){
            fraction[r][n]/=sum;
            if(debug)
                cout<<"Frac"<<r<<"_"<<n<<" "<<fraction[r][n]<<"\n";
        }
    }

    double power[N_reactors][N_detectors];
    double NumberOfFission[N_reactors][N_detectors];
    double ReactorShape[N_reactors][N_detectors][energy_bins];
    for(int r = 0;r<6;r++){
        for(int d = 0;d<8;d++){
            if(Selection==0){
                if(d!=6&&d!=7){
                    power[r][d] = (ThermalPower6AD[r]*0.1756+ThermalPower8AD[r]*(1.-0.1756))*(1.+par[r+powerPullIndexBegin]);
                }else{
                    power[r][d] = ThermalPower8AD[r]*(1.+par[r+powerPullIndexBegin]);
                }
            }else if(Selection==2){
                if(d==0)
                    power[r][d] = (ThermalPower6AD[r]*217.+ThermalPower8AD[r]*1524.)*(1.+par[r+powerPullIndexBegin])/(1524.+217.);
                else if(d!=6&&d!=7){
                    power[r][d] = 0;
                    power[r][d]+=ThermalPower6AD[r]*217.;
                    power[r][d]+=ThermalPower8AD[r]*1524.;
                    power[r][d]+=ThermalPower7AD[r]*217.;
                    power[r][d]*=1/(1958.)*(1.+par[r+powerPullIndexBegin]);
                }else
                    power[r][d] = (ThermalPower7AD[r]*217.+ThermalPower8AD[r]*1524.)*(1.+par[r+powerPullIndexBegin])/(1524.+217.);
            }else if(Selection==1){
                if(d!=6&&d!=7){
                    power[r][d] = (ThermalPower6AD[r]*217.+ThermalPower_nH_period2[r]*404.)*(1.+par[r+powerPullIndexBegin])/(621);
                }else{
                    power[r][d] = ThermalPower_nH_period2[r]*(1.+par[r+powerPullIndexBegin]);
                }
            }
            if(debug)
                cout<<"Power Testing r"<<r+1<<" "<<d+1<<"="<<power[r][d]<<"\n";
            NumberOfFission[r][d] = power[r][d]*1.e3*1.e19/(releaseEnergyPerFission*1.60217*(1.+par[6]));
            if(debug)
                cout<<"NumburOfFission "<<r<<" "<<NumberOfFission[r][d]<<"\n";
            for(int bin = 0;bin<energy_bins;bin++){
		        double _e = ( 0.5 + bin ) * energy_step + energy_min;
                ReactorShape[r][d][bin] = 0;
                for(int n = 0;n<N_nuclears;n++){
                    ReactorShape[r][d][bin] += NumberOfFission[r][d] * fraction[r][n] * nuclear_spectrum[n]->Eval(_e) * energy_step;
                }
                ReactorShape[r][d][bin] *= (1.+par[r+populationPullIndexBegin]+par[r+SNFPullIndexBegin]);
            }
        }
    }
    //////////////////////////////////////////////////////Background calculating
    int AccPullIndexBegin = 43;
    int LiHePullIndexBegin = 51;
    int FastNPullIndexBegin = 54;
    int AmCPullIndexBegin = 57;

    double Bkg[8][2];
    for(int d = 0;d<8;d++){
        int site = 0;
        if(d==0||d==1)
            site = 0;
        else if(d==2||d==7)
            site = 1;
        else
            site = 2;
        Bkg[d][0] = 0;
        Bkg[d][1] = 0;
        Bkg[d][0]+=R_Accidentals[d][0]*(1+par[d+AccPullIndexBegin]);
        Bkg[d][0]+=R_LiHe[site][0]*(1+par[site+LiHePullIndexBegin]);
        Bkg[d][0]+=R_FastN[d][0]*(1+par[site+FastNPullIndexBegin]);
        Bkg[d][0]+=R_AmC[d][0]*(1+par[AmCPullIndexBegin]);
        
        Bkg[d][1]+=pow(R_Accidentals[d][1],2);
        Bkg[d][1]+=pow(R_LiHe[site][1],2);
        Bkg[d][1]+=pow(R_FastN[d][1],2);
        Bkg[d][1]+=pow(R_AmC[d][1],2);
        Bkg[d][1] = sqrt(Bkg[d][1]);
        if(debug)
            cout<<"Total backgrounds rate at AD"<<d+1<<" = "<<Bkg[d][0]<<"\n";
    }
    //////////////////////////////////////////////////////Efficiency calculating
    int epsCorPullIndexBegin = 58;
    int epsUncorrPullIndexBegin = 59;
    int epsMultiPullIndexBegin = 67;

    double eps[8];
    double eps_nH_SamePart[3] = {
        0.1570*0.9496*0.9467*0.9855*0.7504,
        0.9624*0.9123*0.6739*0.8456*0.7504,
        0.8459*0.1581*0.4667*0.8764*0.7504
    };
    double eps_nH[3][8];
    for(int d = 0;d<8;d++){
        eps[d] = (eps_det+par[epsCorPullIndexBegin]+par[d+epsUncorrPullIndexBegin])*(eps_multi[d]+par[d+epsMultiPullIndexBegin]);
        if(debug)
            cout<<"Total efficiency of AD"<<d+1<<" = "<<eps[d]<<"\n";
        if(Selection==1){
            for(int v = 0;v<3;v++)
                eps_nH[v][d] = eps_nH_SamePart[v]*(1.+par[epsCorPullIndexBegin])*(1.+par[d+epsUncorrPullIndexBegin])*(eps_multi[d]+par[d+epsMultiPullIndexBegin]);
        }
    }
    //////////////////////////////////////////////////////Target Protons Calculations
    int GdLSMassPullIndexBegin = 75;
    int LSMassPullIndexBegin = 83;
    int AcrylicMassPullIndexBegin = 91;
    int GdLSdensityPullIndexBegin = 99;
    int LSdensityPullIndexBegin = 100;
    
    double np[8];
    double np_nH[3][8];
    for(int d = 0;d<8;d++){
        np[d]   = (Mass_GdLS[d]+par[d+GdLSMassPullIndexBegin])*(protonDensity_GdLS*(1.+par[GdLSdensityPullIndexBegin]));
        if(Selection==0)
            continue;
        else if(Selection==2){
            np[d]   +=(Mass_LS[d]+par[d+LSMassPullIndexBegin])*(protonDensity_LS)*(1.+par[LSdensityPullIndexBegin]);
            np[d]   +=(Mass_Acrylic[d]*(1+par[d+AcrylicMassPullIndexBegin]))*protonDensity_Acrylic;
        }else if(Selection==1){
            np_nH[0][d] = (Mass_GdLS[d]+par[d+GdLSMassPullIndexBegin])*(protonDensity_GdLS*(1.+par[GdLSdensityPullIndexBegin]));
            np_nH[1][d] =(Mass_LS[d]+par[d+LSMassPullIndexBegin])*(protonDensity_LS)*(1.+par[LSdensityPullIndexBegin]);
            np_nH[2][d] =(Mass_Acrylic[d]*(1+par[d+AcrylicMassPullIndexBegin]))*protonDensity_Acrylic;
        }
    }
    if(debug)
        for(int d = 0;d<8;d++)
            cout<<"This is targetmass times proton density of AD"<<d+1<<" "<<np[d]<<"\n";
    //////////////////////////////////////////////////////Deal with oscillation parameter
    int sin2_theta12PullIndexBegin = 101;
    int deltaM2_21PullIndexBegin = 102;
    int deltaM2_32PullIndexBegin = 103;
    int sin2_2theta13Index = 104;
    double _sin2_theta12 = par[sin2_theta12PullIndexBegin];
    double _sin_theta12 = sqrt(_sin2_theta12);
    double _delta_m2_21 = par[deltaM2_21PullIndexBegin];
    double _delta_m2_32 = par[deltaM2_32PullIndexBegin];
    double _sin2_2theta13 = par[sin2_2theta13Index];
	double _delta_m2_ee = normal_order? _delta_m2_32 + 5.2e-5: _delta_m2_32 - 5.2e-5;
    //////////////////////////////////////////////////////Expected Rate calculating
    for(int d = 0;d<8;d++){
        Expected[d] = 0;
        for(int bin = 0;bin<energy_bins;bin++){
				double _e = ( 0.5 + bin ) * energy_step + energy_min;
            for(int r = 0;r<N_reactors;r++){
				double _prob = survival_prob(_sin_theta12, _sin2_2theta13, _delta_m2_21, _delta_m2_ee, dist_map[d][r], _e);
				double null_prob = survival_prob(_sin_theta12, 0.0, _delta_m2_21, _delta_m2_ee, dist_map[d][r], _e);
                if(Selection==0||Selection==2){
                    Expected[d]+=par[105]*ReactorShape[r][d][bin]*np[d]*_cross_section[bin]*1.e-4/(4*3.14*dist_map[d][r]*dist_map[d][r])*86400*eps[d]*_prob;//0.80 is efficiency
                    nullExpected[d]+=par[105]*ReactorShape[r][d][bin]*np[d]*_cross_section[bin]*1.e-4/(4*3.14*dist_map[d][r]*dist_map[d][r])*86400*eps[d]*null_prob;//0.80 is efficiency
                }else if(Selection==1){
                    for(int v = 0;v<3;v++){
                        Expected[d]+=par[105]*ReactorShape[r][d][bin]*_cross_section[bin]*1.e-4/(4*3.14*dist_map[d][r]*dist_map[d][r])*86400*np_nH[v][d]*eps_nH[v][d]*_prob;//0.80 is efficiency
                        nullExpected[d]+=par[105]*ReactorShape[r][d][bin]*_cross_section[bin]*1.e-4/(4*3.14*dist_map[d][r]*dist_map[d][r])*86400*np_nH[v][d]*eps_nH[v][d]*null_prob;//0.80 is efficiency
                    }
                }
            }
        }
        //if(debug)
        //    cout<<"Expected Daily Rate of AD"<<d+1<<" = "<<Expected[d]<<"\n";
    }
    for(int d = 0;d<8;d++){
        RateError[d] = sqrt(pow(Bkg[d][1],2)+pow(R_candidates[d][1],2));
        Observed[d][0] = R_candidates[d][0]-Bkg[d][0];
        Observed[d][1] = RateError[d];
        cout<<"Final Results AD"<<d+1<<" = "<<100.*Observed[d][0]/nullExpected[d]<<"\n";
    }
}

double eff_prob(double d_ee, double l, double e_nu){
	double _sin2_delta_ee = sin(d_ji(d_ee, l, e_nu)) * sin(d_ji(d_ee, l, e_nu));

	return _sin2_delta_ee;
}

double normalSum(int det){
    const double * par = final_pars;

    const int powerPullIndexBegin = 0;
    const int fracPullIndexBegin = 7;
    const int populationPullIndexBegin=31;
    const int SNFPullIndexBegin=37;

    double fraction[N_reactors][N_nuclears];
    int fracCount = fracPullIndexBegin;
    for(int r = 0;r<N_reactors;r++){
        double sum = 0;
        for(int n = 0;n<N_nuclears;n++){
            fraction[r][n] = Fraction[n]*(1.+par[fracCount++]);        
            sum+=fraction[r][n];
        }
        for(int n = 0;n<N_nuclears;n++){
            fraction[r][n]/=sum;
            if(debug)
                cout<<"Frac"<<r<<"_"<<n<<" "<<fraction[r][n]<<"\n";
        }
    }

    double power[N_reactors][N_detectors];
    double NumberOfFission[N_reactors][N_detectors];
    double ReactorShape[N_reactors][N_detectors][energy_bins];
    for(int r = 0;r<6;r++){
        for(int d = 0;d<8;d++){
            if(Selection==0){
                if(d!=6&&d!=7){
                    power[r][d] = (ThermalPower6AD[r]*0.1756+ThermalPower8AD[r]*(1.-0.1756))*(1.+par[r+powerPullIndexBegin]);
                }else{
                    power[r][d] = ThermalPower8AD[r]*(1.+par[r+powerPullIndexBegin]);
                }
            }else if(Selection==2){
                if(d==0)
                    power[r][d] = (ThermalPower6AD[r]*217.+ThermalPower8AD[r]*1524.)*(1.+par[r+powerPullIndexBegin])/(1524.+217.);
                else if(d!=6&&d!=7){
                    power[r][d] = 0;
                    power[r][d]+=ThermalPower6AD[r]*217.;
                    power[r][d]+=ThermalPower8AD[r]*1524.;
                    power[r][d]+=ThermalPower7AD[r]*217.;
                    power[r][d]*=1/(1958.)*(1.+par[r+powerPullIndexBegin]);
                }else
                    power[r][d] = (ThermalPower7AD[r]*217.+ThermalPower8AD[r]*1524.)*(1.+par[r+powerPullIndexBegin])/(1524.+217.);
            }else if(Selection==1){
                if(d!=6&&d!=7){
                    power[r][d] = (ThermalPower6AD[r]*217.+ThermalPower_nH_period2[r]*404.)*(1.+par[r+powerPullIndexBegin])/(621);
                }else{
                    power[r][d] = ThermalPower_nH_period2[r]*(1.+par[r+powerPullIndexBegin]);
                }
            }
            if(debug)
                cout<<"Power Testing r"<<r+1<<" "<<d+1<<"="<<power[r][d]<<"\n";
            NumberOfFission[r][d] = power[r][d]*1.e3*1.e19/(releaseEnergyPerFission*1.60217*(1.+par[6]));
            if(debug)
                cout<<"NumburOfFission "<<r<<" "<<NumberOfFission[r][d]<<"\n";
            for(int bin = 0;bin<energy_bins;bin++){
		        double _e = ( 0.5 + bin ) * energy_step + energy_min;
                ReactorShape[r][d][bin] = 0;
                for(int n = 0;n<N_nuclears;n++){
                    ReactorShape[r][d][bin] += NumberOfFission[r][d] * fraction[r][n] * nuclear_spectrum[n]->Eval(_e) * energy_step;
                }
                ReactorShape[r][d][bin] *= (1.+par[r+populationPullIndexBegin]+par[r+SNFPullIndexBegin]);
            }
        }
    }
    //////////////////////////////////////////////////////Background calculating
    int AccPullIndexBegin = 43;
    int LiHePullIndexBegin = 51;
    int FastNPullIndexBegin = 54;
    int AmCPullIndexBegin = 57;

    double Bkg[8][2];
    for(int d = 0;d<8;d++){
        int site = 0;
        if(d==0||d==1)
            site = 0;
        else if(d==2||d==7)
            site = 1;
        else
            site = 2;
        Bkg[d][0] = 0;
        Bkg[d][1] = 0;
        Bkg[d][0]+=R_Accidentals[d][0]*(1+par[d+AccPullIndexBegin]);
        Bkg[d][0]+=R_LiHe[site][0]*(1+par[site+LiHePullIndexBegin]);
        Bkg[d][0]+=R_FastN[d][0]*(1+par[site+FastNPullIndexBegin]);
        Bkg[d][0]+=R_AmC[d][0]*(1+par[AmCPullIndexBegin]);
        
        Bkg[d][1]+=pow(R_Accidentals[d][1],2);
        Bkg[d][1]+=pow(R_LiHe[site][1],2);
        Bkg[d][1]+=pow(R_FastN[d][1],2);
        Bkg[d][1]+=pow(R_AmC[d][1],2);
        Bkg[d][1] = sqrt(Bkg[d][1]);
        if(debug)
            cout<<"Total backgrounds rate at AD"<<d+1<<" = "<<Bkg[d][0]<<"\n";
    }
    //////////////////////////////////////////////////////Efficiency calculating
    int epsCorPullIndexBegin = 58;
    int epsUncorrPullIndexBegin = 59;
    int epsMultiPullIndexBegin = 67;

    double eps[8];
    double eps_nH_SamePart[3] = {
        0.1570*0.9496*0.9467*0.9855*0.7504,
        0.9624*0.9123*0.6739*0.8456*0.7504,
        0.8459*0.1581*0.4667*0.8764*0.7504
    };
    double eps_nH[3][8];
    for(int d = 0;d<8;d++){
        eps[d] = (eps_det+par[epsCorPullIndexBegin]+par[d+epsUncorrPullIndexBegin])*(eps_multi[d]+par[d+epsMultiPullIndexBegin]);
        if(debug)
            cout<<"Total efficiency of AD"<<d+1<<" = "<<eps[d]<<"\n";
        if(Selection==1){
            for(int v = 0;v<3;v++)
                eps_nH[v][d] = eps_nH_SamePart[v]*(1.+par[epsCorPullIndexBegin])*(1.+par[d+epsUncorrPullIndexBegin])*(eps_multi[d]+par[d+epsMultiPullIndexBegin]);
        }
    }
    //////////////////////////////////////////////////////Target Protons Calculations
    int GdLSMassPullIndexBegin = 75;
    int LSMassPullIndexBegin = 83;
    int AcrylicMassPullIndexBegin = 91;
    int GdLSdensityPullIndexBegin = 99;
    int LSdensityPullIndexBegin = 100;
    
    double np[8];
    double np_nH[3][8];
    for(int d = 0;d<8;d++){
        np[d]   = (Mass_GdLS[d]+par[d+GdLSMassPullIndexBegin])*(protonDensity_GdLS*(1.+par[GdLSdensityPullIndexBegin]));
        if(Selection==0)
            continue;
        else if(Selection==2){
            np[d]   +=(Mass_LS[d]+par[d+LSMassPullIndexBegin])*(protonDensity_LS)*(1.+par[LSdensityPullIndexBegin]);
            np[d]   +=(Mass_Acrylic[d]*(1+par[d+AcrylicMassPullIndexBegin]))*protonDensity_Acrylic;
        }else if(Selection==1){
            np_nH[0][d] = (Mass_GdLS[d]+par[d+GdLSMassPullIndexBegin])*(protonDensity_GdLS*(1.+par[GdLSdensityPullIndexBegin]));
            np_nH[1][d] =(Mass_LS[d]+par[d+LSMassPullIndexBegin])*(protonDensity_LS)*(1.+par[LSdensityPullIndexBegin]);
            np_nH[2][d] =(Mass_Acrylic[d]*(1+par[d+AcrylicMassPullIndexBegin]))*protonDensity_Acrylic;
        }
    }
    if(debug)
        for(int d = 0;d<8;d++)
            cout<<"This is targetmass times proton density of AD"<<d+1<<" "<<np[d]<<"\n";
    //////////////////////////////////////////////////////Deal with oscillation parameter
    int sin2_theta12PullIndexBegin = 101;
    int deltaM2_21PullIndexBegin = 102;
    int deltaM2_32PullIndexBegin = 103;
    int sin2_2theta13Index = 104;
    double _sin2_theta12 = par[sin2_theta12PullIndexBegin];
    double _sin_theta12 = sqrt(_sin2_theta12);
    double _delta_m2_21 = par[deltaM2_21PullIndexBegin];
    double _delta_m2_32 = par[deltaM2_32PullIndexBegin];
    double _sin2_2theta13 = par[sin2_2theta13Index];
	double _delta_m2_ee = normal_order? _delta_m2_32 + 5.2e-5: _delta_m2_32 - 5.2e-5;
    //////////////////////////////////////////////////////Expected Rate calculating
    double sum = 0;
    for(int bin = 0;bin<energy_bins;bin++){
			double _e = ( 0.5 + bin ) * energy_step + energy_min;
        for(int r = 0;r<N_reactors;r++){
			double _prob = eff_prob(_delta_m2_ee, dist_map[det][r], _e);
            if(norm){
                sum+=par[105]*ReactorShape[r][det][bin]*np[det]*_cross_section[bin]*1.e-4/(4*3.14*dist_map[det][r]*dist_map[det][r])*86400*eps[det]*_prob;//0.80 is efficiency
            }
        }
    }
    return sum;
}

double effSum(double L, int det){
    const double * par = final_pars;

    const int powerPullIndexBegin = 0;
    const int fracPullIndexBegin = 7;
    const int populationPullIndexBegin=31;
    const int SNFPullIndexBegin=37;

    double fraction[N_reactors][N_nuclears];
    int fracCount = fracPullIndexBegin;
    for(int r = 0;r<N_reactors;r++){
        double sum = 0;
        for(int n = 0;n<N_nuclears;n++){
            fraction[r][n] = Fraction[n]*(1.+par[fracCount++]);        
            sum+=fraction[r][n];
        }
        for(int n = 0;n<N_nuclears;n++){
            fraction[r][n]/=sum;
            if(debug)
                cout<<"Frac"<<r<<"_"<<n<<" "<<fraction[r][n]<<"\n";
        }
    }

    double power[N_reactors][N_detectors];
    double NumberOfFission[N_reactors][N_detectors];
    double ReactorShape[N_reactors][N_detectors][energy_bins];
    for(int r = 0;r<6;r++){
        for(int d = 0;d<8;d++){
            if(Selection==0){
                if(d!=6&&d!=7){
                    power[r][d] = (ThermalPower6AD[r]*0.1756+ThermalPower8AD[r]*(1.-0.1756))*(1.+par[r+powerPullIndexBegin]);
                }else{
                    power[r][d] = ThermalPower8AD[r]*(1.+par[r+powerPullIndexBegin]);
                }
            }else if(Selection==2){
                if(d==0)
                    power[r][d] = (ThermalPower6AD[r]*217.+ThermalPower8AD[r]*1524.)*(1.+par[r+powerPullIndexBegin])/(1524.+217.);
                else if(d!=6&&d!=7){
                    power[r][d] = 0;
                    power[r][d]+=ThermalPower6AD[r]*217.;
                    power[r][d]+=ThermalPower8AD[r]*1524.;
                    power[r][d]+=ThermalPower7AD[r]*217.;
                    power[r][d]*=1/(1958.)*(1.+par[r+powerPullIndexBegin]);
                }else
                    power[r][d] = (ThermalPower7AD[r]*217.+ThermalPower8AD[r]*1524.)*(1.+par[r+powerPullIndexBegin])/(1524.+217.);
            }else if(Selection==1){
                if(d!=6&&d!=7){
                    power[r][d] = (ThermalPower6AD[r]*217.+ThermalPower_nH_period2[r]*404.)*(1.+par[r+powerPullIndexBegin])/(621);
                }else{
                    power[r][d] = ThermalPower_nH_period2[r]*(1.+par[r+powerPullIndexBegin]);
                }
            }
        }
    }
    //////////////////////////////////////////////////////Background calculating
    int AccPullIndexBegin = 43;
    int LiHePullIndexBegin = 51;
    int FastNPullIndexBegin = 54;
    int AmCPullIndexBegin = 57;

    double Bkg[8][2];
    for(int d = 0;d<8;d++){
        int site = 0;
        if(d==0||d==1)
            site = 0;
        else if(d==2||d==7)
            site = 1;
        else
            site = 2;
        Bkg[d][0] = 0;
        Bkg[d][1] = 0;
        Bkg[d][0]+=R_Accidentals[d][0]*(1+par[d+AccPullIndexBegin]);
        Bkg[d][0]+=R_LiHe[site][0]*(1+par[site+LiHePullIndexBegin]);
        Bkg[d][0]+=R_FastN[d][0]*(1+par[site+FastNPullIndexBegin]);
        Bkg[d][0]+=R_AmC[d][0]*(1+par[AmCPullIndexBegin]);
        
        Bkg[d][1]+=pow(R_Accidentals[d][1],2);
        Bkg[d][1]+=pow(R_LiHe[site][1],2);
        Bkg[d][1]+=pow(R_FastN[d][1],2);
        Bkg[d][1]+=pow(R_AmC[d][1],2);
        Bkg[d][1] = sqrt(Bkg[d][1]);
        if(debug)
            cout<<"Total backgrounds rate at AD"<<d+1<<" = "<<Bkg[d][0]<<"\n";
    }
    //////////////////////////////////////////////////////Efficiency calculating
    int epsCorPullIndexBegin = 58;
    int epsUncorrPullIndexBegin = 59;
    int epsMultiPullIndexBegin = 67;

    double eps[8];
    double eps_nH_SamePart[3] = {
        0.1570*0.9496*0.9467*0.9855*0.7504,
        0.9624*0.9123*0.6739*0.8456*0.7504,
        0.8459*0.1581*0.4667*0.8764*0.7504
    };
    double eps_nH[3][8];
    for(int d = 0;d<8;d++){
        eps[d] = (eps_det+par[epsCorPullIndexBegin]+par[d+epsUncorrPullIndexBegin])*(eps_multi[d]+par[d+epsMultiPullIndexBegin]);
        if(debug)
            cout<<"Total efficiency of AD"<<d+1<<" = "<<eps[d]<<"\n";
        if(Selection==1){
            for(int v = 0;v<3;v++)
                eps_nH[v][d] = eps_nH_SamePart[v]*(1.+par[epsCorPullIndexBegin])*(1.+par[d+epsUncorrPullIndexBegin])*(eps_multi[d]+par[d+epsMultiPullIndexBegin]);
        }
    }
    //////////////////////////////////////////////////////Target Protons Calculations
    int GdLSMassPullIndexBegin = 75;
    int LSMassPullIndexBegin = 83;
    int AcrylicMassPullIndexBegin = 91;
    int GdLSdensityPullIndexBegin = 99;
    int LSdensityPullIndexBegin = 100;
    
    double np[8];
    double np_nH[3][8];
    for(int d = 0;d<8;d++){
        np[d]   = (Mass_GdLS[d]+par[d+GdLSMassPullIndexBegin])*(protonDensity_GdLS*(1.+par[GdLSdensityPullIndexBegin]));
        if(Selection==0)
            continue;
        else if(Selection==2){
            np[d]   +=(Mass_LS[d]+par[d+LSMassPullIndexBegin])*(protonDensity_LS)*(1.+par[LSdensityPullIndexBegin]);
            np[d]   +=(Mass_Acrylic[d]*(1+par[d+AcrylicMassPullIndexBegin]))*protonDensity_Acrylic;
        }else if(Selection==1){
            np_nH[0][d] = (Mass_GdLS[d]+par[d+GdLSMassPullIndexBegin])*(protonDensity_GdLS*(1.+par[GdLSdensityPullIndexBegin]));
            np_nH[1][d] =(Mass_LS[d]+par[d+LSMassPullIndexBegin])*(protonDensity_LS)*(1.+par[LSdensityPullIndexBegin]);
            np_nH[2][d] =(Mass_Acrylic[d]*(1+par[d+AcrylicMassPullIndexBegin]))*protonDensity_Acrylic;
        }
    }
    if(debug)
        for(int d = 0;d<8;d++)
            cout<<"This is targetmass times proton density of AD"<<d+1<<" "<<np[d]<<"\n";
    //////////////////////////////////////////////////////Deal with oscillation parameter
    int sin2_theta12PullIndexBegin = 101;
    int deltaM2_21PullIndexBegin = 102;
    int deltaM2_32PullIndexBegin = 103;
    int sin2_2theta13Index = 104;
    double _sin2_theta12 = par[sin2_theta12PullIndexBegin];
    double _sin_theta12 = sqrt(_sin2_theta12);
    double _delta_m2_21 = par[deltaM2_21PullIndexBegin];
    double _delta_m2_32 = par[deltaM2_32PullIndexBegin];
    double _sin2_2theta13 = par[sin2_2theta13Index];
	double _delta_m2_ee = normal_order? _delta_m2_32 + 5.2e-5: _delta_m2_32 - 5.2e-5;
    //////////////////////////////////////////////////////Expected Rate calculating
    double sum = 0;
    for(int bin = 0;bin<energy_bins;bin++){
			double _e = ( 0.5 + bin ) * energy_step + energy_min;
        for(int r = 0;r<N_reactors;r++){
			double _prob = eff_prob(_delta_m2_ee, L, _e);
            if(Selection==0||Selection==2)
                sum+=par[105]*ReactorShape[r][det][bin]*np[det]*_cross_section[bin]*1.e-4/(4*3.14*dist_map[det][r]*dist_map[det][r])*86400*eps[det]*_prob;//0.80 is efficiency
            else if(Selection==1)
                for(int v = 0;v<3;v++)
                    sum+=par[105]*ReactorShape[r][det][bin]*_cross_section[bin]*1.e-4/(4*3.14*dist_map[det][r]*dist_map[det][r])*86400*np_nH[v][det]*eps_nH[v][det]*_prob;//0.80 is efficiency
        }
    }
    return sum;
}

double oscSum(double L){
    int det = 0;
    const double * par = final_pars;

    const int powerPullIndexBegin = 0;
    const int fracPullIndexBegin = 7;
    const int populationPullIndexBegin=31;
    const int SNFPullIndexBegin=37;

    double fraction[N_reactors][N_nuclears];
    int fracCount = fracPullIndexBegin;
    for(int r = 0;r<N_reactors;r++){
        double sum = 0;
        for(int n = 0;n<N_nuclears;n++){
            fraction[r][n] = Fraction[n]*(1.+par[fracCount++]);        
            sum+=fraction[r][n];
        }
        for(int n = 0;n<N_nuclears;n++){
            fraction[r][n]/=sum;
            if(debug)
                cout<<"Frac"<<r<<"_"<<n<<" "<<fraction[r][n]<<"\n";
        }
    }

    double power[N_reactors][N_detectors];
    double NumberOfFission[N_reactors][N_detectors];
    double ReactorShape[N_reactors][N_detectors][energy_bins];
    for(int r = 0;r<6;r++){
        for(int d = 0;d<8;d++){
            if(Selection==0){
                if(d!=6&&d!=7){
                    power[r][d] = (ThermalPower6AD[r]*0.1756+ThermalPower8AD[r]*(1.-0.1756))*(1.+par[r+powerPullIndexBegin]);
                }else{
                    power[r][d] = ThermalPower8AD[r]*(1.+par[r+powerPullIndexBegin]);
                }
            }else if(Selection==2){
                if(d==0)
                    power[r][d] = (ThermalPower6AD[r]*217.+ThermalPower8AD[r]*1524.)*(1.+par[r+powerPullIndexBegin])/(1524.+217.);
                else if(d!=6&&d!=7){
                    power[r][d] = 0;
                    power[r][d]+=ThermalPower6AD[r]*217.;
                    power[r][d]+=ThermalPower8AD[r]*1524.;
                    power[r][d]+=ThermalPower7AD[r]*217.;
                    power[r][d]*=1/(1958.)*(1.+par[r+powerPullIndexBegin]);
                }else
                    power[r][d] = (ThermalPower7AD[r]*217.+ThermalPower8AD[r]*1524.)*(1.+par[r+powerPullIndexBegin])/(1524.+217.);
            }else if(Selection==1){
                if(d!=6&&d!=7){
                    power[r][d] = (ThermalPower6AD[r]*217.+ThermalPower_nH_period2[r]*404.)*(1.+par[r+powerPullIndexBegin])/(621);
                }else{
                    power[r][d] = ThermalPower_nH_period2[r]*(1.+par[r+powerPullIndexBegin]);
                }
            }
            if(debug)
                cout<<"Power Testing r"<<r+1<<" "<<d+1<<"="<<power[r][d]<<"\n";
            NumberOfFission[r][d] = power[r][d]*1.e3*1.e19/(releaseEnergyPerFission*1.60217*(1.+par[6]));
            if(debug)
                cout<<"NumburOfFission "<<r<<" "<<NumberOfFission[r][d]<<"\n";
            for(int bin = 0;bin<energy_bins;bin++){
		        double _e = ( 0.5 + bin ) * energy_step + energy_min;
                ReactorShape[r][d][bin] = 0;
                for(int n = 0;n<N_nuclears;n++){
                    ReactorShape[r][d][bin] += NumberOfFission[r][d] * fraction[r][n] * nuclear_spectrum[n]->Eval(_e) * energy_step;
                }
                ReactorShape[r][d][bin] *= (1.+par[r+populationPullIndexBegin]+par[r+SNFPullIndexBegin]);
            }
        }
    }
    //////////////////////////////////////////////////////Background calculating
    int AccPullIndexBegin = 43;
    int LiHePullIndexBegin = 51;
    int FastNPullIndexBegin = 54;
    int AmCPullIndexBegin = 57;

    double Bkg[8][2];
    for(int d = 0;d<8;d++){
        int site = 0;
        if(d==0||d==1)
            site = 0;
        else if(d==2||d==7)
            site = 1;
        else
            site = 2;
        Bkg[d][0] = 0;
        Bkg[d][1] = 0;
        Bkg[d][0]+=R_Accidentals[d][0]*(1+par[d+AccPullIndexBegin]);
        Bkg[d][0]+=R_LiHe[site][0]*(1+par[site+LiHePullIndexBegin]);
        Bkg[d][0]+=R_FastN[d][0]*(1+par[site+FastNPullIndexBegin]);
        Bkg[d][0]+=R_AmC[d][0]*(1+par[AmCPullIndexBegin]);
        
        Bkg[d][1]+=pow(R_Accidentals[d][1],2);
        Bkg[d][1]+=pow(R_LiHe[site][1],2);
        Bkg[d][1]+=pow(R_FastN[d][1],2);
        Bkg[d][1]+=pow(R_AmC[d][1],2);
        Bkg[d][1] = sqrt(Bkg[d][1]);
        if(debug)
            cout<<"Total backgrounds rate at AD"<<d+1<<" = "<<Bkg[d][0]<<"\n";
    }
    //////////////////////////////////////////////////////Efficiency calculating
    int epsCorPullIndexBegin = 58;
    int epsUncorrPullIndexBegin = 59;
    int epsMultiPullIndexBegin = 67;

    double eps[8];
    double eps_nH_SamePart[3] = {
        0.1570*0.9496*0.9467*0.9855*0.7504,
        0.9624*0.9123*0.6739*0.8456*0.7504,
        0.8459*0.1581*0.4667*0.8764*0.7504
    };
    double eps_nH[3][8];
    for(int d = 0;d<8;d++){
        eps[d] = (eps_det+par[epsCorPullIndexBegin]+par[d+epsUncorrPullIndexBegin])*(eps_multi[d]+par[d+epsMultiPullIndexBegin]);
        if(debug)
            cout<<"Total efficiency of AD"<<d+1<<" = "<<eps[d]<<"\n";
        if(Selection==1){
            for(int v = 0;v<3;v++)
                eps_nH[v][d] = eps_nH_SamePart[v]*(1.+par[epsCorPullIndexBegin])*(1.+par[d+epsUncorrPullIndexBegin])*(eps_multi[d]+par[d+epsMultiPullIndexBegin]);
        }
    }
    //////////////////////////////////////////////////////Target Protons Calculations
    int GdLSMassPullIndexBegin = 75;
    int LSMassPullIndexBegin = 83;
    int AcrylicMassPullIndexBegin = 91;
    int GdLSdensityPullIndexBegin = 99;
    int LSdensityPullIndexBegin = 100;
    
    double np[8];
    double np_nH[3][8];
    for(int d = 0;d<8;d++){
        np[d]   = (Mass_GdLS[d]+par[d+GdLSMassPullIndexBegin])*(protonDensity_GdLS*(1.+par[GdLSdensityPullIndexBegin]));
        if(Selection==0)
            continue;
        else if(Selection==2){
            np[d]   +=(Mass_LS[d]+par[d+LSMassPullIndexBegin])*(protonDensity_LS)*(1.+par[LSdensityPullIndexBegin]);
            np[d]   +=(Mass_Acrylic[d]*(1+par[d+AcrylicMassPullIndexBegin]))*protonDensity_Acrylic;
        }else if(Selection==1){
            np_nH[0][d] = (Mass_GdLS[d]+par[d+GdLSMassPullIndexBegin])*(protonDensity_GdLS*(1.+par[GdLSdensityPullIndexBegin]));
            np_nH[1][d] =(Mass_LS[d]+par[d+LSMassPullIndexBegin])*(protonDensity_LS)*(1.+par[LSdensityPullIndexBegin]);
            np_nH[2][d] =(Mass_Acrylic[d]*(1+par[d+AcrylicMassPullIndexBegin]))*protonDensity_Acrylic;
        }
    }
    if(debug)
        for(int d = 0;d<8;d++)
            cout<<"This is targetmass times proton density of AD"<<d+1<<" "<<np[d]<<"\n";
    //////////////////////////////////////////////////////Deal with oscillation parameter
    int sin2_theta12PullIndexBegin = 101;
    int deltaM2_21PullIndexBegin = 102;
    int deltaM2_32PullIndexBegin = 103;
    int sin2_2theta13Index = 104;
    double _sin2_theta12 = par[sin2_theta12PullIndexBegin];
    double _sin_theta12 = sqrt(_sin2_theta12);
    double _delta_m2_21 = par[deltaM2_21PullIndexBegin];
    double _delta_m2_32 = par[deltaM2_32PullIndexBegin];
    double _sin2_2theta13 = par[sin2_2theta13Index];
	double _delta_m2_ee = normal_order? _delta_m2_32 + 5.2e-5: _delta_m2_32 - 5.2e-5;
    //////////////////////////////////////////////////////Expected Rate calculating
    double sum = 0;
    double sum2 = 0;
    for(int bin = 0;bin<energy_bins;bin++){
			double _e = ( 0.5 + bin ) * energy_step + energy_min;
        for(int r = 0;r<N_reactors;r++){
				double _prob = survival_prob(_sin_theta12, _sin2_2theta13, _delta_m2_21, _delta_m2_ee, L, _e);
				double _prob2 = survival_prob(_sin_theta12, 0.0, _delta_m2_21, _delta_m2_ee, L, _e);
            if(Selection==0||Selection==2){
                sum+=par[105]*ReactorShape[r][det][bin]*np[det]*_cross_section[bin]*1.e-4/(4*3.14*L*L)*86400*eps[det]*_prob;//0.80 is efficiency
                sum2+=par[105]*ReactorShape[r][det][bin]*np[det]*_cross_section[bin]*1.e-4/(4*3.14*L*L)*86400*eps[det]*_prob2;//0.80 is efficiency
            }else if(Selection==1){
                for(int v = 0;v<3;v++){
                    sum +=par[105]*ReactorShape[r][det][bin]*np_nH[v][det]*_cross_section[bin]*1.e-4/(4*3.14*L*L)*86400*eps_nH[v][det]*_prob;//0.80 is efficiency
                    sum2+=par[105]*ReactorShape[r][det][bin]*np_nH[v][det]*_cross_section[bin]*1.e-4/(4*3.14*L*L)*86400*eps_nH[v][det]*_prob2;//0.80 is efficiency

                }
            }
        }
    }
    return sum/sum2;
}

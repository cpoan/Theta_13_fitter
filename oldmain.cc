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
using namespace std;

bool debug = false;
bool draw_contour = false;

double eps_mult[8];
const double osc_par_inits[4][3] = {
	{8.60e-2, 0, 0.2},
	{3.07e-1, 0, 1.0},
	{7.53e-5, 0, 0},
	{2.43e-3, 0, 0}
};

const size_t npars = 24;

const int scan_pts = 30;
const double scan_dev = 5;

double minuit_print = 1;
double minuit_strategy = 1;

double simplex_args[2] = {1e5, 0.01};
double seek_args[2] = {1e4, 1.0};
double minimize_args[2] = {1e6, 1e-3};
const double minos_args[1] = {1e5};

const size_t energy_bins = 10;

const double energy_min = 2;
const double energy_max = 8;
const double energy_step = (energy_max - energy_min) / energy_bins;

const size_t reactor_nuclears = 3;
ROOT::Math::Interpolator *r_spec[reactor_nuclears];
ROOT::Math::Interpolator *cross_section;

double _cross_section[energy_bins];

double _weighted_reactor_power[8][6];

double ibd_rate[8][2], bkg_rate[8][2], obs_rate[8][2];

double cs_min, cs_max;

const size_t n_det = 8;
const size_t n_core = 6;

const double dist_map[n_det][n_core] = {
	{362.38, 371.76, 903.47, 817.16, 1353.62, 1265.32},
	{357.94, 368.41, 903.35, 816.90, 1354.23, 1265.89},
	{1332.48, 1358.15, 467.57, 489.58, 557.58, 499.21},
	{1337.43, 1362.88, 472.97, 495.35, 558.71, 501.07},
	{1919.63, 1894.34, 1533.18, 1533.63, 1551.38, 1524.94},
	{1917.52, 1891.98, 1534.92, 1535.03, 1554.77, 1528.05},
	{1925.26, 1899.86, 1538.93, 1539.47, 1556.34, 1530.08},
	{1923.15, 1897.51, 1540.67, 1540.87, 1559.72, 1533.18},
};

const double reactor_power_6ad[n_core] = {
	2.082, 2.874, 2.516, 2.554, 2.825, 1.976
};

const double reactor_power_8ad[n_core] = {
	2.514, 2.447, 2.566, 2.519, 2.519, 2.550
};

const double lt_6ad = 217;

const bool normal_order = true;

const double sin2_theta12[2] = {0.307, 0.013};
const double delta_m2_21[2] = {7.53e-5, 0.18e-5}; 
const double delta_m2_32[2] = {2.43e-3, 0.07e-3}; 


const char *reactor_data[reactor_nuclears] = {
	"../data/Pu239.dat",
	"../data/Pu241.dat",
	"../data/U235.dat"
};

const double det_eff[3] = {0.62165, 0.00545, 0.00523};

double pull(double x, double mean, double stddev);
double relu(double x);
double d_ji(double d_m2_ji, double l, double e_nu);
double bound_par(double x, double low, double up);
double survival_prob(double sin_theta12, double sin2_2theta13, double d_21, double d_ee, double l, double e_nu);
TGraph *minuit_profile(TFitter &minuit, int par_no, int npoints, double ndevs);
void parseInput(double ibd_rate[8][2], double bkg_rate[8][2], double obs_rate[8][2]);
void parseCrossSection(ROOT::Math::Interpolator *&f, double &e_min, double &e_max);
void parseReactorSpectrum(const char *data, ROOT::Math::Interpolator *&f);
double chi2_stats(
	double _sin2_2theta13,
	double _sin_theta12,
	double _delta_m2_21,
	double _delta_m2_ee,
	double *_reactor_shape,
	double _weighted_reactor_power[n_det][n_core],
	double *_cross_section,
	size_t energy_bins,
	double *_det_eff_uncorr,
	double _det_eff_corr,
	double *_bkg_uncorr,
	double _norm,
	double obs_rate[n_det][2],
	bool debug,
	bool print
	);
double chi2_syst(
	double _sin2_theta12,
	double _delta_m2_21,
	double _delta_m2_32,
	double *_det_eff_uncorr,
	double _det_eff_corr,
	double *_bkg_uncorr,
	bool print
	);
void parseParameters(
	double *par,
	double &_sin2_2theta13, 
	double &_sin2_theta12,
	double &_sin_theta12, 
	double &_delta_m2_21,
	double &_delta_m2_32, 
	double &_delta_m2_ee,
	double _reactor_shape[energy_bins],
	double _det_eff_uncorr[8],
	double &_det_eff_corr,
	double _bkg_uncorr[8],
	double &_norm);
void chi2_wrap(int &npar, double *g, double &result, double *par, int flag);
double chi2(double *par, bool print);
int main(int argc, char** argv){
    //double livetime[8] = {1544,2313,2323,2131,2320,2320,2320,2129};
    TCanvas* c1 = new TCanvas("c1","c1",1600,900);
	parseInput(ibd_rate, bkg_rate, obs_rate);		
	parseCrossSection(cross_section, cs_min, cs_max);
	for(size_t bin = 0; bin < energy_bins; ++bin){
		double _e = ( 0.5 + bin ) * energy_step + energy_min;
		_cross_section[bin] = cross_section->Eval(_e);
	}
    
	for(size_t r = 0; r < reactor_nuclears; ++r)
		parseReactorSpectrum(reactor_data[r], r_spec[r]);

	for(size_t det = 0; det < 8; ++det){
		for(size_t core = 0; core < 6; ++core){
            double a = 0;
			if(det == 3 || det == 7)
				_weighted_reactor_power[det][core] = reactor_power_8ad[core];
            else if(det == 0){
                a = lt_6ad*1./1123;
				_weighted_reactor_power[det][core] = a * reactor_power_6ad[core] + (1 - a) * reactor_power_8ad[core]; 
            }else{
                a = lt_6ad*1./1120;
				_weighted_reactor_power[det][core] = a * reactor_power_6ad[core] + (1 - a) * reactor_power_8ad[core]; 
            }
		}
	}


	TFitter minuit(npars);
	minuit.ExecuteCommand("SET PRINTOUT", &minuit_print, 1);
	minuit.ExecuteCommand("SET STRATEGY", &minuit_strategy, 1);
	minuit.SetFCN(chi2_wrap);

	int idx = 0;
	char buf[255];
	for(int p=0;p<4;++p){
		sprintf(buf, "OSC_PAR%d", p);
		minuit.SetParameter(idx++, buf, osc_par_inits[p][0], 0.01 * osc_par_inits[p][0], osc_par_inits[p][1], osc_par_inits[p][2]);
	}
	for(int p=0;p<reactor_nuclears-1;++p){
		sprintf(buf, "REACTOR_SHAPE%d", p);
		minuit.SetParameter(idx++, buf, 1, 0.1, 0, 100);
	}
	minuit.FixParameter(4);
	minuit.FixParameter(5);
	for(int p=0;p<8;++p){
		sprintf(buf, "DET_EFF_UNCORR%d", p);
		minuit.SetParameter(idx++, buf, 0, 1e-4, 0, 0);
	}
	minuit.SetParameter(idx++, "DET_EFF_CORR", 0, 1e-4, 0, 0);
	for(int p=0;p<8;++p){
		sprintf(buf, "BKG_UNCORR%d", p);
		minuit.SetParameter(idx++, buf, bkg_rate[p][0], bkg_rate[p][1], 0, 0);
	}
	minuit.SetParameter(idx++, "NORM", 4.90e5, 1e3, 0, 0);
    minuit.FixParameter(23);
	if(!debug){
		//minuit.ExecuteCommand("SEEK", seek_args, 2);
		//minuit.ExecuteCommand("SIMPLEX", simplex_args, 2);
		minuit.ExecuteCommand("MINIMIZE", minimize_args, 2);
	}else{
		double __tmp = 1;
		minuit.ExecuteCommand("SEEK", &__tmp, 1);
	}
	if(draw_contour){
		TGraph *g = (TGraph*)gMinuit->Contour(100, 0, 3);
		g->Draw();
		c1->SaveAs("test.png");
	}
	double pars[npars][2];
	for(size_t n=0;n<npars;++n){
		pars[n][0] = minuit.GetParameter(n);
		pars[n][1] = minuit.GetParError(n);
	}

	TGraph *g = minuit_profile(minuit, 0, scan_pts, scan_dev);
	g->SetMinimum(0);
	g->Draw();

	for(int dev=1; dev < 5; ++dev){
		int line_style = 9;
		double _c = pars[0][0];
		double _d = dev * pars[0][1];
		double _dev2 = dev * dev;
		TLine *_x1 = new TLine(_c - _d, 0, _c - _d, _dev2);
		_x1->SetLineStyle(line_style);
		TLine *_x2 = new TLine(_c + _d, 0, _c + _d, _dev2);
		_x2->SetLineStyle(line_style);
		TLine *_y1 = new TLine(_c - _d, _dev2, _c + _d, _dev2);
		_y1->SetLineStyle(line_style);
		_x1->Draw("same");
		_x2->Draw("same");
		_y1->Draw("same");	
	}

	c1->SaveAs("profile.png");
    return 1;
}

double pull(double x, double mean, double stddev){
	double tmp = ( x - mean ) / stddev;
	return tmp * tmp;
}

double relu(double x){
	return x > 0 ? x : 0;
}

double bound_par(double x, double low, double up){
	if(x < low)
		return low;
	else if(x > up)
		return up;
	else 
		return x;
}

double d_ji(double d_m2_ji, double l, double e_nu){
	return 1.267 * d_m2_ji * l / e_nu;
}

double survival_prob(double sin_theta12, double sin2_2theta13, double d_21, double d_ee, double l, double e_nu){
	double _theta13 = asin(sqrt(sin2_2theta13)) / 2;
    if(isnan(_theta13)){
        _theta13=0;
    }
	double _cos_theta13 = cos(_theta13);
	double _cos4_theta13 = _cos_theta13 * _cos_theta13 * _cos_theta13 * _cos_theta13;
	double _theta12 = asin(sin_theta12); 
	double _sin2_2theta12 = sin(2 * _theta12) * sin(2 * _theta12);
	double _sin2_delta_21 = sin(d_ji(d_21, l, e_nu)) * sin(d_ji(d_21, l, e_nu));
	double _sin2_delta_ee = sin(d_ji(d_ee, l, e_nu)) * sin(d_ji(d_ee, l, e_nu));
	//cout << _cos4_theta13 * _sin2_2theta12 * _sin2_delta_21 << " " << sin2_2theta13 * _sin2_delta_ee << " " << _cos4_theta13 * _sin2_2theta12 * _sin2_delta_21 / (sin2_2theta13 * _sin2_delta_ee) << endl;
	return 1 - _cos4_theta13 * _sin2_2theta12 * _sin2_delta_21 - sin2_2theta13 * _sin2_delta_ee;
}

TGraph *minuit_profile(TFitter &minuit, int par_no, int npoints, double ndevs){
	cout << "STARTING TO PROFILE PARAMETER" << par_no << endl;
	gMinuit->SetPrintLevel(-1);

	double _p0 = minuit.GetParameter(par_no);
	double _p0_dev = minuit.GetParError(par_no);	
	
	double amin, edm, errdef;
	int nvpar, nparx;
	minuit.GetStats(amin, edm, errdef, nvpar, nparx);

	double chi2_min = amin;

	minuit.FixParameter(par_no);
	vector<double> x, y;
	double scan_step = scan_dev / scan_pts;
	for(int pt = - (npoints - 1); pt < npoints; ++pt){
		cout << "\nSCANNING PAR" << par_no << " AT " << pt * scan_step << "DEV" << endl;
		double _tmp[2] = {par_no*1. + 1, _p0 + _p0_dev * scan_step * pt};
		minuit.ExecuteCommand("SET PARAMETER", _tmp, 2);
		minuit.ExecuteCommand("MINIMIZE", 0, 0);
		minuit.GetStats(amin, edm, errdef, nvpar, nparx);
		x.push_back(_tmp[1]);
		y.push_back(amin - chi2_min);
	}
	minuit.ReleaseParameter(par_no);	

	TGraph *g = new TGraph(x.size(), &x[0], &y[0]);
	return g;
}

void parseInput(double ibd_rate[8][2], double bkg_rate[8][2], double obs_rate[8][2]){
    ifstream inputText;
    inputText.open("../data/unified.p15a.modified.txt");
	double _nus[8];
	for(int i=0;i<8;++i) inputText >> _nus[i];
	double _lt[8];
	for(int i=0;i<8;++i) inputText >> _lt[i];
	double _eps_mu[8];
	for(int i=0;i<8;++i) inputText >> _eps_mu[i];
	double _eps_mult[8];
    for(int i=0;i<8;++i){
        inputText >> _eps_mult[i];
        eps_mult[i] = _eps_mult[i];
    }
	double _bkgs[6][8][2];
	for(int b=0;b<6;++b) 
		for(int i=0;i<8;++i){ 
			inputText >> _bkgs[b][i][0] >> _bkgs[b][i][1];
			bkg_rate[i][0] += _bkgs[b][i][0];
			bkg_rate[i][1] += ( _bkgs[b][i][1] * _bkgs[b][i][1] );
		}
	for(int i=0;i<8;++i){
		bkg_rate[i][1] = sqrt(bkg_rate[i][1]);
		cout << bkg_rate[i][0] << " " << bkg_rate[i][1] << " ";
		double _tmp = _lt[i] * _eps_mu[i];// * _eps_mult[i];
		obs_rate[i][0] = _nus[i]/ _tmp;
		obs_rate[i][1] = sqrt(_nus[i]) / _tmp;
		cout << obs_rate[i][0] << " " << obs_rate[i][1] << endl;
        cout << (obs_rate[i][0]-bkg_rate[i][0]) << endl;
	}
	for(int i=0;i<8;++i) inputText >> ibd_rate[i][0] >> ibd_rate[i][1];
}

void parseCrossSection(ROOT::Math::Interpolator *&f, double &e_min, double &e_max){
	char *cross_section_data = "../data/cross_section";
	cout << "PARSING CROSS SECTION: " << cross_section_data << endl;
	ifstream ifile(cross_section_data);
	double _e, _cs;
	vector<double> x, y;
	while( ifile >> _e >> _cs ){
		x.push_back(_e);
		y.push_back(_cs * 1e44);
		e_min = min(_e, e_min);
		e_max = max(_e, e_max);
	}
	f = new ROOT::Math::Interpolator(x, y, ROOT::Math::Interpolation::kLINEAR);
}

void parseReactorSpectrum(const char *data, ROOT::Math::Interpolator *&f){
	cout << "PARSING REACTOR SPECTRUM: " << data << endl;
	ifstream ifile(data);
	string str;
	vector<double> x, y;
	while(getline(ifile, str)){
		if(str[0] == '#') 
			continue;

		istringstream tmp(str);
		
		double _x, _y;
		tmp >> _x >> _y;
	//	cout << _x << " " << _y << endl;			
		x.push_back(_x);
		y.push_back(_y);
	}
	f = new ROOT::Math::Interpolator(x, y, ROOT::Math::Interpolation::kLINEAR);
}

double chi2_stats(
	double _sin2_2theta13,
	double _sin_theta12,
	double _delta_m2_21,
	double _delta_m2_ee,
	double *_reactor_shape,
	double _weighted_reactor_power[n_det][n_core],
	double *_cross_section,
	size_t energy_bins,
	double *_det_eff_uncorr,
	double _det_eff_corr,
	double *_bkg_uncorr,
	double _norm,
	double obs_rate[n_det][2],
	bool debug,
	bool print
	){

	
	double _chi2 = 0;
	for(size_t det = 0; det < n_det; ++det){
		double expected = 0;
		double expected2 = 0;
		for(size_t core = 0; core < n_core; ++core){
			double _l = dist_map[det][core];
			double _density = 1 / ( _l * _l );
			double _sum = 0;
			double _sum2 = 0;
			for(size_t bin = 0; bin < energy_bins; ++bin){
				double _e = ( 0.5 + bin ) * energy_step + energy_min;
				double _prob = survival_prob(_sin_theta12, _sin2_2theta13, _delta_m2_21, _delta_m2_ee, _l, _e);
				double _prob2 = survival_prob(_sin_theta12, 0, _delta_m2_21, _delta_m2_ee, _l, _e);
				double _n_bin = _prob * _reactor_shape[bin] * _cross_section[bin];
				double _n_bin2 = _prob2 * _reactor_shape[bin] * _cross_section[bin];
				if(debug)
					cout << "bin:" << bin << " e:" << _e << " prob:" << _prob << " r_shape:" << _reactor_shape[bin] << " c_sect:" << _cross_section[bin] << " n_bin:" << _n_bin << endl;
				_sum += _n_bin;
				expected += _prob * _density * _weighted_reactor_power[det][core] * _reactor_shape[bin] * _norm * _cross_section[bin] * energy_step;
				_sum2 += _n_bin2;
				expected2 += _prob2 * _density * _weighted_reactor_power[det][core] * _reactor_shape[bin] * _norm * _cross_section[bin] * energy_step;
			}
			_sum *= (energy_step * _density * _norm);
			_sum2 *= (energy_step * _density * _norm);

			if(debug)
				cout << "sum:" << _sum << endl;

			expected += _sum; 
			expected2 += _sum2; 
		}

		expected *= ((det_eff[0] + _det_eff_corr + _det_eff_uncorr[det]))*eps_mult[det];
		expected2 *= ((det_eff[0] + _det_eff_corr + _det_eff_uncorr[det]))*eps_mult[det];

		if(debug)
			cout << "exp:"<< expected << " obs:" << obs_rate[det][0] << " pull:" << pull(expected + _bkg_uncorr[det], obs_rate[det][0], obs_rate[det][1]) << endl;

		if(0){
			printf("det: %d exp: %10.5f bkg: %10.5f obs: %10.5f/%10.5f chi2: %10.5f\n", det, expected, _bkg_uncorr[det], obs_rate[det][0], obs_rate[det][1], pull(expected + _bkg_uncorr[det], obs_rate[det][0], obs_rate[det][1]));
			fflush(stdout);
		}
        if(0){
            printf("det: %d discrepency: %10.5f; data: %10.5f,expect %10.5f, expect2: %10.5f\n",det,(obs_rate[det][0]-_bkg_uncorr[det])/expected2,obs_rate[det][0]-_bkg_uncorr[det],expected,expected2);
            fflush(stdout);
        }
		_chi2 += pull(expected + _bkg_uncorr[det], obs_rate[det][0], obs_rate[det][1]);

	}	
	return _chi2;	
}

double chi2_syst(
	double _sin2_theta12,
	double _delta_m2_21,
	double _delta_m2_32,
	double *_det_eff_uncorr,
	double _det_eff_corr,
	double *_bkg_uncorr,
	bool print
	){
	
	double _chi2 = 0;
	for(size_t det = 0; det < 8; ++det){
		_chi2 += pull(_bkg_uncorr[det], bkg_rate[det][0], bkg_rate[det][1]);
		_chi2 += pull(_det_eff_uncorr[det], 0, det_eff[2]);
		if(print)
			printf("det: %d bkg_uncorr: %10.5f det_eff_uncorr: %10.5f\n", det, pull(_bkg_uncorr[det], bkg_rate[det][0], bkg_rate[det][1]), pull(_det_eff_uncorr[det], 0, det_eff[2]));
	}
	_chi2 += pull(_det_eff_corr, 0, det_eff[1]);
	_chi2 += pull(_sin2_theta12, sin2_theta12[0], sin2_theta12[1]);
	_chi2 += pull(_delta_m2_21, delta_m2_21[0], delta_m2_21[1]);
	_chi2 += pull(_delta_m2_32, delta_m2_32[0], delta_m2_32[1]);
	if(print){
		printf("det_eff_corr: %10.5f\n", pull(_det_eff_corr, 0, det_eff[1]));
		printf("sin2_theta12: %10.5f\n", pull(_sin2_theta12, sin2_theta12[0], sin2_theta12[1]));
		printf("delta_m2_21: %10.5f\n", pull(_delta_m2_21, delta_m2_21[0], delta_m2_21[1]));
		printf("delta_m2_32: %10.5f\n", pull(_delta_m2_32, delta_m2_32[0], delta_m2_32[1]));
		fflush(stdout);
	}
	return _chi2;	
}


void parseParameters(
	double *par,
	double &_sin2_2theta13, 
	double &_sin2_theta12,
	double &_sin_theta12, 
	double &_delta_m2_21,
	double &_delta_m2_32, 
	double &_delta_m2_ee,
	double _reactor_shape[energy_bins],
	double _det_eff_uncorr[8],
	double &_det_eff_corr,
	double _bkg_uncorr[8],
	double &_norm){

	int par_idx = 0;
	//osc. pars
	_sin2_2theta13 = par[par_idx++];
	_sin2_theta12 = par[par_idx++];
	_sin_theta12 = sqrt(_sin2_theta12);
	_delta_m2_21 = par[par_idx++];
	_delta_m2_32 = par[par_idx++];
	_delta_m2_ee = normal_order? _delta_m2_32 + 5.2e-5: _delta_m2_32 - 5.2e-5;

	//reactor pars
	double _r_frac[reactor_nuclears];
	double _r_frac_sum = 1;
	_r_frac[reactor_nuclears-1] = 1;
	for(size_t n = 0; n < reactor_nuclears-1; ++n){
		double _f = par[par_idx++];
		_r_frac[n] = _f;
		_r_frac_sum += _f;
	}
	for(size_t n = 0; n < reactor_nuclears; ++n){
		_r_frac[n] /= _r_frac_sum;
		if(debug)
			cout << "_r_frac" << n << " " << _r_frac[n] << endl;
	}

	for(size_t bin = 0; bin < energy_bins; ++bin){
		double _e = ( 0.5 + bin ) * energy_step + energy_min;
		double _sum = 0;
		for(size_t n = 0; n < reactor_nuclears; ++n)
			_sum += _r_frac[n] * r_spec[n]->Eval(_e);
		_reactor_shape[bin] = _sum;
	}


	//detector pars
	for(int i=0;i<8;++i)
		_det_eff_uncorr[i] = par[par_idx++];
	_det_eff_corr = par[par_idx++];


	//bkg pars
	for(int i=0;i<8;++i)
		_bkg_uncorr[i] = par[par_idx++];

	_norm = par[par_idx++];

}

void chi2_wrap(int &npar, double *g, double &result, double *par, int flag){

	result = chi2(par, false);
			
	printf("chi2: %20.10f\r", result);
	fflush(stdout);
}

double chi2(double *par, bool print){
	double _sin2_2theta13;
	double _sin2_theta12;
	double _sin_theta12;
	double _delta_m2_21;
	double _delta_m2_32;
	double _delta_m2_ee;
	double _reactor_shape[energy_bins];
	double _det_eff_uncorr[8];
	double _det_eff_corr;
	double _bkg_uncorr[8];
	double _norm;
	
	parseParameters(
	 par,
	 _sin2_2theta13, 
	 _sin2_theta12,
	 _sin_theta12, 
	 _delta_m2_21,
	 _delta_m2_32, 
	 _delta_m2_ee,
	 _reactor_shape,
	 _det_eff_uncorr,
	 _det_eff_corr,
	 _bkg_uncorr,
	 _norm);


	double _chi2_stats = chi2_stats(
		_sin2_2theta13, 
		_sin_theta12, 
		_delta_m2_21,
		_delta_m2_ee, 
		_reactor_shape, 
		_weighted_reactor_power, 
		_cross_section, 
		energy_bins, 
		_det_eff_uncorr, 
		_det_eff_corr,
		_bkg_uncorr, 
		_norm, 
		obs_rate, 
		debug, 
		print);

	double _chi2_syst = chi2_syst(
		_sin2_theta12, 
		_delta_m2_21, 
		_delta_m2_32,
		_det_eff_uncorr, 
		_det_eff_corr, 
		_bkg_uncorr, 
		print);

	if(print){
		printf("chi2_stats: %10.5f chi2_syst: %10.5f chi2_sum: %10.5f\n", _chi2_stats, _chi2_syst, _chi2_stats + _chi2_syst);
		fflush(stdout);
	}
	
	return _chi2_stats + _chi2_syst;
}


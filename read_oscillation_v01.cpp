#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<vector>
#include<map>
#include<set>

#include "WCPLEEANA/TOsc.h"

#include "WCPLEEANA/Configure_Osc.h"

#include "TApplication.h"

//#include <chrono> // timer
//auto time_start = chrono::high_resolution_clock::now();
//auto time_stop = chrono::high_resolution_clock::now();
//auto time_duration = chrono::duration_cast<chrono::seconds>(time_stop - time_start);
//cout<<endl<<" ---> check time duration "<<time_duration.count()<<endl<<endl;
//milliseconds, minutes

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// MAIN //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
	TString roostr = "";

	cout<<endl<<" ---> A story ..."<<endl<<endl;

	int ifile = 1;
	double scaleF_POT_BNB  = 1;
	double scaleF_POT_NuMI = 1;
	int display = 0;

	int it14 = 0;
	int idm2 = 0;
	int it24 = 0;
  
	for(int i=1; i<argc; i++) {
		if( strcmp(argv[i],"-f")==0 ) {
			stringstream convert( argv[i+1] );
			if(  !( convert>>ifile ) ) { cerr<<" ---> Error ifile !"<<endl; exit(1); }
		}
		if( strcmp(argv[i],"-pbnb")==0 ) {
			stringstream convert( argv[i+1] );
			if(  !( convert>>scaleF_POT_BNB ) ) { cerr<<" ---> Error scaleF_POT_BNB !"<<endl; exit(1); }
		}
		if( strcmp(argv[i],"-pnumi")==0 ) {
			stringstream convert( argv[i+1] );
			if(  !( convert>>scaleF_POT_NuMI ) ) { cerr<<" ---> Error scaleF_POT_NuMI !"<<endl; exit(1); }
		}    
		if( strcmp(argv[i],"-d")==0 ) {
			stringstream convert( argv[i+1] );
			if(  !( convert>>display ) ) { cerr<<" ---> Error display !"<<endl; exit(1); }
		}
		if( strcmp(argv[i],"-it14")==0 ) {
			stringstream convert( argv[i+1] );
			if(  !( convert>>it14 ) ) { cerr<<" ---> Error it14 !"<<endl; exit(1); }
		}
		if( strcmp(argv[i],"-idm2")==0 ) {
			stringstream convert( argv[i+1] );
			if(  !( convert>>idm2 ) ) { cerr<<" ---> Error idm2 !"<<endl; exit(1); }
		}
		if( strcmp(argv[i],"-it24")==0 ) {
			stringstream convert( argv[i+1] );
			if(  !( convert>>it24 ) ) { cerr<<" ---> Error it24 !"<<endl; exit(1); }
		}    
	}

	///////////////////////////////////////////////////////////
  
	if( !display ) {
		gROOT->SetBatch( 1 );
	}
  
	TApplication theApp("theApp",&argc,argv);
  
	/////////////////////////////////////////////////////////// Draw style

	gStyle->SetOptStat(0);
	//gStyle->SetPalette(kBird);

	double snWidth = 2;

	// use medium bold lines and thick markers
	gStyle->SetLineWidth(snWidth);
	gStyle->SetFrameLineWidth(snWidth);
	gStyle->SetHistLineWidth(snWidth);
	gStyle->SetFuncWidth(snWidth);
	gStyle->SetGridWidth(snWidth);
	gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
	gStyle->SetMarkerStyle(20);
	gStyle->SetMarkerSize(1.2);
	gStyle->SetEndErrorSize(4);
	gStyle->SetEndErrorSize(0);

	///////////////////////////////////////////////////////////

	TOsc *osc_test = new TOsc();

	///////

	osc_test->scaleF_POT_BNB  = scaleF_POT_BNB;
	osc_test->scaleF_POT_NuMI = scaleF_POT_NuMI;
     
	///////

	osc_test->flag_syst_dirt   = Configure_Osc::flag_syst_dirt;
	osc_test->flag_syst_mcstat = Configure_Osc::flag_syst_mcstat;
	osc_test->flag_syst_flux   = Configure_Osc::flag_syst_flux;
	osc_test->flag_syst_geant  = Configure_Osc::flag_syst_geant;
	osc_test->flag_syst_Xs     = Configure_Osc::flag_syst_Xs;
	osc_test->flag_syst_det    = Configure_Osc::flag_syst_det;
  
	///////

	osc_test->flag_NuMI_nueCC_from_intnue       = Configure_Osc::flag_NuMI_nueCC_from_intnue;
	osc_test->flag_NuMI_nueCC_from_overlaynumu  = Configure_Osc::flag_NuMI_nueCC_from_overlaynumu;
	osc_test->flag_NuMI_nueCC_from_appnue       = Configure_Osc::flag_NuMI_nueCC_from_appnue;
	osc_test->flag_NuMI_nueCC_from_appnumu      = Configure_Osc::flag_NuMI_nueCC_from_appnumu;
	osc_test->flag_NuMI_nueCC_from_overlaynueNC = Configure_Osc::flag_NuMI_nueCC_from_overlaynueNC;
	osc_test->flag_NuMI_nueCC_from_overlaynumuNC= Configure_Osc::flag_NuMI_nueCC_from_overlaynumuNC;
  
	osc_test->flag_NuMI_numuCC_from_overlaynumu   = Configure_Osc::flag_NuMI_numuCC_from_overlaynumu;
	osc_test->flag_NuMI_numuCC_from_overlaynue    = Configure_Osc::flag_NuMI_numuCC_from_overlaynue;
	osc_test->flag_NuMI_numuCC_from_appnue        = Configure_Osc::flag_NuMI_numuCC_from_appnue;
	osc_test->flag_NuMI_numuCC_from_appnumu       = Configure_Osc::flag_NuMI_numuCC_from_appnumu;
	osc_test->flag_NuMI_numuCC_from_overlaynumuNC = Configure_Osc::flag_NuMI_numuCC_from_overlaynumuNC;
	osc_test->flag_NuMI_numuCC_from_overlaynueNC  = Configure_Osc::flag_NuMI_numuCC_from_overlaynueNC;  
  
	osc_test->flag_NuMI_CCpi0_from_overlaynumu  = Configure_Osc::flag_NuMI_CCpi0_from_overlaynumu;
	osc_test->flag_NuMI_CCpi0_from_appnue       = Configure_Osc::flag_NuMI_CCpi0_from_appnue;
	osc_test->flag_NuMI_CCpi0_from_overlaynumuNC= Configure_Osc::flag_NuMI_CCpi0_from_overlaynumuNC;
	osc_test->flag_NuMI_CCpi0_from_overlaynueNC = Configure_Osc::flag_NuMI_CCpi0_from_overlaynueNC;
  
	osc_test->flag_NuMI_NCpi0_from_overlaynumu  = Configure_Osc::flag_NuMI_NCpi0_from_overlaynumu;
	osc_test->flag_NuMI_NCpi0_from_appnue       = Configure_Osc::flag_NuMI_NCpi0_from_appnue;
	osc_test->flag_NuMI_NCpi0_from_overlaynumuNC= Configure_Osc::flag_NuMI_NCpi0_from_overlaynumuNC;
	osc_test->flag_NuMI_NCpi0_from_overlaynueNC = Configure_Osc::flag_NuMI_NCpi0_from_overlaynueNC;


	///////
  
	osc_test->flag_BNB_nueCC_from_intnue       = Configure_Osc::flag_BNB_nueCC_from_intnue;
	osc_test->flag_BNB_nueCC_from_overlaynumu  = Configure_Osc::flag_BNB_nueCC_from_overlaynumu;
	osc_test->flag_BNB_nueCC_from_appnue       = Configure_Osc::flag_BNB_nueCC_from_appnue;
	osc_test->flag_BNB_nueCC_from_appnumu      = Configure_Osc::flag_BNB_nueCC_from_appnumu;
	osc_test->flag_BNB_nueCC_from_overlaynueNC = Configure_Osc::flag_BNB_nueCC_from_overlaynueNC;
	osc_test->flag_BNB_nueCC_from_overlaynumuNC= Configure_Osc::flag_BNB_nueCC_from_overlaynumuNC;
  
	osc_test->flag_BNB_numuCC_from_overlaynumu   = Configure_Osc::flag_BNB_numuCC_from_overlaynumu;
	osc_test->flag_BNB_numuCC_from_overlaynue    = Configure_Osc::flag_BNB_numuCC_from_overlaynue;
	osc_test->flag_BNB_numuCC_from_appnue        = Configure_Osc::flag_BNB_numuCC_from_appnue;
	osc_test->flag_BNB_numuCC_from_appnumu       = Configure_Osc::flag_BNB_numuCC_from_appnumu;
	osc_test->flag_BNB_numuCC_from_overlaynumuNC = Configure_Osc::flag_BNB_numuCC_from_overlaynumuNC;
	osc_test->flag_BNB_numuCC_from_overlaynueNC  = Configure_Osc::flag_BNB_numuCC_from_overlaynueNC;  
  
	osc_test->flag_BNB_CCpi0_from_overlaynumu  = Configure_Osc::flag_BNB_CCpi0_from_overlaynumu;
	osc_test->flag_BNB_CCpi0_from_appnue       = Configure_Osc::flag_BNB_CCpi0_from_appnue;
	osc_test->flag_BNB_CCpi0_from_overlaynumuNC= Configure_Osc::flag_BNB_CCpi0_from_overlaynumuNC;
	osc_test->flag_BNB_CCpi0_from_overlaynueNC = Configure_Osc::flag_BNB_CCpi0_from_overlaynueNC;
  
	osc_test->flag_BNB_NCpi0_from_overlaynumu  = Configure_Osc::flag_BNB_NCpi0_from_overlaynumu;
	osc_test->flag_BNB_NCpi0_from_appnue       = Configure_Osc::flag_BNB_NCpi0_from_appnue;
	osc_test->flag_BNB_NCpi0_from_overlaynumuNC= Configure_Osc::flag_BNB_NCpi0_from_overlaynumuNC;
	osc_test->flag_BNB_NCpi0_from_overlaynueNC = Configure_Osc::flag_BNB_NCpi0_from_overlaynueNC;
  
	/////// set only one time
  
	osc_test->Set_default_cv_cov(Configure_Osc::default_cv_file,
				     Configure_Osc::default_dirtadd_file,
				     Configure_Osc::default_mcstat_file,
				     Configure_Osc::default_fluxXs_dir,
				     Configure_Osc::default_detector_dir);
  
	osc_test->Set_oscillation_base(Configure_Osc::default_eventlist_dir);
  
	/////// Set_oscillation_pars(double val_dm2_41, double val_sin2_2theta_14, double val_sin2_theta_24, double val_sin2_theta_34)
  
	double val_dm2_41         = 7.3;
	double val_sin2_2theta_14 = 0.36;
	double val_sin2_theta_24  = 0;
	double val_sin2_theta_34  = 0;
  
	/////////////////////////////////////////////////////////////////////////////////

	/// standard order
	//  val_dm2_41         = 0;
	//  val_sin2_2theta_14 = 0.1;
	//  val_sin2_theta_24  = 0.1;
  
	//  osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);  
	//  osc_test->Apply_oscillation();  
	//  osc_test->Set_apply_POT();// meas, CV, COV: all ready
  
	// osc_test->Set_meas2fitdata();
  
	// osc_test->Set_asimov2fitdata();

	// osc_test->FCN(const double *par);
 
	// osc_test->Minimization_OscPars_FullCov(dm2, t14, t24, t34, "str_flag_fixpar");

	// osc_test->Set_toy_variations(int num_toys);
	// osc_test->Set_toy2fitdata(int itoy)

	//



	//
	// The statistic test is implemented by the comination of the above functions
	//
	//
	if (0) {
		// Jesse Mendez
		// For Real this time
		// lets go
		// Constants

		// Profiling constants to adjust
		// 3v and small mix angles for profiling
		// (0.0105925,0.000107978,0.0235,0);
		// dm2 t14 t24 t34
		double pars_3v_small[4] = {0, 0.1, 0.1, 0};
		double step_size = 0.0005; // profiling step size
		int steps = 2000; // profiling steps
		const int ntoys = 500; 
	  
		// Profiling data structures
		double profiled_sin2_theta24 = 0;
		double test_sin2_theta24 = 0;
		double chi2_min = pow(10,6);	// set very large to start

		// 4v Hypothesis Data structures
		// Static size to allow writing to root file
		double truth_grid_meas_4v_asimov_chi2[ntoys];
		double truth_3v_meas_4v_asimov_chi2[ntoys];

		// 3v Hypothesis Data Structures
		double truth_grid_meas_3v_asimov_chi2[ntoys];
		double truth_3v_meas_3v_asimov_chi2[ntoys];

		// Delta chisquare Data Structures
		double delta_chi2_grid_toy[ntoys];
		double delta_chi2_3v_toy[ntoys];
	  	  	  
		// Convert command line argument bin numbers to actual values
		// log scale dm^2 80 bins
		TH1D *hist_dm2 = new TH1D("h1_dm2","h1_dm2",80,-2,2);
		// log schale sin squared 3 thetamue 60 bins
		TH1D *hist_sin2_2thetamue = new TH1D("h1_sin2_2thetamue","h1_sin2_2thetamue",60,-4,0);
		double sin2_2theta_mue_grid = hist_sin2_2thetamue->GetBinCenter(it14);
		double dm2_41_grid = hist_dm2->GetBinCenter(idm2);
		sin2_2theta_mue_grid = pow(10.0, sin2_2theta_mue_grid);
		dm2_41_grid = pow(10.0, dm2_41_grid);
	  
		// Profile Section
		// Prepare Asimov data for profiling calculation
		osc_test->Set_oscillation_pars(pars_3v_small[0], pars_3v_small[1], pars_3v_small[2], pars_3v_small[3]);
		osc_test->Apply_oscillation();
		osc_test->Set_apply_POT();
		osc_test->Set_asimov2fitdata();

		// Grid scan
		for(int i = 0; i < steps; i++) {
			test_sin2_theta24 = i * step_size;
			if( sin2_2theta_mue_grid > test_sin2_theta24 ) {
				continue; // No idea why I need this...
			}	
			double pars_4nu[4] = {dm2_41_grid, sin2_2theta_mue_grid, test_sin2_theta24, 0};// dm2, t14, t24, t34
			double chi2 = osc_test->FCN( pars_4nu );
			if(chi2 < chi2_min) {
				chi2_min = chi2;
				profiled_sin2_theta24 = test_sin2_theta24;
			}
		}
		cout << "Profiled sin^2 theta24: " << profiled_sin2_theta24;
		cout << "Chi Square: " << chi2_min;

		// 4v Hypothesis
		double pars_4v_grid[4] ={dm2_41_grid, sin2_2theta_mue_grid, profiled_sin2_theta24, 0};
	  

		// Prepare Asimov data for grid psuedo experiments
		osc_test->Set_oscillation_pars(pars_4v_grid[0], pars_4v_grid[1], pars_4v_grid[2], pars_4v_grid[3]);
		osc_test->Apply_oscillation();
		osc_test->Set_apply_POT();
		osc_test->Set_asimov2fitdata();

		// Generate Asimov Data
		osc_test->Set_toy_variations(ntoys);

		// Calculate chi2
		for(int i=0; i<ntoys; i++) {
			osc_test->Set_toy2fitdata(i+1);
			truth_grid_meas_4v_asimov_chi2[i] = osc_test->FCN(pars_4v_grid);
			truth_3v_meas_4v_asimov_chi2[i] = osc_test->FCN(pars_3v_small);
		}

		// 3v Hypothesis
		// Prepare Asimov data for grid psuedo experiments
		osc_test->Set_oscillation_pars(pars_3v_small[0], pars_3v_small[1], pars_3v_small[2], pars_3v_small[3]);
		osc_test->Apply_oscillation();
		osc_test->Set_apply_POT();
		osc_test->Set_asimov2fitdata();

		// Generate Asimov Data
		osc_test->Set_toy_variations(ntoys);

		// Calculate chi2
		for(int i=0; i<ntoys; i++) {
			osc_test->Set_toy2fitdata(i+1);
			truth_grid_meas_3v_asimov_chi2[i] = osc_test->FCN(pars_4v_grid);
			truth_3v_meas_3v_asimov_chi2[i] = osc_test->FCN(pars_3v_small);
		}

		// Calculate the delta chi squares
		// Collect histograms for Sanity Check
		// 4v Hypothesis Histograms
		TH1F* hist_grid_4v = new TH1F("hist_grid_4v", "Chi2 Distributions Truth grid Meas 4v Asimov", 100, 0, 1 );
		hist_grid_4v->SetCanExtend(TH1::kXaxis); // Change the binning on the fly
		TH1F* hist_3v_4v = new TH1F("hist_3v_4v", "Chi2 Distributions Truth 3v Meas 3v Asimov", 100, 0, 1 );
		hist_3v_4v->SetCanExtend(TH1::kXaxis); // Change the binning on the fly

		// 3v Hypothesis Histograms
		TH1F* hist_grid_3v = new TH1F("hist_grid_3v", "Chi2 Distributions Truth grid Meas 3v Asimov", 100, 0, 1 );
		hist_grid_3v->SetCanExtend(TH1::kXaxis); // Change the binning on the fly
		TH1F* hist_3v_3v = new TH1F("hist_3v_3v", "Chi2 Distributions Truth 3v Meas 3v Asimov", 100, 0, 1 );
		hist_3v_3v->SetCanExtend(TH1::kXaxis); // Change the binning on the fly

		for (int i = 0; i < ntoys; i++) {
			delta_chi2_grid_toy[i] = truth_grid_meas_4v_asimov_chi2[i] - truth_3v_meas_4v_asimov_chi2[i];
	    
			delta_chi2_3v_toy[i] = truth_grid_meas_3v_asimov_chi2[i] - truth_3v_meas_3v_asimov_chi2[i];
			// hist_grid_4v->Fill(truth_grid_meas_4v_asimov_chi2[i]);
			// hist_3v_4v->Fill(truth_3v_meas_4v_asimov_chi2[i]);
			// hist_grid_3v->Fill(truth_grid_meas_3v_asimov_chi2[i]);
			// hist_3v_3v->Fill(truth_3v_meas_3v_asimov_chi2[i]);
	    
			// cout << "Delta chi2 grid toy: " << delta_chi2_grid_toy[i] <<" " <<"Delta chi2 3v toy "<<  delta_chi2_3v_toy[i]<< endl;
	    
		}

		// Save all this to a root file for plotting and stuff in case I am wrong
		roostr = TString::Format("./mendez_osc_analysis_full_data_CLs_dm2_%02d_tt_%02d.root", idm2, it14);
		TFile *rootfile = new TFile(roostr, "recreate");
		TTree *ttree = new TTree("ttree", "ttree");
		ttree->Branch( "idm2", &idm2, "idm2/I" );
		ttree->Branch( "it14", &it14, "it14/I" );
		ttree->Branch( "truth_grid_meas_4v_asimov_chi2",  &truth_grid_meas_4v_asimov_chi2, Form("truth_grid_meas_4v_asimov_chi2[%d]/D", ntoys));
		ttree->Branch( "truth_3v_meas_4v_asimov_chi2",  &truth_3v_meas_4v_asimov_chi2, Form("truth_3v_meas_4v_asimov_chi2[%d]/D", ntoys));
		ttree->Branch( "truth_grid_meas_3v_asimov_chi2",  &truth_grid_meas_3v_asimov_chi2, Form("truth_grid_meas_3v_asimov_chi2[%d]/D", ntoys));
		ttree->Branch( "truth_3v_meas_3v_asimov_chi2",  &truth_3v_meas_3v_asimov_chi2, Form("truth_3v_meas_3v_asimov_chi2[%d]/D", ntoys));
		ttree->Branch( "delta_chi2_grid_toy",  &delta_chi2_grid_toy, Form("delta_chi2_grid_toy[%d]/D", ntoys));
		ttree->Branch( "delta_chi2_3v_toy",  &delta_chi2_3v_toy, Form("delta_chi2_3v_toy[%d]/D", ntoys));

		ttree->Fill();
		rootfile->Write();
		rootfile->Close();
		delete rootfile;
	  
	  
	  
	  
	  

		// // Hanyu Sanity Check
		// hist_grid_4v->SetLineColor(kRed);
		// hist_3v_4v->SetLineColor(kBlue);
		// hist_grid_3v->SetLineColor(kRed);
		// hist_3v_3v->SetLineColor(kBlue);
	  

	  

       
	} // end Jesse Osc Script

	if (1) { // Calculate "data" delta chi square
		
	}
  

// if (0){ //Azimov data set for CL method
	
// 	TH1D* h1_dm2 = new TH1D("h1_dm2","h1_dm2",80,-2,2);//log scale for dm^2
// 	TH1D* h1_sin2_theta14 = new TH1D("h1_sin2_theta14","h1_sin2_theta14",60,-4,1); // log scale for sin^2thet14
// 	TH1D* h1_sin2_theta24 = new TH1D("h1_sin2_btheta24","h1_sin2_theta24",60,-4,1);// //	double val_idm2 = h1_dm2->GetBinCenter(idm2);
// //	double val_sin2_theta14 = h1_sin2_theta14->GetBinCenter(it14);
// //	double val_sin2_theta24 = h1_sin2_theta24->GetBinCenter(it24);
 
// //	for(int t24=1;t24<61;++t24){
// 		for(int dm2=1;dm2<81;++dm2){
// 		double val_idm2 = h1_dm2->GetBinCenter(dm2);
//         	double val_sin2_theta14 = h1_sin2_theta14->GetBinCenter(it14);
//         	double val_sin2_theta24 = h1_sin2_theta24->GetBinCenter(it24);
// 	val_idm2 = pow(10.0,val_idm2);
// 		val_sin2_theta14 = pow(10.0,val_sin2_theta14);
// 		val_sin2_theta24 = pow(10.0,val_sin2_theta24);
// 		osc_test->Set_oscillation_pars(val_idm2, val_sin2_theta14, val_sin2_theta24, 0.0);
// 		osc_test->Apply_oscillation();
// 		osc_test->Set_apply_POT(); //4nu prediction is ready
// 		osc_test->Set_asimov2fitdata(); // 4nu hypothesis is ready
// 		//deltachi^2 = chi^2_4nu - chi^2_3nu
// 		double pars_4nu[4] = {val_idm2,val_sin2_theta14, val_sin2_theta24, 0};// dm2, t14, t24, t34
//     		double chi2_4nu_with4nudata = osc_test->FCN( pars_4nu );
// 		double pars_3nu[4] = {0.0,0.1,0.1,0.0};
// 		double chi2_3nu_with4nudata = osc_test->FCN( pars_3nu );
// 		double deltachi2_with4nudata = chi2_4nu_with4nudata-chi2_3nu_with4nudata; // CL_{s+b}


// 		osc_test->Set_oscillation_pars(0.0, val_sin2_theta14, val_sin2_theta24, 0.0);
//                 osc_test->Apply_oscillation();
//                 osc_test->Set_apply_POT(); //4nu prediction is ready
//                 osc_test->Set_asimov2fitdata(); // 4nu hypothesis is ready
//                 //deltachi^2 = chi^2_4nu - chi^2_3nu
//                 //double pars_4nu[4] = {val_idm2,val_sin2_theta14, val_sin2_theta24, 0};// dm2, t14, t24, t34
//                 double chi2_4nu_with3nudata = osc_test->FCN( pars_4nu );
//                	//double pars_3nu[4] = {0.0,0.1,0.1,0.0};
//                 double chi2_3nu_with3nudata = osc_test->FCN( pars_3nu );
//                 double deltachi2_with3nudata = chi2_4nu_with3nudata-chi2_3nu_with3nudata; // CL_{b}
// 		std::cout<<"dm2="<<dm2<<"; chi2_4nu_with4nudata="<<chi2_4nu_with4nudata<<"; chi2_3nu_with4nudata="<<chi2_3nu_with4nudata<<"; deltachi2_with4nudata="<<deltachi2_with4nudata<<std::endl;

// 		std::cout<<"dm2="<<dm2<<"; chi2_4nu_with3nudata="<<chi2_4nu_with3nudata<<"; chi2_3nu_with3nudata="<<chi2_3nu_with3nudata<<"; deltachi2_with3nudata="<<deltachi2_with3nudata<<std::endl;
// }
// //	}
// }
// // for asimov
// // if (0){ //asimov sample calculation gives delta chi square this is what I am loo// king for to generate the single delta chi-square
// // // This set of instructions performs psuedo experiments

// //     	vector<double> vec_chi2_4v_on_4vPseudo;
// //     	vector<double> vec_chi2_3v_on_4vPseudo;
// //     	vector<double> vec_chi2_4v_on_3vPseudo;
// //     	vector<double> vec_chi2_3v_on_3vPseudo;

// //     //	roostr = TString::Format("/direct/lbne+u/smartynen/job_test/Profiling/ToAll_basic_nospeedup/test_dm2_%02d_tt_%02d.root", idm2, it14);
// // roostr = TString::Format("/home/jmendez/code/TOsc_nueapp_prof_edit_sergey/out/asimov_CLs_dm2_%02d_tt_%02d.root", idm2, it14);    	
// // TFile *subroofile_a = new TFile(roostr, "recreate");
// //     	TTree *tree_a = new TTree("tree_a", "tree_a");

// //     	tree_a->Branch( "idm2",          &idm2,     "idm2/I" );
// //     	tree_a->Branch( "it14",          &it14,     "it14/I" );
// //     	tree_a->Branch( "vec_chi2_4v_on_3vAsimov",  &vec_chi2_4v_on_3vPseudo );
	
// // 	TH1D* h1_dm2 = new TH1D("h1_dm2","h1_dm2",80,-2,2);//log scale for dm^2
// //         TH1D* h1_sin2_theta14 = new TH1D("h1_sin2_theta14","h1_sin2_theta14",60,-4,0); // log scale for sin^2thet14
// //         TH1D* h1_sin2_theta24 = new TH1D("h1_sin2_theta24","h1_sin2_theta24",60,-4,0);
// // 	double val_idm2 = h1_dm2->GetBinCenter(idm2);
// //         double val_sin2_theta14 = h1_sin2_theta14->GetBinCenter(it14);
// // 	double val_sin2_theta24 = h1_sin2_theta24->GetBinCenter(it24);
// //         val_idm2 = pow(10.0,val_idm2);
// //         val_sin2_theta14 = pow(10.0,val_sin2_theta14);

// // 	osc_test->Set_oscillation_pars(0.0, val_sin2_theta14, val_sin2_theta24, 0.0);
// //         osc_test->Apply_oscillation();
// //         osc_test->Set_apply_POT(); //3nu prediction is ready
// //         	osc_test->Set_asimov2fitdata();
// // //		osc_test->Set_meas2fitdata();
// //                 //deltachi^2 = chi^2_4nu - chi^2_3nu
// //                 double pars_4nu[4] = {val_idm2,val_sin2_theta14, val_sin2_theta24, 0};// dm2, t14, t24, t34
// //                 double chi2_4nu_with3nudata = osc_test->FCN( pars_4nu );
// //                 vec_chi2_4v_on_3vPseudo.push_back(chi2_4nu_with3nudata);
	
// // 	subroofile_input->Close();	
// // 	tree_a->Fill();
// // 	subroofile_a->Write();
// // 	subroofile_a->Close();
// // //	delete subroofile;
// // //	delete tree;
	

// // }



// if (0){ //Real data calculation 


//     	vector<double> vec_chi2_4v_on_4vPseudo;
//     	vector<double> vec_chi2_3v_on_4vPseudo;
//     	vector<double> vec_chi2_4v_on_3vPseudo;
//     	vector<double> vec_chi2_3v_on_3vPseudo;

// 	roostr = TString::Format("/gpfs/mnt/gpfs01/lbne/users/smartynen/job_output/final_setup/prof/nueApp/data_as_norm30/grid/data_CLs_dm2_%02d_tt_%02d.root", idm2, it14);    
// //roostr = TString::Format("/gpfs/mnt/gpfs01/lbne/users/smartynen/job_output/final_setup/prof/nueApp/data_numi_as_bnb/grid/data_CLs_dm2_%02d_tt_%02d.root", idm2, it14);
// 	TFile *subroofile_d = new TFile(roostr, "recreate");
//     	TTree *tree_d = new TTree("tree_d", "tree_d");

//     	tree_d->Branch( "idm2",          &idm2,     "idm2/I" );
//     	tree_d->Branch( "it14",          &it14,     "it14/I" );
//     	tree_d->Branch( "vec_chi2_4v_on_3vRealData",  &vec_chi2_4v_on_3vPseudo );
// 	tree_d->Branch( "vec_chi2_3v_on_3vRealData",  &vec_chi2_3v_on_3vPseudo );

// 	TH1D* h1_dm2 = new TH1D("h1_dm2","h1_dm2",80,-2,2);
// 	TH1D* h1_sin2_theta14 = new TH1D("h1_sin2_theta14","h1_sin2_theta14",60,-4,0); 
//         TH1D* h1_sin2_theta24 = new TH1D("h1_sin2_theta24","h1_sin2_theta24",60,-4,0);
// 	double val_idm2 = h1_dm2->GetBinCenter(idm2);
//         double val_sin2_theta14 = h1_sin2_theta14->GetBinCenter(it14);
//         val_idm2 = pow(10.0,val_idm2);
//         val_sin2_theta14 = pow(10.0,val_sin2_theta14);

// 	roostr = TString::Format("/gpfs/mnt/gpfs01/lbne/users/smartynen/job_output/final_setup/prof/nueApp/fit_res_as_norm30/fitres_dat_CLs_dm2_%02d_tt_%02d.root", idm2, it14);
//         TFile *subroofile_input = new TFile(roostr, "read");
//         TTree *input = (TTree*)subroofile_input->Get("tree");
//         double result;
//         input->SetBranchAddress("theta24",  &result);
//         input->GetEntry(0);
//         double val_sin2_theta24 = result;
	
// 	osc_test->Set_oscillation_pars(0.0, val_sin2_theta14, val_sin2_theta24, 0.0);
//         osc_test->Apply_oscillation();
//         osc_test->Set_apply_POT(); 
// //        osc_test->Set_asimov2fitdata();
    
//   	osc_test->Set_meas2fitdata();
// /*for(int i=1;i<h1_dm2->GetXaxis()->GetNbins()+1;++i){
// for(int j=1;j<h1_sin2_theta14->GetXaxis()->GetNbins()+1;++j){
// std::cout<<i<<" "<<j<<std::endl;
// val_idm2 = h1_dm2->GetBinCenter(i);
// val_sin2_theta14 = h1_sin2_theta14->GetBinCenter(j);
//  val_idm2 = pow(10.0,val_idm2);
// val_sin2_theta14 = pow(10.0,val_sin2_theta14);
// */
// std::cout<<"t24="<<val_sin2_theta24<<std::endl;
//           double pars_4nu[4] = {val_idm2,val_sin2_theta14, val_sin2_theta24, 0};
//                 double chi2_4nu_with3nudata = osc_test->FCN( pars_4nu );
//                 vec_chi2_4v_on_3vPseudo.push_back(chi2_4nu_with3nudata);
// //	}
// //a}
//                 double pars_3nu[4] = {0.0,val_sin2_theta14, val_sin2_theta24, 0};
//                 double chi2_3nu_with3nudata = osc_test->FCN( pars_3nu );
//                 vec_chi2_3v_on_3vPseudo.push_back(chi2_3nu_with3nudata);
// 	subroofile_input->Close();	
// 	tree_d->Fill();
//         subroofile_d->Write();
//         subroofile_d->Close();
// }


// if (0){//fit for theta24 Profiling


//         double result;
// 	int conv;
// 	double chi2null;
// //roostr = TString::Format("/direct/lbne+u/smartynen/job_test/Profiling/ToAll_basic_nospeedup/fitres_CLs_dm2_%02d_tt_%02d.root", idm2, it14);
// roostr = TString::Format("/home/jmendez/code/TOsc_nueapp_prof_edit_sergey/out/fitres_dat_CLs_dm2_%02d_tt_%02d.root", idm2, it14); //sets output file
// //roostr = TString::Format("/gpfs/mnt/gpfs01/lbne/users/smartynen/job_output/final_setup/prof/nueApp/fit_res_as/test_%02d_tt_%02d.root", idm2, it14);        

// TFile *subroofile = new TFile(roostr, "recreate");
//         TTree *tree = new TTree("tree", "tree");

//         tree->Branch( "idm2",          &idm2,     "idm2/I" );
//         tree->Branch( "it14",          &it14,     "it14/I" );
// 	tree->Branch( "conv",&conv,"conv/I");
//         tree->Branch( "theta24",  &result, "theta24/D" );
// 	tree->Branch( "chi2null",&chi2null,"chi2null/D");
// 	// input -idm2 * -it14 *
// 	TH1D* h1_dm2 = new TH1D("h1_dm2","h1_dm2",80,-2,2);//log scale for dm^2
//         TH1D* h1_sin2_2thetamue = new TH1D("h1_sin2_2thetamue","h1_sin2_2thetamue",60,-4,0); // log scale for sin^22thetmue
// 	double val_idm2 = h1_dm2->GetBinCenter(idm2);
//         double val_sin2_2thetamue = h1_sin2_2thetamue->GetBinCenter(it14);
//         val_idm2 = pow(10.0,val_idm2);
//         val_sin2_2thetamue = pow(10.0,val_sin2_2thetamue);
// //for sensetivity curve we use 3nu asimov data
//   	osc_test->Set_oscillation_pars(0,0.1,0.1,0);//(0.0105925,0.000107978,0.0235,0);// dm2, t14, t24, t34
//     	osc_test->Apply_oscillation();// Oscillation formula: Prob_oscillaion(double Etrue, double baseline, int strflag_osc);
//     	osc_test->Set_apply_POT();// meas, CV, COV ---> all ready
// 	osc_test->Set_asimov2fitdata();
// 	// osc_test->Set_meas2fitdata();// set the measured data xs
	
// 	//minimization part 
// 	//cource grid scan [to find initial pararmeters] for parameter sin^2theta24 [0,1] for BNB 0.005 is the step, called degenerecy point
// 	//define to regions [0,0.01] near the degenerecy point with ~20 grids and the other one [0.01,1] with grid ~200 far from the point 
// 	double chi2_min = pow(10,6);	
// 	double min_theta = 0;
//         std::vector<double> chi2_t24(5000,-10);
//         tree->Branch( "chi2_t24",  &chi2_t24);
// 	for(int i=1;i<=2000;++i){
// 	double t24=0.0005*i;
// 	if( val_sin2_2thetamue > t24 ) {
//               continue;
//             }	
// 		//std::cout<<i<<" "<<val_idm2<<" "<<val_sin2_2thetamue<<" "<<t24<<std::endl;
//                 double pars_4nu[4] = {val_idm2,val_sin2_2thetamue, t24, 0};// dm2, t14, t24, t34
//                 double chi2 = osc_test->FCN( pars_4nu );
// chi2_t24[i-1]=chi2;
// 		//std::cout<<chi2<<std::endl;
// 	if(chi2<chi2_min){
//                 chi2_min=chi2;
//                 min_theta=t24;
//         }
// std::cout<<i<<" "<<chi2_min<<" "<<min_theta<<std::endl;
// 	}

//  /*  double initial_dm2 = val_idm2;
//     double initial_t14 = val_sin2_2thetamue;
//     double initial_t24 = min_theta;
//     osc_test->Minimization_OscPars_FullCov(initial_dm2, initial_t14, initial_t24, 0, "t14_dm2");
// std::cout<<osc_test->minimization_status<<std::endl;//if it is not zero -> minimization failed sometimes it happens need to account for it 
// conv=osc_test->minimization_status;
// if(conv==0){
// result = osc_test->minimization_sin2_theta_24_val;
// }else{
// result=initial_t24;
// }
// */
// result=min_theta;
// chi2null=chi2_min;
// //input for CL pseudo experiments 4nu hypothesis would be val_idm2, val_sin2_2thetamue, min_sin2_theta_24_val, 0
//         tree->Fill();
//         subroofile->Write();
//         subroofile->Close();
// }





// if (0){ //Pseudo experiments for profiling
	
//   	vector<double> vec_chi2_4v_on_4vPseudo;
//     	vector<double> vec_chi2_3v_on_4vPseudo;
//     	vector<double> vec_chi2_4v_on_3vPseudo;
//     	vector<double> vec_chi2_3v_on_3vPseudo;
// ///gpfssmnt/gpfs01/lbne/users/smartynen/job_output/final_setup/prof/nueApp/brazil_bnbnumi
// //roostr = TString::Format("/gpfs/mnt/gpfs01/lbne/users/smartynen/job_output/final_setup/prof/nueApp/brazil_bnbnumi/sub_frequentist_CLs_dm2_%02d_tt_%02d.root", idm2, it14);    
// 	// roostr = TString::Format("/gpfs/mnt/gpfs01/lbne/users/smartynen/job_output/final_setup/prof/nueApp/data_as_norm30/pseudo/sub_frequentist_CLs_dm2_%02d_tt_%02d.root", idm2, it14);
// 	roostr = TString::Format("/home/jmendez/code/TOsc_nueapp_prof_edit_sergey/out/sub_frequentist_CLs_dm2_%02d_tt_%02d.root", idm2, it14);
// //roostr = TString::Format("/direct/lbne+u/smartynen/job_test/Profiling/ToAll_basic_nospeedup/test_CLs_dm2_%02d_tt_%02d.root", idm2, it14);  
//   	TFile *subroofile = new TFile(roostr, "recreate");
//     	TTree *tree = new TTree("tree", "tree");

//     	tree->Branch( "idm2",          &idm2,     "idm2/I" );
//     	tree->Branch( "it14",          &it14,     "it14/I" );
//     	tree->Branch( "vec_chi2_4v_on_4vPseudo",  &vec_chi2_4v_on_4vPseudo );
//     	tree->Branch( "vec_chi2_3v_on_4vPseudo",  &vec_chi2_3v_on_4vPseudo );
//     	tree->Branch( "vec_chi2_4v_on_3vPseudo",  &vec_chi2_4v_on_3vPseudo );
//     	tree->Branch( "vec_chi2_3v_on_3vPseudo",  &vec_chi2_3v_on_3vPseudo );
	
// 	TH1D* h1_dm2 = new TH1D("h1_dm2","h1_dm2",80,-2,2);//log scale for dm^2
//         TH1D* h1_sin2_theta14 = new TH1D("h1_sin2_theta14","h1_sin2_theta14",60,-4,0); // log scale for sin^2thet14
//         TH1D* h1_sin2_theta24 = new TH1D("h1_sin2_theta24","h1_sin2_theta24",60,-4,0);
// 	double val_idm2 = h1_dm2->GetBinCenter(idm2);
//         double val_sin2_theta14 = h1_sin2_theta14->GetBinCenter(it14);
// 	roostr = TString::Format("/home/jmendez/code/TOsc_nueapp_prof_edit_sergey/out/fitres_dat_CLs_dm2_%02d_tt_%02d.root", idm2, it14);
// 	TFile *subroofile_input = new TFile(roostr, "read");
// 	TTree *input = (TTree*)subroofile_input->Get("tree");
// 	double result;
// 	input->SetBranchAddress("theta24",  &result);
// 	input->GetEntry(0);
//        	double val_sin2_theta24 = result;
//         val_idm2 = pow(10.0,val_idm2);
//         val_sin2_theta14 = pow(10.0,val_sin2_theta14);
//         //val_sin2_theta24 = pow(10.0,val_sin2_theta24);
// 	osc_test->Set_oscillation_pars(val_idm2, val_sin2_theta14, val_sin2_theta24, 0.0);
//         osc_test->Apply_oscillation();
//         osc_test->Set_apply_POT(); //4nu prediction is ready
// //TMatrixD saved_data = osc_test->matrix_eff_newworld_pred;     
//    int ntoys=1000;
// 	osc_test->Set_toy_variations( ntoys );
	
// for(int toy=1;toy<ntoys+1;++toy){
// //for(int i=0;i<26*7;++i){osc_test->map_matrix_toy_pred[toy](0,i)=saved_data(0,i);}
// //for(int i=26*9;i<26*14;++i){osc_test->map_matrix_toy_pred[toy](0,i)=saved_data(0,i);} 
//                osc_test->Set_toy2fitdata( toy );
//                 double pars_4nu[4] = {val_idm2,val_sin2_theta14, val_sin2_theta24, 0};// dm2, t14, t24, t34
//                 double chi2_4nu_with4nudata = osc_test->FCN( pars_4nu );
//                 double pars_3nu[4] = {0.0,0.1,0.1,0.0};
//                 double chi2_3nu_with4nudata = osc_test->FCN( pars_3nu );
//                 vec_chi2_4v_on_4vPseudo.push_back(chi2_4nu_with4nudata);
//                 vec_chi2_3v_on_4vPseudo.push_back(chi2_3nu_with4nudata);
//         }

//         osc_test->Set_oscillation_pars(0.0, val_sin2_theta14, val_sin2_theta24, 0.0);
//         osc_test->Apply_oscillation();
//         osc_test->Set_apply_POT(); 
// //saved_data = osc_test->matrix_eff_newworld_pred;
//         osc_test->Set_toy_variations(ntoys);
//         for(int toy=1;toy<ntoys+1;++toy){
// //for(int i=0;i<26*7;++i){osc_test->map_matrix_toy_pred[toy](0,i)=saved_data(0,i);}
// //for(int i=26*9;i<26*14;++i){osc_test->map_matrix_toy_pred[toy](0,i)=saved_data(0,i);} 
//                 osc_test->Set_toy2fitdata(toy);
//                 double pars_4nu[4] = {val_idm2,val_sin2_theta14, val_sin2_theta24, 0};// dm2, t14, t24, t34
//                 double chi2_4nu_with3nudata = osc_test->FCN( pars_4nu );
//                 double pars_3nu[4] = {0.0,0.1,0.1,0.0};
//                 double chi2_3nu_with3nudata = osc_test->FCN( pars_3nu );
//                 vec_chi2_4v_on_3vPseudo.push_back(chi2_4nu_with3nudata);
//                 vec_chi2_3v_on_3vPseudo.push_back(chi2_3nu_with3nudata);
//         }

// subroofile_input->Close();
// 	//delete input;
	
// 	tree->Fill();
// 	subroofile->Write();
// 	subroofile->Close();
// }


	cout<<endl;
	cout<<" ------------------------------ check at the final step ------------------------------"<<endl;
	cout<<" ------------------------------ check at the final step ------------------------------"<<endl;

	cout<<endl;
	cout<<TString::Format(" ---> display(-d) %d, ifile(-f) %d, scaleF_POT_BNB(-pbnb) %6.4f, scaleF_POT_NuMI(-pnumi) %6.4f, dm2(-idm2) %d, theta14(-it14) %d, theta24(-it24) %d",
			      display, ifile, osc_test->scaleF_POT_BNB, osc_test->scaleF_POT_NuMI, idm2, it14, it24)<<endl;

	cout<<endl;
	cout<<TString::Format(" ---> Files/Directories: ")<<endl;
	cout<<TString::Format("      ---> default_cv_file       %-50s", Configure_Osc::default_cv_file.Data())<<endl;
	cout<<TString::Format("      ---> default_dirtadd_file  %-50s", Configure_Osc::default_dirtadd_file.Data())<<endl;
	cout<<TString::Format("      ---> default_mcstat_file   %-50s", Configure_Osc::default_mcstat_file.Data())<<endl;
	cout<<TString::Format("      ---> default_fluxXs_dir    %-50s", Configure_Osc::default_fluxXs_dir.Data())<<endl;
	cout<<TString::Format("      ---> default_detector_dir  %-50s", Configure_Osc::default_detector_dir.Data())<<endl;
	cout<<TString::Format("      ---> default_eventlist_dir %-50s", Configure_Osc::default_eventlist_dir.Data())<<endl;
  
	cout<<endl;
	cout<<TString::Format(" ---> flag_syst_dirt    %d", osc_test->flag_syst_dirt)<<endl;
	cout<<TString::Format(" ---> flag_syst_mcstat  %d", osc_test->flag_syst_mcstat)<<endl;
	cout<<TString::Format(" ---> flag_syst_flux    %d", osc_test->flag_syst_flux)<<endl;
	cout<<TString::Format(" ---> flag_syst_geant   %d", osc_test->flag_syst_geant)<<endl;
	cout<<TString::Format(" ---> flag_syst_Xs      %d", osc_test->flag_syst_Xs)<<endl;
	cout<<TString::Format(" ---> flag_syst_det     %d", osc_test->flag_syst_det)<<endl;

	cout<<endl;
	cout<<" ---> Finished sucessfully"<<endl;
  
	cout<<endl;
	if (display) {
		cout << " Enter Ctrl+c to end the program" << endl;
		cout << " Enter Ctrl+c to end the program" << endl;
		cout << endl;
		theApp.Run();
	}

	return 0;
}

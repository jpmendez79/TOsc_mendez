#include <iostream>
#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>

using namespace std;

int main(int argc, char *argv[]) 
{
	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " path_to_your_file.root min max" << std::endl;
        return 1;
    }

	const char *file_path = argv[1];
	TFile *rootfile = TFile::Open(file_path);
	
	TTree *tree = (TTree*)rootfile->Get("ttree");
	TLeaf *leaf = tree->GetLeaf("truth_grid_meas_4v_asimov_chi2");
	int array_size = leaf->GetLen();
	cout << array_size <<endl;

	double *truth_grid_meas_4v_asimov_chi2 = new double[array_size];
	double *truth_3v_meas_4v_asimov_chi2 = new double[array_size];
	double *truth_grid_meas_3v_asimov_chi2 = new double[array_size];
	double *truth_3v_meas_3v_asimov_chi2 = new double[array_size];
	double *delta_chi2_grid_toy = new double[array_size];
	double *delta_chi2_3v_toy = new double[array_size];

	int idm2;
	int it14;
	
	tree->SetBranchAddress("idm2", &idm2);
	tree->SetBranchAddress("it14", &it14);
	
		
	tree->SetBranchAddress("truth_grid_meas_4v_asimov_chi2", truth_grid_meas_4v_asimov_chi2);
	tree->SetBranchAddress("truth_3v_meas_4v_asimov_chi2", truth_3v_meas_4v_asimov_chi2);
	tree->SetBranchAddress("truth_grid_meas_3v_asimov_chi2", truth_grid_meas_3v_asimov_chi2);
	tree->SetBranchAddress("truth_3v_meas_3v_asimov_chi2", truth_3v_meas_3v_asimov_chi2);
	tree->SetBranchAddress("delta_chi2_grid_toy", delta_chi2_grid_toy);
	tree->SetBranchAddress("delta_chi2_3v_toy", delta_chi2_3v_toy);
	
	tree->GetEntry(0);

	TH1F *delta_chi2_grid = new TH1F("hist_delta_grid", "Grid 4v asimov",100, 0, 3000);
	delta_chi2_grid->SetCanExtend(TH1::kAllAxes);
	TH1F *delta_chi2_3v = new TH1F("hist_delta_3v", "3v 3v Asimov",100, 0, 3000);
	delta_chi2_3v->SetCanExtend(TH1::kAllAxes);

	// Sets the initial minimum value of the first value
	
	
	for (int i=0; i<array_size; i++) {
		delta_chi2_grid->Fill(delta_chi2_grid_toy[i]);
		delta_chi2_3v->Fill(delta_chi2_3v_toy[i]);
		// // Test if the initial min is smaller than the current grid toy array
		// plot_min = min(plot_min, delta_chi2_grid_toy[0]);
		// // Test if the new min is smaller thant he current 3v toy
		// plot_min = min(plot_min ,delta_chi2_3v_toy[0]);
		// // Test if the initial min is smaller than the current grid toy array
		// plot_max = max(plot_max, delta_chi2_grid_toy[0]);
		// // Test if the new min is smaller thant he current 3v toy
		// plot_max = max(plot_max, delta_chi2_3v_toy[0]);

				
	}
	
	TCanvas* c1 = new TCanvas("c1", "Canvas", 800, 600);
	double plot_min = min(delta_chi2_grid->GetXaxis()->GetXmin(), delta_chi2_3v->GetXaxis()->GetXmin());
	double plot_max = max(delta_chi2_grid->GetXaxis()->GetXmax(), delta_chi2_3v->GetXaxis()->GetXmax());
	
	if (argc == 4) {
		plot_min = stof(argv[2]);
		plot_max = stof(argv[3]);

	}
	

	delta_chi2_grid->GetXaxis()->SetRangeUser(plot_min, plot_max);
	delta_chi2_3v->GetXaxis()->SetRangeUser(plot_min, plot_max);
	delta_chi2_grid->SetLineColor(kRed);
	delta_chi2_3v->SetLineColor(kBlue);
	delta_chi2_grid->SetTitle("Delta Chi2 Distribution");
	delta_chi2_grid->Draw();
	delta_chi2_3v->Draw("SAME");
	TString hist_png_name = TString::Format("./delta_chi_dm2_%02d_tt_%02d.png", idm2, it14);
	c1->SaveAs(hist_png_name);
	
	

	
	
	
}

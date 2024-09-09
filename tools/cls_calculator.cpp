#include <iostream>
#include <TH1D.h>
#include <cstring>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>

using namespace std;


int main(int argc, char *argv[]) 
{
	// Open file
	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " path_to_your_file.root min max" << std::endl;
		return 1;
	}

	const char *file_path = argv[1];
	TFile *rootfile = TFile::Open(file_path);
	
	TTree *tree = (TTree*)rootfile->Get("ttree");
	TLeaf *leaf = tree->GetLeaf("truth_grid_meas_4v_asimov_chi2");
	int array_size = leaf->GetLen();
	// cout << array_size <<endl;

	double *truth_grid_meas_4v_asimov_chi2 = new double[array_size];
	double *truth_3v_meas_4v_asimov_chi2 = new double[array_size];
	double *truth_grid_meas_3v_asimov_chi2 = new double[array_size];
	double *truth_3v_meas_3v_asimov_chi2 = new double[array_size];
	double *delta_chi2_grid_toy = new double[array_size];
	double *delta_chi2_3v_toy = new double[array_size];

	double delta_chisquare_data;
	int idm2;
	int it14;
	
	tree->SetBranchAddress("idm2", &idm2);
	tree->SetBranchAddress("it14", &it14);
	tree->SetBranchAddress("delta_chisquare_data", &delta_chisquare_data);
	
		
	tree->SetBranchAddress("truth_grid_meas_4v_asimov_chi2", truth_grid_meas_4v_asimov_chi2);
	tree->SetBranchAddress("truth_3v_meas_4v_asimov_chi2", truth_3v_meas_4v_asimov_chi2);
	tree->SetBranchAddress("truth_grid_meas_3v_asimov_chi2", truth_grid_meas_3v_asimov_chi2);
	tree->SetBranchAddress("truth_3v_meas_3v_asimov_chi2", truth_3v_meas_3v_asimov_chi2);
	tree->SetBranchAddress("delta_chi2_grid_toy", delta_chi2_grid_toy);
	tree->SetBranchAddress("delta_chi2_3v_toy", delta_chi2_3v_toy);
	
	tree->GetEntry(0);

	double cl_s = 0;
	double cl_sb = 0;

	for (int i=0; i<array_size; i++) {
		if (delta_chi2_grid_toy[i] > delta_chisquare_data) {
			cl_sb++;
			
		}
		if (delta_chi2_3v_toy[i] > delta_chisquare_data) {
			cl_s++;
			
		}
		
	}

	cl_sb = cl_sb/array_size;
	cl_s = cl_s/array_size;
	
	cout << "CL_s+b: " << cl_sb << endl;
	cout << "CL_s: " << cl_s <<endl;
	


}

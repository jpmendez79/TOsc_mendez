#include <iostream>
#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLine.h>
#include <vector>
#include <cmath>
#include <TGraph.h>

using namespace std;

int main(int argc, char *argv[]) {

  const char *file_path = argv[1];
  TFile *rootfile = TFile::Open(file_path);

  TTree *tree = (TTree *)rootfile->Get("ttree");

  vector<double>* xdata = nullptr;
  vector<double>* ydata = nullptr;

  tree->SetBranchAddress("curve_sin2_2theta_ee", &xdata);
  tree->SetBranchAddress("curve_delta_msquare_41", &ydata);
  tree->GetEntry(0);

  // Create a canvas
    TCanvas *c1 = new TCanvas("c1", "Canvas with Log Scale", 800, 600);
    
    // Set logarithmic scale on X and Y axes
    // c1->SetLogx();
    // c1->SetLogy();
    TGraph *graph = new TGraph(xdata->size());
    
  for (int i=0; i<xdata->size(); i++) {
	  // cout << log10(xdata->at(i))<< "," << log10(ydata->at(i)) << endl;
	  graph->SetPoint(i, xdata->at(i), ydata->at(i));	  
  }
    // TLine *line = new TLine(xdata, ydata);
    // line->Draw();
  graph->Draw();
    graph->SetMarkerStyle(20);
    c1->Update();
    c1->SaveAs("curve.png");
    
    
}

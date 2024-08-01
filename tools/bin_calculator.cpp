#include <iostream>
#include <TH1D.h>
#include <cstring>
#include <cmath>
using namespace std;

int main(int argc, char *argv[]) 
{

	  // log scale dm^2 80 bins
	  TH1D *hist_dm2 = new TH1D("h1_dm2","h1_dm2",80,-2,2);
	  // log schale sin squared 3 thetamue 60 bins
	  TH1D *hist_theta = new TH1D("h1_thet","h1_theta",60,-4,0);

	  if(strcmp(argv[1], "bin2val") == 0) {
		  // Convert bin to value
		  double theta = stof(argv[2]);
		  double dm2 = stof(argv[3]);
		  		  
		  double theta_grid = hist_theta->GetBinCenter(theta);
		  double dm2_grid = hist_dm2->GetBinCenter(dm2);
		  theta_grid = pow(10.0, theta_grid);
		  dm2_grid = pow(10.0, dm2_grid);

		  cout << "theta: " << theta_grid << endl;
		  cout << "dm2: " << dm2_grid << endl;
		  
	  } else if (strcmp(argv[1],"val2bin") == 0){
		  // Convert value to bin
		  double theta = stof(argv[2]);
		  double dm2 = stof(argv[3]);

		  theta = log10(theta);
		  dm2 = log10(dm2);
		  int theta_bin = hist_theta->FindBin(theta);
		  int dm2_bin = hist_dm2->FindBin(dm2);

		  cout << "theta: " << theta_bin << endl;
		  cout << "dm2: " << dm2_bin << endl;

		  
	  }
	
}


#include <string>
#include <sstream>
#include <iostream>

#include "XeNeuSimsAnalysisEnvironment.hh"

#include "TMinuit.h"
#include "TMath.h"

XeNeuSimsAnalysisEnvironment x;

TMinuit * minuit = new TMinuit(10);

// Define the function to minimize
void FCN(int &npar, double *grad, double &fval, double *par, int flag) {
  fval = 0.;
  int ch_num = TMath::Nint(par[2]);
  x.GenerateYieldDist(par[0],par[1],0.18,0.95,0.95,ch_num);
  fval = x.ComputeNegativeLogLikelihood( ch_num );
  std::cout << "par[0]: " << par[0] << " par[1] " << par[1] << 
               " fval: " << fval << std::endl;
}



int main() {

 x.AddDataHistFiles( 102 );

 std::stringstream simsfile;
 simsfile << "/g/g20/lenardo1/Simulations/" <<
             "tunl_sims_conical_ls_hits_reduced_high_stats_split.root";
 x.AddReducedEventsFile( simsfile.str() ); 

x.GenerateYieldDist(0.1,7.,0.18,0.95,0.95,13);

 minuit->SetFCN(&FCN);

 minuit->DefineParameter(0,"Scale",0.1,0.5,0.,0.5);
 minuit->DefineParameter(1,"Yield",7.,500.,2.,15.);
 minuit->DefineParameter(2,"ChanNum",8.,0.1,0.,16.);

 minuit->FixParameter(2);

 minuit->SetPrintLevel(1);  
 minuit->Migrad(); 

 return 0;

}

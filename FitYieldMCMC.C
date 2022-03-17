#include <string>
#include <sstream>
#include <iostream>

#include "XeNeuSimsAnalysisEnvironment.hh"

#include "TMath.h"
#include "TRandom3.h"
#include "TFile.h"

XeNeuSimsAnalysisEnvironment x;


TRandom3 r;



//////////////////////////////////////////////////////////////////////////////
void MetropolisHastingsSampler( int ch_num, int nsteps );
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char * argv[]) {

 if( argc != 5 ) {
     printf("Give channel number and yes/no gaussian fit (1/0) as argument.\n");
     printf("Usage:\n");
     printf("\t./FitYieldMCMC <channel_num> <gauss fit on/off> <fit_min> <fit_max>\n");
     return -1;
 }

 int chan_num = atoi(argv[1]);
 if( chan_num != 5 &&
     chan_num != 6 &&
     chan_num != 8 &&
     chan_num != 9 &&
     chan_num != 10 &&
     chan_num != 11 &&
     chan_num != 12 &&
     chan_num != 13 &&
     chan_num != 14 &&
     chan_num != 15 ) {
        printf("Invalid channel number.\n");
        return -1;
 }

//x.GenerateYieldDist(7.,1.,13);
x.SetELifetime( 158. ); // lifetime in microseconds
x.SetEEE( 0.958 ); // electron extraction efficiency
x.SetDriftVelocity( 1.75 ); // drift speed in mm/us
double se_size = 57.0;
x.SetSESize( se_size );
x.SetSEWidth( 11.6 ); // Absolute width of the pulse area 
                      // distribution for single electrons

 x.AddDataHistFiles( 101 );

 std::stringstream simsfile;
 simsfile << "/g/g20/lenardo1/Simulations/" <<
//             "AllOldFits_1-22-2019/reduced_data/" << 
//             "Fit_1-22-2019_18degRot_NoOffset/reduced_data_02/" << 
//             "tunl_sims_conical_ls_hits_reduced_7mmOffset_split_fidcut.root";
//             "tunl_sims_conical_ls_hits_reduced_18degRot_split_2mmOffsetfidcut_02.root";
//               "tunl_sims_conical_ls_hits_reduced_18degRot_split_fidcut.root";
//             "tunl_sims_conical_ls_hits_reduced_high_stats_batch1-17_fid_cut_split.root";
 "tunl_sims_conical_ls_hits_reduced_high_stats_and_18degRot_split_fidcut_5mmOffset_03.root";


 x.AddReducedEventsFile( simsfile.str() ); 


if( atoi(argv[2])==1 )
  x.SetDoGaussianFit(true);
else
  x.SetDoGaussianFit(false);

if( chan_num == 8 ) {
  x.SetApplyScalePrior(true);
  x.SetScaleMean( 0.10 );
  x.SetScaleSig(  0.015 );
}
if( chan_num == 15 ) {
  x.SetApplyScalePrior(true);
  x.SetScaleMean( 0.10 );
  x.SetScaleSig(  0.015 );
}

printf("Setting fit range:\n");
printf("Min: %f\t Max: %f\n",atof(argv[3]),atof(argv[4]));
x.SetFitMin( atof(argv[3])*se_size );
x.SetFitMax( atof(argv[4])*se_size );

MetropolisHastingsSampler( chan_num, 5000 );

 return 0;

}
///////////////////////////////////////////////////////////////////////////////
void MetropolisHastingsSampler( int ch_num, int nsteps ) {
   
   // Define variables, initial errors.
   double scale = 0.2;
   double scale_err = 0.01;
   double yield = 10.;
   double yield_err = 0.3;
   double width = 1.;
   double width_err = 0.2;

   if( ch_num == 15 || ch_num == 8 ) yield_err = 0.7;


   r.SetSeed(0);

   // Create the output file
//   char outfilename[100];
   std::stringstream outfilename;
   outfilename << "MCMC_samples_101_ch_" << ch_num << "_";
   if( x.GetDoGaussianFit() )
       outfilename << "GaussFit_";
       //sprintf(outfilename,"MCMC_samples_ch_%d_%s_%dsteps.root",ch_num,"GaussFit",nsteps);
   else
       outfilename << "PoissFit_";
       //sprintf(outfilename,"MCMC_samples_ch_%d_%s_%dsteps.root",ch_num,"PoissFit",nsteps);
   if( x.GetApplyScalePrior() )
       outfilename << "ScalePrior_";    
   else
       outfilename << "NoScalePrior_";
   outfilename << "fidcut_spect3_high_stats_and_18degRot_5mmOffset";
   outfilename << nsteps << "steps_03.root";

   TFile * outfile = new TFile(outfilename.str().c_str(),"recreate");

   // Create the output tree
   TTree * samples = new TTree();
   samples->SetName("samples");
   int step_num;
   double likelihood;


   samples->Branch("step_num",&step_num,"step_num/I");
   samples->Branch("scale",&scale,"scale/D");
   samples->Branch("yield",&yield,"yield/D");
   samples->Branch("width",&width,"width/D");   
   samples->Branch("likelihood",&likelihood,"likelihood/D");

   // First, get an initial likelihood value
   likelihood = (-1.)*x.ComputeNegativeLogLikelihood(scale,yield,width,ch_num);
   step_num = 0;
   samples->Fill();

   double currentLL = likelihood, 
          newLL,
          logProb,
          jumpProb,
          newYield,
          newScale,
          newWidth;

   int stepsAccepted = 0;
   // start stepping (counting from 1)
   for(int i=1; i<nsteps; i++) {
       if( i % 1000 == 0 ) { 
            printf("Step %d\t",i);
            printf("Accepted: %d\t",stepsAccepted);
            printf("Acceptance fraction: %f\n",(double)stepsAccepted/(double)i);
            printf("yield: %f newYield: %f scale: %f oldScale: %f oldLL: %f, newLL: %f \n",
               yield,newYield,scale,newScale,currentLL,newLL);
       }
       step_num = i;
       // Proposal is symmetric, uncorrelated gaussian
       newYield = r.Gaus(yield,yield_err);
       newScale = r.Gaus(scale,scale_err);
       if( x.GetDoGaussianFit() ) 
           newWidth = r.Gaus(width,width_err);
       else
           newWidth = width;        
       // Generate the new distribution, calculate the likelihood.
       currentLL = (-1.)*x.ComputeNegativeLogLikelihood(scale,yield,width,ch_num);       
       newLL = (-1.)*x.ComputeNegativeLogLikelihood(newScale,newYield,newWidth,ch_num);
       if( newLL > currentLL ){
           likelihood = newLL;
           currentLL = newLL;
           yield = newYield;
           scale = newScale;
           width = newWidth;
           stepsAccepted++;
       } else {
           logProb = newLL - currentLL;
           jumpProb = TMath::Exp(logProb);
           if( r.Uniform() < jumpProb ) {
              likelihood = newLL;
              currentLL = newLL;
              yield = newYield;
              scale = newScale;
              width = newWidth;
              stepsAccepted++;
           } 
       }
       samples->Fill();
   }

   samples->Write();
   outfile->Close();
   delete outfile;
}

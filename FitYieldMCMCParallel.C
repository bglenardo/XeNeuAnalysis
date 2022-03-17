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
void MetropolisHastingsSampler( int ch_num, int nsteps, int dataset, char* offset, 
                                char* drift_velocity, char* elifetime, char* EEE );
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char * argv[]) {

 if( argc != 10 ) {
     printf("Give channel number and yes/no gaussian fit (1/0) as argument.\n");
     printf("Usage:\n");
     printf("\t./FitYieldMCMC <channel_num> <gauss fit on/off> <fit_min> <fit_max> <dataset> <offset> <drift velocity> <elifetime> <EEE>\n");
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

 // Drift velocity check
 if( atof(argv[7]) < 1. || atof(argv[7]) > 3. ){
     printf("Drift velocity %fmm/us outside allowed range. Check the inputs.\n",atof(argv[7]));
     return -1;
 }
 // Lifetime check
 if( atof(argv[8]) < 70. && atof(argv[8]) > 1000. ){
     printf("Electron lifetime %fus outside allowed range. Check the inputs.\n",atof(argv[8]));
 }
 // Poisson/gauss fit check
 if( atoi(argv[2]) != 1 && atoi(argv[2]) != 0 ){
     printf("Invalid choice of %d for Poisson vs. Gaussian fit. Check the inputs\n",atoi(argv[2]));
 } 


//x.GenerateYieldDist(7.,1.,13);
x.SetEEE( atof(argv[9]) ); // electron extraction efficiency
double se_size = 57.0;
x.SetSESize( se_size );
x.SetSEWidth( 11.6 ); // Absolute width of the pulse area 
                      // distribution for single electrons

x.AddDataHistFiles( atoi(argv[5]) );
x.SetELifetime( atof(argv[8]) ); // lifetime in microseconds
x.SetDriftVelocity( atof(argv[7]) ); // drift speed in mm/us

 std::stringstream simsfile;
 simsfile << "/g/g20/lenardo1/Simulations/" <<
//             "AllOldFits_1-22-2019/reduced_data/" << 
//             "Fit_1-22-2019_18degRot_NoOffset/reduced_data_02/" << 
//             "tunl_sims_conical_ls_hits_reduced_7mmOffset_split_fidcut.root";
//             "tunl_sims_conical_ls_hits_reduced_18degRot_split_2mmOffsetfidcut_02.root";
//               "tunl_sims_conical_ls_hits_reduced_18degRot_split_fidcut.root";
//             "tunl_sims_conical_ls_hits_reduced_high_stats_batch1-17_fid_cut_split.root";
 "tunl_sims_conical_ls_hits_reduced_high_stats_and_18degRot_split_fidcut_" << argv[6] << "Offset_03.root";
 printf("Input file:\n %s\n",simsfile.str().c_str());

 x.AddReducedEventsFile( simsfile.str() ); 


if( atoi(argv[2])==1 )
  x.SetDoGaussianFit(true);
else
  x.SetDoGaussianFit(false);

//if( chan_num == 8 ) {
//  x.SetApplyScalePrior(true);
//  if( atoi(argv[5]) == 100 ) {
//      x.SetScaleMean( 0.0371 );
//      x.SetScaleSig(  0.0033 );
//  } else if( atoi(argv[5]) == 101 ) {
//      x.SetScaleMean( 0.0721 );
//      x.SetScaleSig(  0.0085 );
//  } else if( atoi(argv[5]) == 102 ) {
//      x.SetScaleMean( 0.0712 );
//      x.SetScaleSig(  0.0083 );
//  } else if( atoi(argv[5]) == 103 ) {
//      x.SetScaleMean( 0.0519 );
//      x.SetScaleSig(  0.0062 );
//  }
//}
//if( chan_num == 15 ) {
//  x.SetApplyScalePrior(true);
//  if( atoi(argv[5]) == 100 ) {
//      x.SetScaleMean( 0.0371 );
//      x.SetScaleSig(  0.0033 );
//  } else if( atoi(argv[5]) == 101 ) {
//      x.SetScaleMean( 0.0721 );
//      x.SetScaleSig(  0.0085 );
//  } else if( atoi(argv[5]) == 102 ) {
//      x.SetScaleMean( 0.0712 );
//      x.SetScaleSig(  0.0083 );
//  } else if( atoi(argv[5]) == 103 ) {
//      x.SetScaleMean( 0.0519 );
//      x.SetScaleSig(  0.0062 );
//  }
//}
if( chan_num == 8 ) {
  x.SetApplyScalePrior(true);
  if( atoi(argv[5]) == 100 ) {
      x.SetScaleMean( 0.0393 );
      x.SetScaleSig(  0.0033 );
  } else if( atoi(argv[5]) == 101 ) {
      x.SetScaleMean( 0.0793 );
      x.SetScaleSig(  0.0071 );
  } else if( atoi(argv[5]) == 102 ) {
      x.SetScaleMean( 0.0797 );
      x.SetScaleSig(  0.0057 );
  } else if( atoi(argv[5]) == 103 ) {
      x.SetScaleMean( 0.0580 );
      x.SetScaleSig(  0.0037 );
  }
}
if( chan_num == 15 ) {
  x.SetApplyScalePrior(true);
  if( atoi(argv[5]) == 100 ) {
      x.SetScaleMean( 0.0353 );
      x.SetScaleSig(  0.0019 );
  } else if( atoi(argv[5]) == 101 ) {
      x.SetScaleMean( 0.0668 );
      x.SetScaleSig(  0.0046 );
  } else if( atoi(argv[5]) == 102 ) {
      x.SetScaleMean( 0.0650 );
      x.SetScaleSig(  0.0021 );
  } else if( atoi(argv[5]) == 103 ) {
      x.SetScaleMean( 0.0473 );
      x.SetScaleSig(  0.0026 );
  }
}


printf("Setting fit range:\n");
printf("Min: %f\t Max: %f\n",atof(argv[3]),atof(argv[4]));
x.SetFitMin( atof(argv[3])*se_size );
x.SetFitMax( atof(argv[4])*se_size );

MetropolisHastingsSampler( chan_num, 5000, atoi(argv[5]), argv[6], argv[7], argv[8], argv[9] );

 return 0;

}
///////////////////////////////////////////////////////////////////////////////
void MetropolisHastingsSampler( int ch_num, int nsteps, int dataset, char* offset, 
                                char* drift_velocity, char * elifetime, char* EEE ) {
   
   // Define variables, initial errors.
   double scale = 0.06;
   double scale_err = 0.01;
   if( ch_num == 8 || ch_num == 15)   scale_err = 0.001;
   double yield = 4.7;
   double yield_err = 0.1;
   if( ch_num==8 || ch_num == 5 ) yield_err = 0.4;
   double width = 0.95;
   double width_err = 0.2;

   if( ch_num == 15 || ch_num == 8 ) yield_err = 0.5;


   r.SetSeed(0);

   // Create the output file
//   char outfilename[100];
   std::stringstream outfilename;
   outfilename << "MCMC_" << dataset 
               << "_ch_" << ch_num << "_";
   if( x.GetDoGaussianFit() )
       outfilename << "GaussFit_";
       //sprintf(outfilename,"MCMC_samples_ch_%d_%s_%dsteps.root",ch_num,"GaussFit",nsteps);
   else
       outfilename << "PoissFit_";
       //sprintf(outfilename,"MCMC_samples_ch_%d_%s_%dsteps.root",ch_num,"PoissFit",nsteps);
   if( x.GetApplyScalePrior() )
       outfilename << "ScalePriorLR_";    
   else
       outfilename << "NoScalePrior_";
   outfilename << "off_" << offset <<
                  "_DV_" << drift_velocity <<
                  "_EL_" << elifetime <<
                  "_EEE_" << EEE; 
   outfilename << "_03.root";
//   outfilename << nsteps << "steps_03.root";

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
       if( i % 100 == 0 ) { 
            printf("Step %d\t",i);
            printf("Accepted: %d\t",stepsAccepted);
            printf("Acceptance fraction: %f\n",(double)stepsAccepted/(double)i);
            printf("yield: %f newYield: %f scale: %f oldScale: %f oldLL: %f, newLL: %f \n",
               yield,newYield,scale,newScale,currentLL,newLL);
            fflush(stdout);
       }
       step_num = i;
       // Proposal is symmetric, uncorrelated gaussian
       newYield = r.Gaus(yield,yield_err);
       newScale = r.Gaus(scale,scale_err);
       if( x.GetDoGaussianFit() ) 
           // In the lowest energy channels, the width and 
           // yield are highly correlated for the gaussian
           // fit. So, we'll sample assuming a relative
           // slope of width = -0.5*yield_jump + yield
           if( ch_num == 8 || ch_num == 15){
             double sqrt5 = TMath::Sqrt(5.);
             double long_dim = r.Gaus(0.,0.7);
             double short_dim = r.Gaus(0.,0.1);
             newWidth = width - 1./sqrt5*long_dim + 2./sqrt5*short_dim;
             newYield = yield + 2./sqrt5*long_dim - 1./sqrt5*short_dim;
           } else { 
             newWidth = r.Gaus(width,width_err);
           }
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

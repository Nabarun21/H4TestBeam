#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TApplication.h"
#include "TChain.h"
#include "TTreeFormula.h"
#include "TCanvas.h"
#include "TVector.h"
#include "TLine.h"
#include "TGraph.h"

#include <iostream>
#include <fstream>

//---- from Event.hpp
#define MAX_ADC_CHANNELS 200
#define MAX_DIGI_SAMPLES 100000
#define MAX_TDC_CHANNELS 200
#define MAX_SCALER_WORDS 16
#define MAX_PATTERNS 16
#define MAX_PATTERNS_SHODO 16
#define SMALL_HODO_X_NFIBERS 8
#define SMALL_HODO_Y_NFIBERS 8
#define MAX_TRIG_WORDS 32
#define MAX_RO 10

//---- from TBtree Shashlik ----
#include "include/TBEvent.h"
#include "include/HodoCluster.h"
#include "include/TBRecHit.h"
#include "include/Mapper.h"

#include "include/CaloCluster.h"


#define hodoX1 0
#define hodoY1 1
#define hodoX2 2
#define hodoY2 3

//func to calculate mean

float mean(vector<float> x){

  float mean =0; 
  double sum =0;
  
  if (x.size()>0){
    
    for (int i=0; i<x.size(); i++){
      sum+=x[i];
    }
    
    mean=sum/x.size();
    
  }
    return mean;
} 



void mean_stdev(double& sd,double& mean,vector<float> x){
  //"""function to calculate mean and standard deviation of a vector"""

  double sum =0;
  double st_sum =0;

  if (x.size()>0){
    
    for (int i=0; i<x.size(); i++){
      sum+=x[i];
    }
    
    mean=sum/x.size();
    
    for (int i=0; i<x.size(); i++){
      st_sum+=pow((mean -x[i]),2);
    }
    
    sd=sqrt(st_sum/x.size());
  }
}


//---- draw shashlik matrix

void DrawShashlikModule(TPad* cc){
 
 cc->cd();
 for (int i=0; i<5; i++) {
//   TString name = Form ("vert_%d",i);
//   TLine* vert = new TLine (name.Data(),i*14-28,-28,i*14-28,28);
  TLine* vert = new TLine (i*14-28,-28,i*14-28,28);
  vert->SetLineColor(kRed);
  vert->DrawClone("same");
 }

 for (int i=0; i<5; i++) {
//   TString name = Form ("oriz_%d",i);
//   TLine* oriz = new TLine (name.Data(),-28,i*14-28,28,i*14-28);
  TLine* oriz = new TLine (-28,i*14-28,28,i*14-28);
  oriz->SetLineColor(kRed);
  oriz->DrawClone("same");
 }
 
}



//---- boost
// #include "boost/program_options.hpp"
// #include "boost/program_options/options_backescription.hpp"

// namespace po = boost::program_options;

#include <sstream> 


//---- transform map into vectors
void TransformFibers(std::map<std::pair<int,int>, int > fibers, std::vector <int>& fibers_X1, std::vector <int>& fibers_X2, std::vector <int>& fibers_Y1, std::vector <int>& fibers_Y2){
  
 std::pair<int,int> fibers_mappairY1;
 fibers_mappairY1.first  = hodoY1;
 std::pair<int,int> fibers_mappairY2;
 fibers_mappairY2.first  = hodoY2;
 
 //---- Y direction is inverted !?!?
 for(int iBinY=63;iBinY>=0;iBinY--){
//   for(int iBinY=0;iBinY<64;iBinY++){
  fibers_mappairY1.second = iBinY;
  fibers_mappairY2.second = iBinY;
  fibers_Y1.push_back( fibers[fibers_mappairY1] );
  fibers_Y2.push_back( fibers[fibers_mappairY2] );
 }

 std::pair<int,int> fibers_mappairX1;
 fibers_mappairX1.first  = hodoX1;
 std::pair<int,int> fibers_mappairX2;
 fibers_mappairX2.first  = hodoX2;
 for(int iBinX=0;iBinX<64;iBinX++){
  fibers_mappairX1.second = iBinX;
  fibers_mappairX2.second = iBinX;
  fibers_X1.push_back( fibers[fibers_mappairX1] );
  fibers_X2.push_back( fibers[fibers_mappairX2] );
 }
 
}



//---- Hodoscope clusters
std::vector<HodoCluster*> getHodoClusters( std::vector<int> hodo) {
 float fibreWidth = 0.5;
 int nClusterMax = 10;
 float Cut = 0;
 
 std::vector<HodoCluster*> clusters;
 HodoCluster* currentCluster = new HodoCluster( hodo.size(), fibreWidth );
 for ( unsigned i=0; i<hodo.size(); ++i ) {
  if ( hodo[i] > Cut) { // hit
   if ( currentCluster->getSize() < nClusterMax ) {
    currentCluster->addFibre( i );
   } else {
    clusters.push_back( currentCluster ); // store old one
    currentCluster = new HodoCluster( hodo.size(), fibreWidth ); // create a new one
    currentCluster->addFibre( i ); // get that fibre!
   }
  } else { // as soon as you find a hole
   if ( currentCluster->getSize() > 0 ) {
    clusters.push_back( currentCluster ); // store old one
    currentCluster = new HodoCluster( hodo.size(), fibreWidth ); // create a new one
   }
  }
 } // for fibres
 if ( currentCluster->getSize()>0 )
  clusters.push_back( currentCluster ); // store last cluster
  return clusters;
}


//---- Reconstruct Hodoscope clusters
 void doHodoReconstruction( std::vector<int> input_values, std::vector<int>& nFibres, std::vector<float>& cluster_position, float table) {
 std::vector<HodoCluster*> clusters = getHodoClusters( input_values );
 for( unsigned i=0; i<clusters.size(); ++i ) {
  nFibres.push_back( clusters[i]->getSize() );
  cluster_position.push_back( clusters[i]->getPosition() - table );
 }
}




//---- distance definition in calorimeter position
float DR (float x1, float x2, float y1, float y2) {
 return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

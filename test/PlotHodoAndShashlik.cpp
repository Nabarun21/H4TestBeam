#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TApplication.h"
#include "TChain.h"
#include "TTreeFormula.h"
#include "TCanvas.h"

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



#define hodoX1 0
#define hodoY1 1
#define hodoX2 2
#define hodoY2 3




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
 
 for(int iBinY=0;iBinY<64;iBinY++){
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
void doHodoReconstruction( std::vector<int> input_values, std::vector<int>& nFibres, std::vector<float>& cluster_position) {
 std::vector<HodoCluster*> clusters = getHodoClusters( input_values );
 for( unsigned i=0; i<clusters.size(); ++i ) {
  nFibres.push_back( clusters[i]->getSize() );
  cluster_position.push_back( clusters[i]->getPosition() );
 }
}



//---- Reconstruct Calorimeter clusters
void doCalorimeterReconstruction( Mapper* mapper,  std::vector<TBRecHit>* rechits, std::vector<float>& caloCluster_position_X, std::vector<float>& caloCluster_position_Y, std::vector<float>& caloCluster_Energy, int face) {
  
 //----      adc              X      Y
 std::map < float, std::pair<float, float> > map_of_calo_clusters;
 float max;
 int channelID;
 for (Int_t j = 0; j < rechits->size(); j++){
  TBRecHit &hit=rechits->at(j);
//   ped = hit.Pedestal();
//   sig = hit.NoiseRMS();
  max = hit.AMax();
//   maxTime = hit.TRise();
  channelID = hit.GetChannelID();
  
  if (max<0.1) continue; //---- very low threshold, just to remove noisy noise
  
  int moduleID,fiberID;
  mapper->ChannelID2ModuleFiber(channelID,moduleID,fiberID);  // get module and fiber IDs
  
  double x,y;
//   mapper->ModuleXY(moduleID,x,y);
  mapper->FiberXY(fiberID, x, y);
  
  //---- only the wanted face (front [1] or back [-1])
  if (moduleID<0 && face>0) {
   continue;
  }
  if (moduleID>0 && face<0) {
   continue;
  }
  
  std::pair<float, float> xy_pair;
  xy_pair.first = x;
  xy_pair.second = y;
  
  map_of_calo_clusters[-max] = xy_pair;
 }     
 
 //---- dump the map into the std::vector
//  for( std::map < float, std::pair<float, float> >::iterator ii=map_of_calo_clusters.begin(); ii!=map_of_calo_clusters.end(); ii++) {
//   caloCluster_Energy.push_back( ii->first );
//   caloCluster_position_X.push_back( ii->second.first );
//   caloCluster_position_Y.push_back( ii->second.second );
//  }

 
 
 //---- do clustering ----
//  int num_clusters = 0;
 float x_cluster = 0;
 float y_cluster = 0;
 float energy_cluster = 0;
 for( std::map < float, std::pair<float, float> >::iterator ii=map_of_calo_clusters.begin(); ii!=map_of_calo_clusters.end(); ii++) {
  x_cluster = x_cluster + (ii->second.first)  * (- ii->first) ;
  y_cluster = y_cluster + (ii->second.second) * (- ii->first) ;
  energy_cluster = energy_cluster - ii->first;
//   num_clusters++;
 }
 if (energy_cluster != 0) {
  x_cluster = x_cluster / energy_cluster;
  y_cluster = y_cluster / energy_cluster;
 }
 
 if (energy_cluster != 0) {
  caloCluster_Energy.push_back( energy_cluster );
  caloCluster_position_X.push_back( x_cluster );
  caloCluster_position_Y.push_back( y_cluster );
 }
 
//  std::cout << " num_clusters = " << num_clusters << std::endl;
} 
  
  
  
  
  


int main(int argc, char**argv){
 
 std::string input_file;
 int maxEvents = -1;
 
 //---- configuration
 
 int c;
 while ((c = getopt (argc, argv, "i:m:")) != -1)
  switch (c)
  {
   case 'i': //---- input
    input_file = string(optarg);
    break;
   case 'm':
    maxEvents =  atoi(optarg);
    break;
    
   case '?':
    if (optopt == 'i' || optopt == 'm')
     fprintf (stderr, "Option -%c requires an argument.\n", optopt);
    else if (isprint (optopt))
     fprintf (stderr, "Unknown option `-%c'.\n", optopt);
    else
     fprintf (stderr,
              "Unknown option character `\\x%x'.\n",
              optopt);
     return 1;
   default:
    exit (-1);
  }
  
  
  
  //---- get vector of files
  
  std::vector<std::string> input_files_vector;
  std::stringstream ss(input_file);
  
  std::string token_string;
  while(std::getline(ss, token_string, ',')) {
   std::cout << token_string << '\n';
   input_files_vector.push_back(token_string);
  }
  
  std::cout << " input files:" << std::endl;
  for (int i=0; i<input_files_vector.size(); i++) {
   std::cout << " beam: " << input_files_vector.at(i) << std::endl;
  }
  
  //---- configuration (end)
  
  
  //---- read file
  
  TChain* H4tree = new TChain("t1041");
  for (unsigned int i=0; i<input_files_vector.size(); i++) {
   H4tree->Add(input_files_vector.at(i).c_str());
  }
  
  
  //---- read file
  int nEntries = H4tree->GetEntries(); 
  std::cout << " nEntries = " << nEntries << std::endl;
  if (maxEvents != -1) nEntries = maxEvents>nEntries ? nEntries : maxEvents ;
  std::cout << " new nEntries = " << nEntries << std::endl;
  
  
  TBSpill* tbspill = new TBSpill();
  TBEvent* tbevent = new TBEvent();
  
  TBranch *branch_event = H4tree->GetBranch("tbevent");
  branch_event->SetAddress(&tbevent);
  TBranch *branch_spill = H4tree->GetBranch("tbspill");
  branch_spill->SetAddress(&tbspill);
  
  TApplication* gMyRootApp = new TApplication("My ROOT Application", &argc, argv);
  
  TCanvas* cc_X = new TCanvas ("cc_X","",800,600);
  TCanvas* cc_Y = new TCanvas ("cc_Y","",800,600);
  
  TH2F *hHS_HS1_Cal_front_X  = new TH2F("hHS_HS1_Cal_front_X", "Hodoscope 1 vs Cal front X", 64, -32, 32, 64, -32, 32);
  TH2F *hHS_HS2_Cal_front_X  = new TH2F("hHS_HS2_Cal_front_X", "Hodoscope 2 vs Cal front X", 64, -32, 32, 64, -32, 32);  
  TH2F *hHS_HS1_Cal_front_Y  = new TH2F("hHS_HS1_Cal_front_Y", "Hodoscope 1 vs Cal front Y", 64, -32, 32, 64, -32, 32);
  TH2F *hHS_HS2_Cal_front_Y  = new TH2F("hHS_HS2_Cal_front_Y", "Hodoscope 2 vs Cal front Y", 64, -32, 32, 64, -32, 32);
  
  TH2F *hHS_HS1_Cal_back_X  = new TH2F("hHS_HS1_Cal_back_X", "Hodoscope 1 vs Cal back X", 64, -32, 32, 64, -32, 32);
  TH2F *hHS_HS2_Cal_back_X  = new TH2F("hHS_HS2_Cal_back_X", "Hodoscope 2 vs Cal back X", 64, -32, 32, 64, -32, 32);  
  TH2F *hHS_HS1_Cal_back_Y  = new TH2F("hHS_HS1_Cal_back_Y", "Hodoscope 1 vs Cal back Y", 64, -32, 32, 64, -32, 32);
  TH2F *hHS_HS2_Cal_back_Y  = new TH2F("hHS_HS2_Cal_back_Y", "Hodoscope 2 vs Cal back Y", 64, -32, 32, 64, -32, 32);
  
  bool haverechits = false;
  std::vector<TBRecHit> *rechits=0;
  if(H4tree->GetListOfBranches()->FindObject("tbrechits")) {
   std::cout << " found rechits " << std::endl;
   H4tree->SetBranchAddress("tbrechits",&rechits);
   haverechits = true;
  }
  
  Mapper *mapper = Mapper::Instance();
  
  for (int i=0; i<nEntries; i++) {
   
   if ((i%100)==0) {
    std::cout <<  " entry: " << i << "::" << nEntries << std::endl;
   }
   
   H4tree->GetEntry(i);
   
   //---- calorimeter data
   if (i==0) mapper->SetEpoch(tbevent->GetTimeStamp());

   std::vector<float> caloCluster_position_X_front;
   std::vector<float> caloCluster_position_Y_front;
   std::vector<float> caloCluster_Energy_front;
   
   doCalorimeterReconstruction( mapper, rechits, caloCluster_position_X_front, caloCluster_position_Y_front, caloCluster_Energy_front, 1);    
    
   std::vector<float> caloCluster_position_X_back;
   std::vector<float> caloCluster_position_Y_back;
   std::vector<float> caloCluster_Energy_back;
   
   doCalorimeterReconstruction( mapper, rechits, caloCluster_position_X_back, caloCluster_position_Y_back, caloCluster_Energy_back, -1);
    
   
   
   //---- hodoscope data
   Hodoscope hsch = tbevent->GetHSChan();
   //   hsch.Dump();
   std::map<std::pair<int,int>, int > fibers = hsch.GetFibers();
   
   std::vector <int> fibers_X1;
   std::vector <int> fibers_X2;
   std::vector <int> fibers_Y1;
   std::vector <int> fibers_Y2;
   
   TransformFibers(fibers, fibers_X1, fibers_X2, fibers_Y1, fibers_Y2);
     
   
   std::vector <int> n_fibers_X1;
   std::vector <int> n_fibers_X2;
   std::vector <int> n_fibers_Y1;
   std::vector <int> n_fibers_Y2;
   
   std::vector <float> pos_fibers_X1;
   std::vector <float> pos_fibers_X2;
   std::vector <float> pos_fibers_Y1;
   std::vector <float> pos_fibers_Y2;
   
   
   doHodoReconstruction( fibers_X1, n_fibers_X1, pos_fibers_X1 );
   doHodoReconstruction( fibers_X2, n_fibers_X2, pos_fibers_X2 );
   doHodoReconstruction( fibers_Y1, n_fibers_Y1, pos_fibers_Y1 );
   doHodoReconstruction( fibers_Y2, n_fibers_Y2, pos_fibers_Y2 );
   
   
   
   //---- now merge and compare
   if (pos_fibers_X1.size() > 1) {
    
    for (int iCluster = 0; iCluster < pos_fibers_X1.size(); iCluster++) {
     for (int iCalo = 0; iCalo < caloCluster_position_X_front.size(); iCalo++) {
      std::cout << " caloCluster_position_X_front.at(" << iCalo << "), pos_fibers_X1.at(" << iCluster << ")  = " << caloCluster_position_X_front.at(iCalo) << "," <<  pos_fibers_X1.at(iCluster) << std::endl;
      hHS_HS1_Cal_front_X->Fill(caloCluster_position_X_front.at(iCalo), pos_fibers_X1.at(iCluster));
     }
     for (int iCalo = 0; iCalo < caloCluster_position_X_back.size(); iCalo++) {
      hHS_HS1_Cal_back_X->Fill(caloCluster_position_X_back.at(iCalo), pos_fibers_X1.at(iCluster));
     }
    }
    for (int iCluster = 0; iCluster < pos_fibers_X2.size(); iCluster++) {
     for (int iCalo = 0; iCalo < caloCluster_position_X_front.size(); iCalo++) {
      hHS_HS2_Cal_front_X->Fill(caloCluster_position_X_front.at(iCalo), pos_fibers_X2.at(iCluster));
     }
     for (int iCalo = 0; iCalo < caloCluster_position_X_back.size(); iCalo++) {
      hHS_HS2_Cal_back_X->Fill(caloCluster_position_X_back.at(iCalo), pos_fibers_X2.at(iCluster));
     }
    }
    
    for (int iCluster = 0; iCluster < pos_fibers_Y1.size(); iCluster++) {
     for (int iCalo = 0; iCalo < caloCluster_position_Y_front.size(); iCalo++) {
      hHS_HS1_Cal_front_Y->Fill(caloCluster_position_Y_front.at(iCalo), pos_fibers_Y1.at(iCluster));
     }
     for (int iCalo = 0; iCalo < caloCluster_position_Y_back.size(); iCalo++) {
      hHS_HS1_Cal_back_Y->Fill(caloCluster_position_Y_back.at(iCalo), pos_fibers_Y1.at(iCluster));
     }
    }
    for (int iCluster = 0; iCluster < pos_fibers_Y2.size(); iCluster++) {
     for (int iCalo = 0; iCalo < caloCluster_position_Y_front.size(); iCalo++) {
      hHS_HS2_Cal_front_Y->Fill(caloCluster_position_Y_front.at(iCalo), pos_fibers_Y2.at(iCluster));
     }
     for (int iCalo = 0; iCalo < caloCluster_position_Y_back.size(); iCalo++) {
      hHS_HS2_Cal_back_Y->Fill(caloCluster_position_Y_back.at(iCalo), pos_fibers_Y2.at(iCluster));
     }
    }
   }
    
  }
  
  
  TF1* fxy = new TF1 ("fxy","x",-20,20);

  cc_X->Divide(2,2);
  
  cc_X->cd(1)->SetGrid();
  hHS_HS1_Cal_front_X->Draw("colz");
  hHS_HS1_Cal_front_X->GetXaxis()->SetTitle("Calo front");
  hHS_HS1_Cal_front_X->GetYaxis()->SetTitle("X1");
  fxy->Draw("same");
  
  cc_X->cd(2)->SetGrid();
  hHS_HS2_Cal_front_X->Draw("colz");
  hHS_HS2_Cal_front_X->GetXaxis()->SetTitle("Calo front");
  hHS_HS2_Cal_front_X->GetYaxis()->SetTitle("X2");
  fxy->Draw("same");
  
  cc_X->cd(3)->SetGrid();
  hHS_HS1_Cal_back_X->Draw("colz");
  hHS_HS1_Cal_back_X->GetXaxis()->SetTitle("Calo back");
  hHS_HS1_Cal_back_X->GetYaxis()->SetTitle("X1");
  fxy->Draw("same");
  
  cc_X->cd(4)->SetGrid();
  hHS_HS2_Cal_back_X->Draw("colz");
  hHS_HS2_Cal_back_X->GetXaxis()->SetTitle("Calo back");
  hHS_HS2_Cal_back_X->GetYaxis()->SetTitle("X2");
  fxy->Draw("same");
  
  
 
  
  
  cc_Y->Divide(2,2);
  
  cc_Y->cd(1)->SetGrid();
  hHS_HS1_Cal_front_Y->Draw("colz");
  hHS_HS1_Cal_front_Y->GetXaxis()->SetTitle("Calo front");
  hHS_HS1_Cal_front_Y->GetYaxis()->SetTitle("Y1");
  fxy->Draw("same");
  
  cc_Y->cd(2)->SetGrid();
  hHS_HS2_Cal_front_Y->Draw("colz");
  hHS_HS2_Cal_front_Y->GetXaxis()->SetTitle("Calo front");
  hHS_HS2_Cal_front_Y->GetYaxis()->SetTitle("Y2");
  fxy->Draw("same");
  
  cc_Y->cd(3)->SetGrid();
  hHS_HS1_Cal_back_Y->Draw("colz");
  hHS_HS1_Cal_back_Y->GetXaxis()->SetTitle("Calo back");
  hHS_HS1_Cal_back_Y->GetYaxis()->SetTitle("Y1");
  fxy->Draw("same");
  
  cc_Y->cd(4)->SetGrid();
  hHS_HS2_Cal_back_Y->Draw("colz");
  hHS_HS2_Cal_back_Y->GetXaxis()->SetTitle("Calo back");
  hHS_HS2_Cal_back_Y->GetYaxis()->SetTitle("Y2");
  fxy->Draw("same");
  
  gMyRootApp->Run(); 
  
}



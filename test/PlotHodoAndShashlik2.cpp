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
#include <cmath>
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


//----- find 2 lowest numbers in a vector

void find2lowest(vector<float> v,float &a,float &b){

  float temp=v[0];
  int pos=0;
  for(int i=1;i<v.size();i++){
    if(v[i]<temp){
      temp=v[i];
      pos=i;}
  }
  a=temp;
  v[pos]=9999999;
  temp=999999;
  for(int i=1;i<v.size();i++){
    if(v[i]<temp){
      temp=v[i];
    }
  }
  b=temp;
}


//---- draw shashlik matrix


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






  

int main(int argc, char**argv){
 
 std::string input_file;
 int maxEvents = -1;
 int doFiber = 0;
 float table_x_reference = 200; //---- mm
 float table_y_reference = 350; //---- mm
 float table_x = 200; //---- mm
 float table_y = 350; //---- mm
 
 float w0 = 5.0;
 //---- configuration

 int c;
 while ((c = getopt (argc, argv, "i:m:f:x:y:w:")) != -1)
  switch (c)
  {
   case 'i': //---- input
    input_file = string(optarg);
    break;
   case 'm':
    maxEvents =  atoi(optarg);
    break;
   case 'f':
    doFiber =  atoi(optarg);
    break;
   case 'x':
    table_x =  atof(optarg);
    break;
   case 'y':
    table_y =  atof(optarg);
    break;
   case 'w':
    w0 =  atof(optarg);
    break;
    
    
    
   case '?':
    if (optopt == 'i' || optopt == 'm' || optopt == 'f' || optopt == 'x' || optopt == 'y' || optopt == 'w')
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

  std::cout << " Table: " << std::endl;
  std::cout << "   x = " << table_x << " mm " << std::endl;
  std::cout << "   y = " << table_y << " mm " << std::endl;
  
  
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
   std::cout << "  >> file: " << input_files_vector.at(i) << std::endl;
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
  
//   TBranch *branch_event = H4tree->GetBranch("tbevent");
//   branch_event->SetAddress(&tbevent);
//   TBranch *branch_spill = H4tree->GetBranch("tbspill");
//   branch_spill->SetAddress(&tbspill);
  
  H4tree->SetBranchAddress("tbevent", &tbevent);
  H4tree->SetBranchAddress("tbspill", &tbspill);

  
  
  TApplication* gMyRootApp = new TApplication("My ROOT Application", &argc, argv);
  
  TCanvas* cc_hodo_scat = new TCanvas ("cc_hodo_scat","",800,800);
  TCanvas* cc_hodo_scat2 = new TCanvas ("cc_hodo_scat2","",800,400);
  TCanvas* cc_calo_scat = new TCanvas ("cc_calo_scat","",800,400);

  TCanvas* cc_Energy_front = new TCanvas ("cc_Energy_front","",800,800);
  TCanvas* cc_Energy_back  = new TCanvas ("cc_Energy_back", "",800,800);
  
  TH1F *h_energy_front  = new TH1F("h_energy_front", "front log ei/etot", 100, 0, 15);
  TH1F *h_energy_back   = new TH1F("h_energy_back",  "back  log ei/etot", 100, 0, 15);

  TH2F *h_energy_front_size  = new TH2F("h_energy_front_size", "front log ei/etot vs cluster size", 100, 0, 15,17,0,17);
  TH2F *h_energy_back_size   = new TH2F("h_energy_back_size",  "back  log ei/etot vs cluster size", 100, 0, 15,17,0,17);
  
  TH1F *Energy_Cal_back   = new TH1F("Energy_Cal_back",  "Energy Cal back " , 200, 0, 200);
  TH1F *Energy_Cal_front  = new TH1F("Energy_Cal_front", "Energy Cal front" , 200, 0, 200);
  Energy_Cal_back->SetLineColor(kBlue);
  Energy_Cal_front->SetLineColor(kBlue);

  
  TH1F *Energy_Beam   = new TH1F("Energy_Beam",  "Energy Beam " , 200, 0, 200);
  Energy_Beam->SetLineColor(kRed);
  
  
  
  TCanvas* cc_X = new TCanvas ("cc_X","",800,600);
  TCanvas* cc_Y = new TCanvas ("cc_Y","",800,600);
  TCanvas* cc_sd = new TCanvas ("cc_sd","",800,600);
  
  TH2F *hHS_HS1_Cal_front_X  = new TH2F("hHS_HS1_Cal_front_X", "Hodoscope 1 vs Cal front X", 64, -32, 32, 64, -32, 32);
  TH2F *hHS_HS2_Cal_front_X  = new TH2F("hHS_HS2_Cal_front_X", "Hodoscope 2 vs Cal front X", 64, -32, 32, 64, -32, 32);  
  TH2F *hHS_HS1_Cal_front_Y  = new TH2F("hHS_HS1_Cal_front_Y", "Hodoscope 1 vs Cal front Y", 64, -32, 32, 64, -32, 32);
  TH2F *hHS_HS2_Cal_front_Y  = new TH2F("hHS_HS2_Cal_front_Y", "Hodoscope 2 vs Cal front Y", 64, -32, 32, 64, -32, 32);
  
  TH2F *hHS_HS1_Cal_back_X  = new TH2F("hHS_HS1_Cal_back_X", "Hodoscope 1 vs Cal back X", 64, -32, 32, 64, -32, 32);
  TH2F *hHS_HS2_Cal_back_X  = new TH2F("hHS_HS2_Cal_back_X", "Hodoscope 2 vs Cal back X", 64, -32, 32, 64, -32, 32);  
  TH2F *hHS_HS1_Cal_back_Y  = new TH2F("hHS_HS1_Cal_back_Y", "Hodoscope 1 vs Cal back Y", 64, -32, 32, 64, -32, 32);
  TH2F *hHS_HS2_Cal_back_Y  = new TH2F("hHS_HS2_Cal_back_Y", "Hodoscope 2 vs Cal back Y", 64, -32, 32, 64, -32, 32);
  
  TCanvas* cc_DX = new TCanvas ("cc_DX","",800,400);
  
  TH2F *X_hHS_HS1_Cal_front  = new TH2F("X_hHS_HS1_Cal_front", "X Hodoscope 1 vs Cal front ", 64, -32, 32, nEntries, 0, nEntries);
  TH2F *X_hHS_HS2_Cal_front  = new TH2F("X_hHS_HS2_Cal_front", "X Hodoscope 2 vs Cal front ", 64, -32, 32, nEntries, 0, nEntries);
  TH2F *X_hHS_HS1_Cal_back   = new TH2F("X_hHS_HS1_Cal_back",  "X Hodoscope 1 vs Cal back " , 64, -32, 32, nEntries, 0, nEntries);
  TH2F *X_hHS_HS2_Cal_back   = new TH2F("X_hHS_HS2_Cal_back",  "X Hodoscope 2 vs Cal back " , 64, -32, 32, nEntries, 0, nEntries);

  
  TH1F *X_h1_HS1_Cal_front  = new TH1F("X_h1_HS1_Cal_front", "X Hodoscope 1 vs Cal front ", 256, -32, 32);
  TH1F *X_h1_HS2_Cal_front  = new TH1F("X_h1_HS2_Cal_front", "X Hodoscope 2 vs Cal front ", 256, -32, 32);
  TH1F *X_h1_HS1_Cal_back   = new TH1F("X_h1_HS1_Cal_back",  "X Hodoscope 1 vs Cal back " , 256, -32, 32);
  TH1F *X_h1_HS2_Cal_back   = new TH1F("X_h1_HS2_Cal_back",  "X Hodoscope 2 vs Cal back " , 256, -32, 32);

  
  TH2F *sd_h1x  = new TH2F("sd_h1x", "hodo 1 rms v hits X ", 64, 0, 10, 64, 0, 20);
  TH2F *sd_h2x  = new TH2F("sd_h2x", "hodo 2 rms v hits X ", 64, 0, 10, 64, 0, 20);
  TH2F *sd_h1y  = new TH2F("sd_h1y", "hodo 1 rms v hits Y ", 64, 0, 10, 64, 0, 20);
  TH2F *sd_h2y  = new TH2F("sd_h2y", "hodo 2 rms v hits Y ", 64, 0, 10, 64, 0, 20);

  TCanvas* cc_DY = new TCanvas ("cc_DY","",800,400);
  
  TH2F *Y_hHS_HS1_Cal_front  = new TH2F("Y_hHS_HS1_Cal_front", "Y Hodoscope 1 vs Cal front ", 64, -32, 32, nEntries, 0, nEntries);
  TH2F *Y_hHS_HS2_Cal_front  = new TH2F("Y_hHS_HS2_Cal_front", "Y Hodoscope 2 vs Cal front ", 64, -32, 32, nEntries, 0, nEntries);
  TH2F *Y_hHS_HS1_Cal_back   = new TH2F("Y_hHS_HS1_Cal_back",  "Y Hodoscope 1 vs Cal back " , 64, -32, 32, nEntries, 0, nEntries);
  TH2F *Y_hHS_HS2_Cal_back   = new TH2F("Y_hHS_HS2_Cal_back",  "Y Hodoscope 2 vs Cal back " , 64, -32, 32, nEntries, 0, nEntries);
  
  
  TH1F *Y_h1_HS1_Cal_front  = new TH1F("Y_h1_HS1_Cal_front", "Y Hodoscope 1 vs Cal front ", 256, -32, 32);
  TH1F *Y_h1_HS2_Cal_front  = new TH1F("Y_h1_HS2_Cal_front", "Y Hodoscope 2 vs Cal front ", 256, -32, 32);
  TH1F *Y_h1_HS1_Cal_back   = new TH1F("Y_h1_HS1_Cal_back",  "Y Hodoscope 1 vs Cal back " , 256, -32, 32);
  TH1F *Y_h1_HS2_Cal_back   = new TH1F("Y_h1_HS2_Cal_back",  "Y Hodoscope 2 vs Cal back " , 256, -32, 32);
  
  
  
  TCanvas* cc_hodo = new TCanvas ("cc_hodo","",800,800);
  
  TH2F *hHS_HS2_HS1_X  = new TH2F("hHS_HS2_HS1_X", "Hodoscope 2 vs Hodoscope 1 X", 64, -32, 32, 64, -32, 32);
  TH2F *hHS_HS2_HS1_Y  = new TH2F("hHS_HS2_HS1_Y", "Hodoscope 2 vs Hodoscope 1 Y", 64, -32, 32, 64, -32, 32);
    
  TH2F *hHS_HS1  = new TH2F("hHS_HS1", "Hodoscope 1", 64, -32, 32, 64, -32, 32);
  TH2F *hHS_HS2  = new TH2F("hHS_HS2", "Hodoscope 2", 64, -32, 32, 64, -32, 32);
  
  TCanvas* cc_Cal = new TCanvas ("cc_Cal","",800,800);
  
  TH2F *hHS_Cal_front;
  TH2F *hHS_Cal_back;

  TH2F *hHS_Cal_front_module;
  TH2F *hHS_Cal_back_module;
  


  if (doFiber) {
   hHS_Cal_front  = new TH2F("hHS_Cal_front", "Cal front ", 48, -28, 28, 48, -28, 28);
   hHS_Cal_back   = new TH2F("hHS_Cal_back",  "Cal back " , 48, -28, 28, 48, -28, 28);
   
   hHS_Cal_front_module  = new TH2F("hHS_Cal_front_module", "Cal front ", 24, -28, 28, 24, -28, 28);
   hHS_Cal_back_module   = new TH2F("hHS_Cal_back_module",  "Cal back " , 24, -28, 28, 24, -28, 28);
  }
  else {
   hHS_Cal_front  = new TH2F("hHS_Cal_front", "Cal front ", 72, -28, 28, 72, -28, 28);
   hHS_Cal_back   = new TH2F("hHS_Cal_back",  "Cal back " , 72, -28, 28, 72, -28, 28);
   
   hHS_Cal_front_module  = new TH2F("hHS_Cal_front_module", "Cal front ", 12, -28, 28, 12, -28, 28);
   hHS_Cal_back_module   = new TH2F("hHS_Cal_back_module",  "Cal back " , 12, -28, 28, 12, -28, 28);
  }
  
  
  
  
  bool haverechits = false;
  std::vector<TBRecHit> *rechits=0;
  if(H4tree->GetListOfBranches()->FindObject("tbrechits")) {
   std::cout << " found rechits " << std::endl;
   H4tree->SetBranchAddress("tbrechits",&rechits);
   haverechits = true;
  }
  
//   Mapper *mapper = Mapper::Instance();
  
  CaloCluster* caloCluster = new CaloCluster();
  caloCluster->setW0(w0);
  caloCluster->setInterCalibrationConstants("data/InterCalibrationConstants.txt");
  
    //new
  std::vector<float> h1x;
  std::vector<float> h1y;
  std::vector<float> h2x;
  std::vector<float> h2y;

  std::vector<float> calox;
  std::vector<float> caloy;

  std::vector<float> calox_bk;
  std::vector<float> caloy_bk;


  
  ///---- loop ----
  for (int i=0; i<nEntries; i++) {
   
   if ((i%1000)==0) {
     // std::cout <<  " entry: " << i << "::" << nEntries << std::endl;
   }
   
   H4tree->GetEntry(i);
   
  
   float table_x_shift = tbspill->GetTableX();
   float table_y_shift = tbspill->GetTableY();
   
   if ( ((table_x_reference - table_x_shift) != table_x) || (i == 0) ) {
     //  std::cout << " Table: " << std::endl;
     // std::cout << "   x = " << table_x_reference - table_x_shift << " mm " << std::endl;
     // std::cout << "   y = " << table_y_reference - table_y_shift << " mm " << std::endl;
   }
   
   table_x = table_x_reference - table_x_shift;
   table_y = table_y_reference - table_y_shift;
   
//    if (i == 0) {
//     std::cout << " Table: " << std::endl;
//     std::cout << "   x = " << table_x << " mm " << std::endl;
//     std::cout << "   y = " << table_y << " mm " << std::endl;
//    }
   
   
   Energy_Beam->Fill(tbspill->GetMomentum());

   
   
   //---- calorimeter data
   if (i==0) caloCluster->setMapperEpoch(tbevent->GetTimeStamp());
   
   
   std::vector<float> caloCluster_position_X_front;
   std::vector<float> caloCluster_position_Y_front;
   std::vector<float> caloCluster_Energy_front;
   std::vector<float> caloCluster_Energies_front;
   
   caloCluster->doCalorimeterReconstruction( rechits, 1, 60, doFiber);
   
   caloCluster_position_X_front.push_back( caloCluster->getPositionX() );
   caloCluster_position_Y_front.push_back( caloCluster->getPositionY() );
   caloCluster_Energy_front.push_back( caloCluster->getEnergy() );
   caloCluster_Energies_front = caloCluster->getCaloClusterComponents();
   
//    std::cout << " caloCluster_Energies_front.size()= " << caloCluster_Energies_front.size() << std::endl;
   for (int iComponent = 0; iComponent<caloCluster_Energies_front.size(); iComponent++) {
    h_energy_front->Fill(-log(caloCluster_Energies_front.at(iComponent) / caloCluster_Energy_front.at(0)));
    h_energy_front_size->Fill(-log(caloCluster_Energies_front.at(iComponent) / caloCluster_Energy_front.at(0)),caloCluster_Energies_front.size());
   }

   Energy_Cal_front->Fill(caloCluster_Energy_front.at(0));
   
   
   
   std::vector<float> caloCluster_position_X_back;
   std::vector<float> caloCluster_position_Y_back;
   std::vector<float> caloCluster_Energy_back;
   std::vector<float> caloCluster_Energies_back;;
   
     
   caloCluster->doCalorimeterReconstruction( rechits, -1, 60, doFiber);
   
   caloCluster_position_X_back.push_back( caloCluster->getPositionX() );
   caloCluster_position_Y_back.push_back( caloCluster->getPositionY() );
   caloCluster_Energy_back.push_back( caloCluster->getEnergy() );
   caloCluster_Energies_back = caloCluster->getCaloClusterComponents();
   
   for (int iComponent = 0; iComponent<caloCluster_Energies_back.size(); iComponent++) {
    h_energy_back->Fill(-log(caloCluster_Energies_back.at(iComponent) / caloCluster_Energy_front.at(0)));
    h_energy_back_size->Fill(-log(caloCluster_Energies_back.at(iComponent) / caloCluster_Energy_front.at(0)),caloCluster_Energies_back.size());
   }
   Energy_Cal_back->Fill(caloCluster_Energy_back.at(0));
   
   
   //---- modular level DR = 5 mm
   std::vector<float> caloCluster_position_X_front_module;
   std::vector<float> caloCluster_position_Y_front_module;
   std::vector<float> caloCluster_Energy_front_module;
   
   caloCluster->doCalorimeterReconstruction( rechits, 1, 5, doFiber);

   caloCluster_position_X_front_module.push_back( caloCluster->getPositionX() );
   caloCluster_position_Y_front_module.push_back( caloCluster->getPositionY() );
   caloCluster_Energy_front_module.push_back( caloCluster->getEnergy() );
      
   
   std::vector<float> caloCluster_position_X_back_module;
   std::vector<float> caloCluster_position_Y_back_module;
   std::vector<float> caloCluster_Energy_back_module;
   
   caloCluster->doCalorimeterReconstruction( rechits, -1, 5, doFiber);
   
   caloCluster_position_X_back_module.push_back( caloCluster->getPositionX() );
   caloCluster_position_Y_back_module.push_back( caloCluster->getPositionY() );
   caloCluster_Energy_back_module.push_back( caloCluster->getEnergy() );


   
   
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
   
   std::vector <float> pos_fibers_X1; //---- units is mm
   std::vector <float> pos_fibers_X2; //---- units is mm
   std::vector <float> pos_fibers_Y1; //---- units is mm
   std::vector <float> pos_fibers_Y2; //---- units is mm
   
   
   //----                                                  table position in mm
   doHodoReconstruction( fibers_X1, n_fibers_X1, pos_fibers_X1, (table_x - 200)); //---- change of coordinate system using numbers from googledoc
   doHodoReconstruction( fibers_X2, n_fibers_X2, pos_fibers_X2, (table_x - 200)); //---- change of coordinate system using numbers from googledoc
   doHodoReconstruction( fibers_Y1, n_fibers_Y1, pos_fibers_Y1, (table_y - 350)); //---- change of coordinate system using numbers from googledoc
   doHodoReconstruction( fibers_Y2, n_fibers_Y2, pos_fibers_Y2, (table_y - 350)); //---- change of coordinate system using numbers from googledoc




   //---- just hodoscope information
   for (int iCluster1 = 0; iCluster1 < pos_fibers_X1.size(); iCluster1++) {
    for (int iCluster2 = 0; iCluster2 < pos_fibers_X2.size(); iCluster2++) {
     hHS_HS2_HS1_X->Fill(pos_fibers_X1.at(iCluster1),pos_fibers_X2.at(iCluster2));
//      std::cout << " pos_fibers_X:: " << pos_fibers_X1.at(iCluster1) << " :: " << pos_fibers_X2.at(iCluster2) << std::endl;
    }
   }
   
   for (int iCluster1 = 0; iCluster1 < pos_fibers_Y1.size(); iCluster1++) {
    for (int iCluster2 = 0; iCluster2 < pos_fibers_Y2.size(); iCluster2++) {
     hHS_HS2_HS1_Y->Fill(pos_fibers_Y1.at(iCluster1),pos_fibers_Y2.at(iCluster2));
    }
   }
   
   for (int iCluster1 = 0; iCluster1 < pos_fibers_X1.size(); iCluster1++) {
    for (int iCluster2 = 0; iCluster2 < pos_fibers_Y1.size(); iCluster2++) {
     hHS_HS1->Fill(pos_fibers_X1.at(iCluster1),pos_fibers_Y1.at(iCluster2));
    }
   }

   for (int iCluster1 = 0; iCluster1 < pos_fibers_X2.size(); iCluster1++) {
    for (int iCluster2 = 0; iCluster2 < pos_fibers_Y2.size(); iCluster2++) {
     hHS_HS2->Fill(pos_fibers_X2.at(iCluster1),pos_fibers_Y2.at(iCluster2));
    }
   }
   
   
  
   
   //---- now merge and compare
   if (pos_fibers_X1.size() ==1 ) {
   
      double sd=0;
       double ave=0;
       mean_stdev(sd,ave,pos_fibers_X1);
    
     //---- X
          
     for (int iCalo = 0; iCalo < caloCluster_position_X_front.size(); iCalo++) {
       
	 hHS_HS1_Cal_front_X->Fill(caloCluster_position_X_front.at(iCalo), ave);
	 X_hHS_HS1_Cal_front->Fill(caloCluster_position_X_front.at(iCalo) - ave, i);
	 X_h1_HS1_Cal_front->Fill(caloCluster_position_X_front.at(iCalo) - ave);
	 sd_h1x->Fill(pos_fibers_X1.size(),sd);
       
	 h1x.push_back(ave);
	 calox.push_back(caloCluster_position_X_front.at(iCalo));
     }   
   

           
     
     
     for (int iCalo = 0; iCalo < caloCluster_position_X_back.size(); iCalo++) {
             
       hHS_HS1_Cal_back_X->Fill(caloCluster_position_X_back.at(iCalo), ave);
       X_hHS_HS1_Cal_back->Fill(caloCluster_position_X_back.at(iCalo) - ave, i);
       X_h1_HS1_Cal_back->Fill(caloCluster_position_X_back.at(iCalo) - ave);
       
       mean_stdev(sd,ave,caloCluster_position_X_back);	  
       calox_bk.push_back(ave);
       
     }
   }
   
   else{
     if (pos_fibers_X1.size()!=0){
       float a=0;
       float b=0;
       double sd=0;
       find2lowest(pos_fibers_X1,a,b);
       double ave=(a+b)/2;
       if(abs(a-b)<1){
	 for (int iCalo = 0; iCalo < caloCluster_position_X_front.size(); iCalo++) {
	   
	   hHS_HS1_Cal_front_X->Fill(caloCluster_position_X_front.at(iCalo), ave);
	   X_hHS_HS1_Cal_front->Fill(caloCluster_position_X_front.at(iCalo) - ave, i);
	   X_h1_HS1_Cal_front->Fill(caloCluster_position_X_front.at(iCalo) - ave);
	   cout<<"1"<<endl;
	   h1x.push_back(ave);
	   calox.push_back(caloCluster_position_X_front.at(iCalo));
	 }
	 
	 
	 
	 for (int iCalo = 0; iCalo < caloCluster_position_X_back.size(); iCalo++) {
      
	   hHS_HS1_Cal_back_X->Fill(caloCluster_position_X_back.at(iCalo), ave);
	   X_hHS_HS1_Cal_back->Fill(caloCluster_position_X_back.at(iCalo) - ave, i);
	   X_h1_HS1_Cal_back->Fill(caloCluster_position_X_back.at(iCalo) - ave);
	   
	   mean_stdev(sd,ave,caloCluster_position_X_back);	  
	   calox_bk.push_back(ave);
	 }
       }
     }
   }
   
   //--------------------------------------------   	
   if (pos_fibers_X2.size()==1)
     {
       double sd=0;       
       double ave=0;
     
     for (int iCalo = 0; iCalo < caloCluster_position_X_front.size(); iCalo++) {

	 mean_stdev(sd,ave,pos_fibers_X2);       
	 sd_h2x->Fill(pos_fibers_X2.size(),sd);       
	 hHS_HS2_Cal_front_X->Fill(caloCluster_position_X_front.at(iCalo), ave);
	 X_hHS_HS2_Cal_front->Fill(caloCluster_position_X_front.at(iCalo) - ave, i);
	 X_h1_HS2_Cal_front->Fill(caloCluster_position_X_front.at(iCalo) - ave);
	 h2x.push_back(ave);
     }
     for (int iCalo = 0; iCalo < caloCluster_position_X_back.size(); iCalo++) {
	 hHS_HS2_Cal_back_X->Fill(caloCluster_position_X_back.at(iCalo), ave);
	 X_hHS_HS2_Cal_back->Fill(caloCluster_position_X_back.at(iCalo) - ave, i);
	 X_h1_HS2_Cal_back->Fill(caloCluster_position_X_back.at(iCalo) - ave);
     }

     }
   else{
       if (pos_fibers_X2.size()!=0){
	 float a=0;
	 float b=0;
	 double sd=0;
	 find2lowest(pos_fibers_X2,a,b);
	 double ave=(a+b)/2;       
	 if(abs(a-b)<1){
	 cout<<"2"<<endl;	   
	   for (int iCalo = 0; iCalo < caloCluster_position_X_front.size(); iCalo++) {

	       hHS_HS2_Cal_front_X->Fill(caloCluster_position_X_front.at(iCalo), ave);
	       X_hHS_HS2_Cal_front->Fill(caloCluster_position_X_front.at(iCalo) - ave, i);
	       X_h1_HS2_Cal_front->Fill(caloCluster_position_X_front.at(iCalo) - ave);
	 h2x.push_back(ave);
	   }

	   for (int iCalo = 0; iCalo < caloCluster_position_X_back.size(); iCalo++) {
	     hHS_HS2_Cal_back_X->Fill(caloCluster_position_X_back.at(iCalo), ave);
	       X_hHS_HS2_Cal_back->Fill(caloCluster_position_X_back.at(iCalo) - ave, i);
	       X_h1_HS2_Cal_back->Fill(caloCluster_position_X_back.at(iCalo) - ave);
	   }
	 }
       }
   }

   
   if (pos_fibers_Y1.size()==1)     {
     
     for (int iCalo = 0; iCalo < caloCluster_position_Y_front.size(); iCalo++) {
       
       double sd=0;
       double ave=0;
       
       mean_stdev(sd,ave,pos_fibers_Y1);       
       sd_h1y->Fill(pos_fibers_Y1.size(),sd);       
       
       if (sd<50){
      
	 hHS_HS1_Cal_front_Y->Fill(caloCluster_position_Y_front.at(iCalo), ave);
	 Y_hHS_HS1_Cal_front->Fill(caloCluster_position_Y_front.at(iCalo) - ave, i);
	 Y_h1_HS1_Cal_front->Fill(caloCluster_position_Y_front.at(iCalo) - ave);
       }
       
       h1y.push_back(ave);
       caloy.push_back(caloCluster_position_Y_front.at(iCalo)); }
     
     for (int iCalo = 0; iCalo < caloCluster_position_Y_back.size(); iCalo++) {
       
       double sd=0;
       double ave=0;
       
       mean_stdev(sd,ave,pos_fibers_Y1);       
       if (sd<50){
	 hHS_HS1_Cal_back_Y->Fill(caloCluster_position_Y_back.at(iCalo), ave);
	 Y_hHS_HS1_Cal_back->Fill(caloCluster_position_Y_back.at(iCalo) - ave, i);
	 Y_h1_HS1_Cal_back->Fill(caloCluster_position_Y_back.at(iCalo) - ave);
       }
       mean_stdev(sd,ave,caloCluster_position_Y_back);	  
        caloy_bk.push_back(ave);
     }
   }
   else{
     if (pos_fibers_Y1.size()!=0){
       float a=0;
       float b=0;
       double sd=0;
       find2lowest(pos_fibers_Y1,a,b);
       double ave=(a+b)/2;
	 if(abs(a-b)<1){
	 cout<<"3"<<endl;	   
	   for (int iCalo = 0; iCalo < caloCluster_position_Y_front.size(); iCalo++) {
	       hHS_HS1_Cal_front_Y->Fill(caloCluster_position_Y_front.at(iCalo), ave);
	       Y_hHS_HS1_Cal_front->Fill(caloCluster_position_Y_front.at(iCalo) - ave, i);
	       Y_h1_HS1_Cal_front->Fill(caloCluster_position_Y_front.at(iCalo) - ave);
	       h1y.push_back(ave);
	       caloy.push_back(caloCluster_position_Y_front.at(iCalo)); }
	   
	     for (int iCalo = 0; iCalo < caloCluster_position_Y_back.size(); iCalo++) {

		 hHS_HS1_Cal_back_Y->Fill(caloCluster_position_Y_back.at(iCalo), ave);
		 Y_hHS_HS1_Cal_back->Fill(caloCluster_position_Y_back.at(iCalo) - ave, i);
		 Y_h1_HS1_Cal_back->Fill(caloCluster_position_Y_back.at(iCalo) - ave);
	       mean_stdev(sd,ave,caloCluster_position_Y_back);	  
	       caloy_bk.push_back(ave);
	     }
	 }
     }
   }
     
   if (pos_fibers_Y2.size()==1){
       double ave=0;
       double sd=0;
       mean_stdev(sd,ave,pos_fibers_Y2);       
     
     for (int iCalo = 0; iCalo < caloCluster_position_Y_front.size(); iCalo++) {
       sd_h2y->Fill(pos_fibers_Y2.size(),sd);       
       hHS_HS2_Cal_front_Y->Fill(caloCluster_position_Y_front.at(iCalo), ave);
       Y_hHS_HS2_Cal_front->Fill(caloCluster_position_Y_front.at(iCalo) - ave, i);
       Y_h1_HS2_Cal_front->Fill(caloCluster_position_Y_front.at(iCalo) - ave);
       h2y.push_back(ave);
	 
       
     }
     for (int iCalo = 0; iCalo < caloCluster_position_Y_back.size(); iCalo++) {
       hHS_HS2_Cal_back_Y->Fill(caloCluster_position_Y_back.at(iCalo), ave);
       Y_hHS_HS2_Cal_back->Fill(caloCluster_position_Y_back.at(iCalo) - ave, i);
       Y_h1_HS2_Cal_back->Fill(caloCluster_position_Y_back.at(iCalo) - ave);
     }
   }
   else{
     if (pos_fibers_Y1.size()!=0){
       float a=0;
       float b=0;
       double sd =0;
       double ave=(a+b)/2;
       find2lowest(pos_fibers_Y1,a,b);
       if(abs(a-b)<1){ 
	 cout<<"4"<<endl;	   
	 for (int iCalo = 0; iCalo < caloCluster_position_Y_front.size(); iCalo++) {
		 hHS_HS2_Cal_front_Y->Fill(caloCluster_position_Y_front.at(iCalo), ave);
		 Y_hHS_HS2_Cal_front->Fill(caloCluster_position_Y_front.at(iCalo) - ave, i);
		 Y_h1_HS2_Cal_front->Fill(caloCluster_position_Y_front.at(iCalo) - ave);
		 h2y.push_back(ave);
	 }
	 for (int iCalo = 0; iCalo < caloCluster_position_Y_back.size(); iCalo++) {
	     hHS_HS2_Cal_back_Y->Fill(caloCluster_position_Y_back.at(iCalo), ave);
	     Y_hHS_HS2_Cal_back->Fill(caloCluster_position_Y_back.at(iCalo) - ave, i);
	     Y_h1_HS2_Cal_back->Fill(caloCluster_position_Y_back.at(iCalo) - ave);
	 }
	 
       } 
     }
   }


   //---- Fill Calorimeter mapping
   for (int iCaloX = 0; iCaloX < caloCluster_position_X_front.size(); iCaloX++) {
    for (int iCaloY = 0; iCaloY < caloCluster_position_Y_front.size(); iCaloY++) {
     hHS_Cal_front->Fill(caloCluster_position_X_front.at(iCaloX),caloCluster_position_Y_front.at(iCaloY));
    }
   }
   for (int iCaloX = 0; iCaloX < caloCluster_position_X_back.size(); iCaloX++) {
    for (int iCaloY = 0; iCaloY < caloCluster_position_Y_back.size(); iCaloY++) {
     hHS_Cal_back->Fill(caloCluster_position_X_back.at(iCaloX),caloCluster_position_Y_back.at(iCaloY));
    }
   }
   
   //---- Fill Calorimeter mapping with module information only
   for (int iCaloX = 0; iCaloX < caloCluster_position_X_front_module.size(); iCaloX++) {
    for (int iCaloY = 0; iCaloY < caloCluster_position_Y_front_module.size(); iCaloY++) {
     hHS_Cal_front_module->Fill(caloCluster_position_X_front_module.at(iCaloX),caloCluster_position_Y_front_module.at(iCaloY));
//      std::cout << " calo cluster position (module) : " << caloCluster_position_X_front_module.at(iCaloX) << " , " << caloCluster_position_Y_front_module.at(iCaloY) << std::endl;
    }
   }
   for (int iCaloX = 0; iCaloX < caloCluster_position_X_back_module.size(); iCaloX++) {
    for (int iCaloY = 0; iCaloY < caloCluster_position_Y_back_module.size(); iCaloY++) {
     hHS_Cal_back_module->Fill(caloCluster_position_X_back_module.at(iCaloX),caloCluster_position_Y_back_module.at(iCaloY));
    }
   }
}


  float *ah1x=&(h1x[0]);
  float *acalox=&(calox[0]);
  TGraph* Tg=new TGraph(h1x.size(),ah1x,acalox);
  
  
  float *ah2x=&(h2x[0]);
  TGraph* Tg2=new TGraph(h2x.size(),ah2x,acalox);
  
  
  float *ah1y=&(h1y[0]);
  float *acaloy=&(caloy[0]);
  TGraph* Tg3=new TGraph(h1y.size(),ah1y,acaloy);
  
  float *ah2y=&(h2y[0]);
  TGraph* Tg4=new TGraph(h2y.size(),ah2y,acaloy);
  
  TGraph* Tg5=new TGraph(h1x.size(),ah1x,ah1y);

  TGraph* Tg6=new TGraph(h1y.size(),ah2x,ah2y);

  float *acaloy_bk=&(caloy_bk[0]);
  float *acalox_bk=&(calox_bk[0]);


  TGraph* Tg7=new TGraph(calox.size(),acalox,acaloy);
  TGraph* Tg8=new TGraph(calox_bk.size(),acalox_bk,acaloy_bk);




  TPad* tempPad;
  
  //---- plot ----
/*     
  cc_calo_scat->Divide(2,1);
  
  cc_calo_scat->cd(1)->SetGrid();
  Tg7->Draw("ap");
  Tg7->GetXaxis()->SetTitle("calo_x_front");  
  Tg7->GetYaxis()->SetTitle("calo_y_front");  
  Tg7->SetTitle("CALOX VS CALOY_front");

  cc_calo_scat->cd(2)->SetGrid();
  Tg8->Draw("ap");
  Tg8->GetXaxis()->SetTitle("calo_x_back");  
  Tg8->GetYaxis()->SetTitle("calo_y_back");  
  Tg8->SetTitle("CALOX VS CALOY_back");

  cc_calo_scat->SaveAs("tempplot1/xcalovycalo.pdf");


  cc_hodo_scat2->Divide(2,1);
  
  cc_hodo_scat2->cd(1)->SetGrid();
  Tg5->Draw("ap");
  Tg5->GetXaxis()->SetTitle("hodo1_x");  
  Tg5->GetYaxis()->SetTitle("hodo1_y");  
  Tg5->SetTitle("HODO1X VS HODO1Y");

  cc_hodo_scat2->cd(2)->SetGrid();
  Tg6->Draw("ap");
  Tg6->GetXaxis()->SetTitle("hodo2_x");  
  Tg6->GetYaxis()->SetTitle("hodo2_y");  
  Tg6->SetTitle("HODO2X VS HODO2Y");

  cc_hodo_scat->Divide(2,2);

  cc_hodo_scat->cd(1)->SetGrid();
  Tg->Draw("ap");
  Tg->GetXaxis()->SetTitle("hodo1_x");  
  Tg->GetYaxis()->SetTitle("calo_x");  
  Tg->SetTitle("HODO1X VS CALOX");

   
  cc_hodo_scat->cd(2)->SetGrid();
  Tg2->Draw("ap");
  Tg2->GetXaxis()->SetTitle("hodo2_x");  
  Tg2->GetYaxis()->SetTitle("calo_x");  
  Tg2->SetTitle("HODO2X VS CALOX");

   
 
  cc_hodo_scat->cd(3)->SetGrid();
  Tg3->Draw("ap");
  Tg3->GetXaxis()->SetTitle("hodo1_y");  
  Tg3->GetYaxis()->SetTitle("calo_y");  
  Tg3->SetTitle("HODO1Y VS CALOY");

   

  cc_hodo_scat->cd(4)->SetGrid();
  Tg4->Draw("ap");
  Tg4->GetXaxis()->SetTitle("hodo2_y");  
  Tg4->GetYaxis()->SetTitle("calo_y");  
  Tg4->SetTitle("HODO2Y VS CALOY");
  
*/

  
    
  TF1* fxy = new TF1 ("fxy","x",-20,20);
  
  cc_hodo->cd(1)->SetGrid();
  hHS_HS2_HS1_X->Draw("colz");
  hHS_HS2_HS1_X->GetXaxis()->SetTitle("X1");
  hHS_HS2_HS1_X->GetYaxis()->SetTitle("X2");
  fxy->Draw("same");
  
  cc_hodo->cd(2)->SetGrid();
  hHS_HS2_HS1_Y->Draw("colz");
  hHS_HS2_HS1_Y->GetXaxis()->SetTitle("Y1");
  hHS_HS2_HS1_Y->GetYaxis()->SetTitle("Y2");
  fxy->Draw("same");
  
  
  cc_hodo->cd(3)->SetGrid();
  hHS_HS1->Draw("colz");
  hHS_HS1->GetXaxis()->SetTitle("X");
  hHS_HS1->GetYaxis()->SetTitle("Y");
  //  tempPad = (TPad*) gPad;
  // DrawShashlikModule(tempPad);
  
  cc_hodo->cd(4)->SetGrid();
  hHS_HS2->Draw("colz");
  hHS_HS2->GetXaxis()->SetTitle("X");
  hHS_HS2->GetYaxis()->SetTitle("Y");
  // tempPad = (TPad*) gPad;
  // DrawShashlikModule(tempPad);
  
 
  
    
  cc_Cal->Divide(2,2);
  cc_Cal->cd(1)->SetGrid();
  hHS_Cal_front->Draw("colz");
  hHS_Cal_front->GetXaxis()->SetTitle("cal X");
  hHS_Cal_front->GetYaxis()->SetTitle("cal Y");
  tempPad = (TPad*) gPad;
  DrawShashlikModule(tempPad);
  
  cc_Cal->cd(2)->SetGrid();
  hHS_Cal_back->Draw("colz");
  hHS_Cal_back->GetXaxis()->SetTitle("cal X");
  hHS_Cal_back->GetYaxis()->SetTitle("cal Y");
  tempPad = (TPad*) gPad;
  DrawShashlikModule(tempPad);
  
  cc_Cal->cd(3);
  hHS_Cal_front_module->Draw("colz");
  hHS_Cal_front_module->GetXaxis()->SetTitle("cal X");
  hHS_Cal_front_module->GetYaxis()->SetTitle("cal Y");
  tempPad = (TPad*) gPad;
  DrawShashlikModule(tempPad);
  
  cc_Cal->cd(4);
  hHS_Cal_back_module->Draw("colz");
  hHS_Cal_back_module->GetXaxis()->SetTitle("cal X");
  hHS_Cal_back_module->GetYaxis()->SetTitle("cal Y");
  tempPad = (TPad*) gPad;
  DrawShashlikModule(tempPad);
 
    
  TF1* fgaus = new TF1("fgaus","gaus(0)+pol2(3)",-30,30);
  fgaus->SetParameter(0,10);
  fgaus->SetParameter(1,0.0);
  fgaus->SetParLimits(1,-10,10);
  fgaus->SetParameter(2,0.5);
  fgaus->SetParLimits(2,0.1,3.0);
  
  
  fgaus->SetParameter(4,0.0);
  fgaus->SetParameter(5,0.0);
  
  std::cout << " =================================== " << std::endl;
  std::cout << " >>> X <<<" << std::endl; 
  
  cc_DX->Divide(4,2);
  cc_DX->cd(1)->SetGrid();
  X_hHS_HS1_Cal_front->Draw("colz");
  X_hHS_HS1_Cal_front->GetXaxis()->SetTitle("calo - hodoscope");
  X_hHS_HS1_Cal_front->GetYaxis()->SetTitle("event number");

  cc_DX->cd(2)->SetGrid();
  X_hHS_HS2_Cal_front->Draw("colz");
  X_hHS_HS2_Cal_front->GetXaxis()->SetTitle("calo - hodoscope");
  X_hHS_HS2_Cal_front->GetYaxis()->SetTitle("event number");

  cc_DX->cd(5)->SetGrid();
  X_hHS_HS1_Cal_back->Draw("colz");
  X_hHS_HS1_Cal_back->GetXaxis()->SetTitle("calo - hodoscope");
  X_hHS_HS1_Cal_back->GetYaxis()->SetTitle("event number");

  cc_DX->cd(6)->SetGrid();
  X_hHS_HS2_Cal_back->Draw("colz");
  X_hHS_HS2_Cal_back->GetXaxis()->SetTitle("calo - hodoscope");
  X_hHS_HS2_Cal_back->GetYaxis()->SetTitle("event number");
  

  std::cout<<"asdasd"<<endl;  
  cc_DX->cd(3)->SetGrid();
  X_h1_HS1_Cal_front->Draw();
  X_h1_HS1_Cal_front->GetXaxis()->SetTitle("calo - hodoscope");
  X_h1_HS1_Cal_front->Fit("fgaus","R");
 
 //  X_h1_HS1_Cal_front->GetFunction("gaus")->SetLineColor(3);


  cc_DX->cd(4)->SetGrid();
  X_h1_HS2_Cal_front->Draw();
  X_h1_HS2_Cal_front->GetXaxis()->SetTitle("calo - hodoscope");
  X_h1_HS2_Cal_front->Fit("fgaus","R");
  //X_h1_HS2_Cal_front->GetFunction("gaus")->SetLineColor(3);

  cc_DX->cd(7)->SetGrid();
  X_h1_HS1_Cal_back->Draw();
  X_h1_HS1_Cal_back->GetXaxis()->SetTitle("calo - hodoscope");
  X_h1_HS1_Cal_back->Fit("fgaus","R");
  //X_h1_HS1_Cal_back->GetFunction("gaus")->SetLineColor(3);
  
  cc_DX->cd(8)->SetGrid();
  X_h1_HS2_Cal_back->Draw();
  X_h1_HS2_Cal_back->GetXaxis()->SetTitle("calo - hodoscope");
  X_h1_HS2_Cal_back->Fit("fgaus","R");
  //X_h1_HS1_Cal_back->GetFunction("gaus")->SetLineColor(3);




  cc_sd->Divide(2,2);
  
  cc_sd->cd(1)->SetGrid();
  sd_h1x->Draw("colz");
  sd_h1x->GetXaxis()->SetTitle("no. of hits");
  sd_h1x->GetYaxis()->SetTitle("rms");
 
  
  cc_sd->cd(2)->SetGrid();
  sd_h2x->Draw("colz");
  sd_h2x->GetXaxis()->SetTitle("no. of hits");
  sd_h2x->GetYaxis()->SetTitle("rms");
 
  cc_sd->cd(3)->SetGrid();
  sd_h1y->Draw("colz");
  sd_h1y->GetXaxis()->SetTitle("no. of hits");
  sd_h1y->GetYaxis()->SetTitle("rms");
 
  
  cc_sd->cd(4)->SetGrid();
  sd_h2y->Draw("colz");
  sd_h2y->GetXaxis()->SetTitle("no. of hits");
  sd_h2y->GetYaxis()->SetTitle("rms");
  

  
 

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
    

 
  std::cout << " =================================== " << std::endl;
  std::cout << " >>> Y <<<" << std::endl; 
  
  cc_DY->Divide(4,2);
  cc_DY->cd(1)->SetGrid();
  Y_hHS_HS1_Cal_front->Draw("colz");
  Y_hHS_HS1_Cal_front->GetXaxis()->SetTitle("calo - hodoscope");
  Y_hHS_HS1_Cal_front->GetYaxis()->SetTitle("event number");
  
  cc_DY->cd(2)->SetGrid();
  Y_hHS_HS2_Cal_front->Draw("colz");
  Y_hHS_HS2_Cal_front->GetXaxis()->SetTitle("calo - hodoscope");
  Y_hHS_HS2_Cal_front->GetYaxis()->SetTitle("event number");
  
  cc_DY->cd(5)->SetGrid();
  Y_hHS_HS1_Cal_back->Draw("colz");
  Y_hHS_HS1_Cal_back->GetXaxis()->SetTitle("calo - hodoscope");
  Y_hHS_HS1_Cal_back->GetYaxis()->SetTitle("event number");
  
  cc_DY->cd(6)->SetGrid();
  Y_hHS_HS2_Cal_back->Draw("colz");
  Y_hHS_HS2_Cal_back->GetXaxis()->SetTitle("calo - hodoscope");
  Y_hHS_HS2_Cal_back->GetYaxis()->SetTitle("event number");
  
  
  cc_DY->cd(3)->SetGrid();
  Y_h1_HS1_Cal_front->Draw();
  Y_h1_HS1_Cal_front->GetXaxis()->SetTitle("calo - hodoscope");
  Y_h1_HS1_Cal_front->Fit("fgaus","R");
  //  Y_h1_HS1_Cal_front->GetFunction("gaus")->SetLineColor(3);

  cc_DY->cd(4)->SetGrid();
  Y_h1_HS2_Cal_front->Draw();
  Y_h1_HS2_Cal_front->GetXaxis()->SetTitle("calo - hodoscope");
  Y_h1_HS2_Cal_front->Fit("fgaus","R");
  // Y_h1_HS2_Cal_front->GetFunction("gaus")->SetLineColor(3);  

  cc_DY->cd(7)->SetGrid();
  Y_h1_HS1_Cal_back->Draw();
  Y_h1_HS1_Cal_back->GetXaxis()->SetTitle("calo - hodoscope");
  Y_h1_HS1_Cal_back->Fit("fgaus","R"); 
  //  Y_h1_HS1_Cal_back->GetFunction("gaus")->SetLineColor(3);
  
  cc_DY->cd(8)->SetGrid();
  Y_h1_HS2_Cal_back->Draw();
  Y_h1_HS2_Cal_back->GetXaxis()->SetTitle("calo - hodoscope");
  Y_h1_HS2_Cal_back->Fit("fgaus","R");
  //Y_h1_HS2_Cal_back->GetFunction("gaus")->SetLineColor(3);  
  
  
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
  
 
/* 
  //---- energy distribution
  cc_Energy_front->Divide(2,2);
  cc_Energy_front->cd(1)->SetGrid();  
  h_energy_front->Draw();
  h_energy_front->GetXaxis()->SetTitle("-log(E_{i}/E_{tot})");
  cc_Energy_front->cd(2)->SetGrid();  
  h_energy_front_size->Draw("colz");
  h_energy_front_size->GetXaxis()->SetTitle("-log(E_{i}/E_{tot})");
  h_energy_front_size->GetYaxis()->SetTitle("cluster size");
  cc_Energy_front->cd(3)->SetGrid();  
  Energy_Cal_front->Draw();
  Energy_Beam->Draw("same");
  Energy_Cal_front->GetXaxis()->SetTitle("cluster energy [GeV]");
  
  
  cc_Energy_back->Divide(2,2);
  cc_Energy_back->cd(1)->SetGrid();  
  h_energy_back->Draw();
  h_energy_back->GetXaxis()->SetTitle("-log(E_{i}/E_{tot})");
  cc_Energy_back->cd(2)->SetGrid();  
  h_energy_back_size->Draw("colz");
  h_energy_back_size->GetXaxis()->SetTitle("-log(E_{i}/E_{tot})");
  h_energy_back_size->GetYaxis()->SetTitle("cluster size");
  cc_Energy_back->cd(3)->SetGrid();  
  Energy_Cal_back->Draw();
  Energy_Beam->Draw("same");
  Energy_Cal_back->GetXaxis()->SetTitle("cluster energy [GeV]");
  
*/
  
  gMyRootApp->Run(); 

}

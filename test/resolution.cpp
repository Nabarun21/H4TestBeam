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
#include "TStyle.h"

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

//TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");



//---- draw shashlik matrix

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

  gStyle->SetOptFit(1);
 
  std::string input_file;
  int maxEvents = -1;
  int doFiber = 0;
  float table_x_reference = 200; //---- mm
  float table_y_reference = 350; //---- mm
  float table_x = 200; //---- mm
  float table_y = 350; //---- mm
 
  float w0 = 4.5;
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
  

  
  TCanvas* cc_hodo_scat = new TCanvas ("cc_hodo_scat","",800,600);
  TCanvas* cc_hodo_scat_y = new TCanvas ("cc_hodo_scat_y","",800,600);

  TCanvas* cc_diff_scat = new TCanvas ("cc_diff_scat","",800,600);
  TCanvas* cc_diff_scat_y = new TCanvas ("cc_diff_scat_y","",800,600);


  
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
  

  
  TH2F *X_hHS_HS1_Cal_front  = new TH2F("X_hHS_HS1_Cal_front", "X Hodoscope 1 vs Cal front ", 64, -32, 32, nEntries, 0, nEntries);
  TH2F *X_hHS_HS2_Cal_front  = new TH2F("X_hHS_HS2_Cal_front", "X Hodoscope 2 vs Cal front ", 64, -32, 32, nEntries, 0, nEntries);
  TH2F *X_hHS_HS1_Cal_back   = new TH2F("X_hHS_HS1_Cal_back",  "X Hodoscope 1 vs Cal back " , 64, -32, 32, nEntries, 0, nEntries);
  TH2F *X_hHS_HS2_Cal_back   = new TH2F("X_hHS_HS2_Cal_back",  "X Hodoscope 2 vs Cal back " , 64, -32, 32, nEntries, 0, nEntries);

  
  TH1F *X_h1_HS1_Cal_front  = new TH1F("X_h1_HS1_Cal_front", "X Hodoscope 1 vs Cal front ", 256, -32, 32);
  TH1F *X_h1_HS2_Cal_front  = new TH1F("X_h1_HS2_Cal_front", "X Hodoscope 2 vs Cal front ", 256, -32, 32);
  TH1F *X_h1_HS1_Cal_back   = new TH1F("X_h1_HS1_Cal_back",  "X Hodoscope 1 vs Cal back " , 256, -32, 32);
  TH1F *X_h1_HS2_Cal_back   = new TH1F("X_h1_HS2_Cal_back",  "X Hodoscope 2 vs Cal back " , 256, -32, 32);

  


  TCanvas* cc_DY = new TCanvas ("cc_DY","",800,800);  TCanvas* cc_DX = new TCanvas ("cc_DX","",800,800);
  
  TH2F *Y_hHS_HS1_Cal_front  = new TH2F("Y_hHS_HS1_Cal_front", "Y Hodoscope 1 vs Cal front ", 64, -32, 32, nEntries, 0, nEntries);
  TH2F *Y_hHS_HS2_Cal_front  = new TH2F("Y_hHS_HS2_Cal_front", "Y Hodoscope 2 vs Cal front ", 64, -32, 32, nEntries, 0, nEntries);
  TH2F *Y_hHS_HS1_Cal_back   = new TH2F("Y_hHS_HS1_Cal_back",  "Y Hodoscope 1 vs Cal back " , 64, -32, 32, nEntries, 0, nEntries);
  TH2F *Y_hHS_HS2_Cal_back   = new TH2F("Y_hHS_HS2_Cal_back",  "Y Hodoscope 2 vs Cal back " , 64, -32, 32, nEntries, 0, nEntries);
  
  
  TH1F *Y_h1_HS1_Cal_front  = new TH1F("Y_h1_HS1_Cal_front", "Y Hodoscope 1 vs Cal front ", 256, -32, 32);
  TH1F *Y_h1_HS2_Cal_front  = new TH1F("Y_h1_HS2_Cal_front", "Y Hodoscope 2 vs Cal front ", 256, -32, 32);
  TH1F *Y_h1_HS1_Cal_back   = new TH1F("Y_h1_HS1_Cal_back",  "Y Hodoscope 1 vs Cal back " , 256, -32, 32);
  TH1F *Y_h1_HS2_Cal_back   = new TH1F("Y_h1_HS2_Cal_back",  "Y Hodoscope 2 vs Cal back " , 256, -32, 32);
  
  
  
  bool haverechits = false;
  std::vector<TBRecHit> *rechits=0;
  if(H4tree->GetListOfBranches()->FindObject("tbrechits")) {
    std::cout << " found rechits " << std::endl;
    H4tree->SetBranchAddress("tbrechits",&rechits);
    haverechits = true;
  }
  
  //   Mapper *mapper = Mapper::Instance();
  
  std::vector<float> h1x;
  std::vector<float> h1y;
  std::vector<float> h2x;
  std::vector<float> h2y;

  std::vector<float> calox;
  std::vector<float> caloy;

  std::vector<float> calox_bk;
  std::vector<float> caloy_bk;

  std::vector<float> l1x;
  std::vector<float> l1y;
  std::vector<float> l2x;
  std::vector<float> l2y;
  std::vector<float> m1x;
  std::vector<float> m1y;
  std::vector<float> m2x;
  std::vector<float> m2y;



  
  CaloCluster* caloCluster = new CaloCluster();
  caloCluster->setW0(4.95);
  caloCluster->setInterCalibrationConstants("data/InterCalibrationConstants3.txt");
    
  
  CaloCluster* caloCluster2 = new CaloCluster();
  caloCluster2->setW0(5);
  caloCluster2->setInterCalibrationConstants("data/InterCalibrationConstants3.txt");
    
  
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
   
    //---- calorimeter data
  
    if (i==0) caloCluster->setMapperEpoch(tbevent->GetTimeStamp());   
     
   
    std::vector<float> caloCluster_position_X_front;
    caloCluster->doCalorimeterReconstruction( rechits, 1, 60, doFiber);
    caloCluster_position_X_front.push_back( caloCluster->getPositionX() );

    caloCluster->doCalorimeterReconstruction( rechits, -1, 60, doFiber);
    std::vector<float> caloCluster_position_X_back;   
    caloCluster_position_X_back.push_back( caloCluster->getPositionX() );
   
    //    delete caloCluster;



    std::vector<float> caloCluster_position_Y_front;
    caloCluster2->doCalorimeterReconstruction( rechits, 1, 60, doFiber);
    caloCluster_position_Y_front.push_back( caloCluster2->getPositionY() );

    std::vector<float> caloCluster_position_Y_back;   
    caloCluster2->doCalorimeterReconstruction( rechits, -1, 60, doFiber);
    caloCluster_position_Y_back.push_back( caloCluster2->getPositionY() );
   
    // delete caloCluster2;
   
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
    doHodoReconstruction( fibers_X1, n_fibers_X1, pos_fibers_X1,0); //(table_x - 200)); //---- change of coordinate system using numbers from googledoc
    doHodoReconstruction( fibers_X2, n_fibers_X2, pos_fibers_X2,0); //(table_x - 200)); //---- change of coordinate system using numbers from googledoc
    doHodoReconstruction( fibers_Y1, n_fibers_Y1, pos_fibers_Y1,0);// (table_y - 350)); //---- change of coordinate system using numbers from googledoc
    doHodoReconstruction( fibers_Y2, n_fibers_Y2, pos_fibers_Y2,0);// (table_y - 350)); //---- change of coordinate system using numbers from googledoc



     
  
   
    //---- now merge and compare
    if (pos_fibers_X1.size() ==1 && pos_fibers_X2.size() ==1)

      {     
	double sd1=0;
	double ave1=0;
	mean_stdev(sd1,ave1,pos_fibers_X1);
     
	double sd2=0;
	double ave2=0;
	mean_stdev(sd2,ave2,pos_fibers_X2);
     
	float c_ave1=mean(caloCluster_position_X_front);    
	float c_ave2=mean(caloCluster_position_X_back);    
     

	if( (ave1-ave2) >=-3.92 && (ave1-ave2) <=-1.92  )  {
	  //---- X

	  h1x.push_back(ave1);
	  h2x.push_back(ave2);
     
	  calox.push_back(c_ave1/0.72);
	  calox_bk.push_back(c_ave2/0.77);


	  l1x.push_back(ave1-c_ave1/0.72);
	  l2x.push_back(ave2-c_ave1/0.72);
	  m1x.push_back(ave1-c_ave2/0.77);
	  m2x.push_back(ave2-c_ave2/0.77);

   
	  for (int iCalo = 0; iCalo < caloCluster_position_X_front.size(); iCalo++) {
          
	    hHS_HS1_Cal_front_X->Fill(caloCluster_position_X_front.at(iCalo)/0.72, ave1);
	    X_hHS_HS1_Cal_front->Fill(caloCluster_position_X_front.at(iCalo)/0.72 - ave1, i);
	    X_h1_HS1_Cal_front->Fill(caloCluster_position_X_front.at(iCalo)/0.72 - ave1);
       
	  }   
	  for (int iCalo = 0; iCalo < caloCluster_position_X_back.size(); iCalo++) {
             
	    hHS_HS1_Cal_back_X->Fill(caloCluster_position_X_back.at(iCalo)/0.77, ave1);
	    X_hHS_HS1_Cal_back->Fill(caloCluster_position_X_back.at(iCalo)/0.77 - ave1, i);
	    X_h1_HS1_Cal_back->Fill(caloCluster_position_X_back.at(iCalo)/0.77 - ave1);

     
	  }
       
 
	  for (int iCalo = 0; iCalo < caloCluster_position_X_front.size(); iCalo++) {
          
	    hHS_HS2_Cal_front_X->Fill(caloCluster_position_X_front.at(iCalo)/0.72, ave2);
	    X_hHS_HS2_Cal_front->Fill(caloCluster_position_X_front.at(iCalo)/0.72 - ave2, i);
	    X_h1_HS2_Cal_front->Fill(caloCluster_position_X_front.at(iCalo)/0.72 - ave2);
       
	  }   
	  for (int iCalo = 0; iCalo < caloCluster_position_X_back.size(); iCalo++) {
             
	    hHS_HS2_Cal_back_X->Fill(caloCluster_position_X_back.at(iCalo)/0.77, ave2);
	    X_hHS_HS2_Cal_back->Fill(caloCluster_position_X_back.at(iCalo)/0.77 - ave2, i);
	    X_h1_HS2_Cal_back->Fill(caloCluster_position_X_back.at(iCalo)/0.77 - ave2);

     
	  }
       
	}
      }
   
   
   
   
    if (pos_fibers_Y1.size() ==1 && pos_fibers_Y2.size() ==1) {
     
     
      double sd1=0;
      double ave1=0;
      mean_stdev(sd1,ave1,pos_fibers_Y1);
     
     
      double sd2=0;    
      double ave2=0;
      mean_stdev(sd2,ave2,pos_fibers_Y2);

      float c_ave1=mean(caloCluster_position_Y_front);    
      float c_ave2=mean(caloCluster_position_Y_back);    
     

      if( (ave1-ave2) >=-0.6035 && (ave1-ave2) <=1.3965  )  {	

	h1y.push_back(ave1);
	h2y.push_back(ave2);
     
	caloy.push_back(c_ave1/0.59);
	caloy_bk.push_back(c_ave2/0.85);

	l1y.push_back(ave1-c_ave1/0.59);
	l2y.push_back(ave2-c_ave1/0.59);
	m1y.push_back(ave1-c_ave2/0.85);
	m2y.push_back(ave2-c_ave2/0.85);

	for (int iCalo = 0; iCalo < caloCluster_position_Y_front.size(); iCalo++) {
	 
	  hHS_HS1_Cal_front_Y->Fill(caloCluster_position_Y_front.at(iCalo)/0.59, ave1);
	  Y_hHS_HS1_Cal_front->Fill(caloCluster_position_Y_front.at(iCalo)/0.59 - ave1, i);
	  Y_h1_HS1_Cal_front->Fill(caloCluster_position_Y_front.at(iCalo)/0.59 - ave1);
	 
	 
	}
	for (int iCalo = 0; iCalo < caloCluster_position_Y_back.size(); iCalo++) {
			 
	  hHS_HS1_Cal_back_Y->Fill(caloCluster_position_Y_back.at(iCalo)/0.85, ave1);
	  Y_hHS_HS1_Cal_back->Fill(caloCluster_position_Y_back.at(iCalo)/0.85 - ave1, i);
	  Y_h1_HS1_Cal_back->Fill(caloCluster_position_Y_back.at(iCalo)/0.85 - ave1);
	 
	}
   

	for (int iCalo = 0; iCalo < caloCluster_position_Y_front.size(); iCalo++) {
	 
	  hHS_HS2_Cal_front_Y->Fill(caloCluster_position_Y_front.at(iCalo)/0.59, ave2);
	  Y_hHS_HS2_Cal_front->Fill(caloCluster_position_Y_front.at(iCalo)/0.59 - ave2, i);
	  Y_h1_HS2_Cal_front->Fill(caloCluster_position_Y_front.at(iCalo)/0.59 - ave2);
	 
	 
	}
	for (int iCalo = 0; iCalo < caloCluster_position_Y_back.size(); iCalo++) {
			 
	  hHS_HS2_Cal_back_Y->Fill(caloCluster_position_Y_back.at(iCalo)/0.85, ave2);
	  Y_hHS_HS2_Cal_back->Fill(caloCluster_position_Y_back.at(iCalo)/0.85 - ave2, i);
	  Y_h1_HS2_Cal_back->Fill(caloCluster_position_Y_back.at(iCalo)/0.85 - ave2);
	 
	}
      }
    }
  }
   
  

  float *ah1x=&(h1x[0]);
  float *ah2x=&(h2x[0]);

  float *ah1y=&(h1y[0]);
  float *ah2y=&(h2y[0]);


  float *acalox=&(calox[0]);
  float *acaloy=&(caloy[0]);

  float *acaloy_bk=&(caloy_bk[0]);
  float *acalox_bk=&(calox_bk[0]);

  float *al1x=&(l1x[0]);
  float *al2x=&(l2x[0]);

  float *al1y=&(l1y[0]);
  float *al2y=&(l2y[0]);

  float *am1x=&(m1x[0]);
  float *am2x=&(m2x[0]);

  float *am1y=&(m1y[0]);
  float *am2y=&(m2y[0]);

  



  TGraph* Tg=new TGraph(h1x.size(),ah1x,acalox);

  TGraph* Tg2=new TGraph(h2x.size(),ah2x,acalox);
  
  TGraph* Tg3=new TGraph(h1x.size(),ah1x,acalox_bk);

  TGraph* Tg4=new TGraph(h2x.size(),ah2x,acalox_bk);
  
  


  TGraph* Tg5=new TGraph(h1y.size(),ah1y,acaloy);
  
  TGraph* Tg6=new TGraph(h2y.size(),ah2y,acaloy);

  TGraph* Tg7=new TGraph(h1y.size(),ah1y,acaloy_bk);
  
  TGraph* Tg8=new TGraph(h2y.size(),ah2y,acaloy_bk);
  

  TGraph* Mg=new TGraph(l1x.size(),ah1x,al1x);

  TGraph* Mg2=new TGraph(l2x.size(),ah2x,al2x);
  
  TGraph* Mg3=new TGraph(l1y.size(),ah1y,al1y);

  TGraph* Mg4=new TGraph(l2y.size(),ah2y,al2y);


  TGraph* Mg5=new TGraph(m1x.size(),ah1x,am1x);

  TGraph* Mg6=new TGraph(m2x.size(),ah2x,am2x);
  
  TGraph* Mg7=new TGraph(m1y.size(),ah1y,am1y);

  TGraph* Mg8=new TGraph(m2y.size(),ah2y,am2y);


  


  TF1* lin = new TF1("lin","[0]*x+[1]",-30,30);
  lin->SetParameter(0,1);
  lin->SetParameter(1,0);

  float mean_xf=X_h1_HS1_Cal_front->GetMean();
  float rms_xf=X_h1_HS1_Cal_front->GetRMS();
  float mean_yf=Y_h1_HS1_Cal_front->GetMean();
  float rms_yf=Y_h1_HS1_Cal_front->GetRMS();

  float mean_xb=X_h1_HS1_Cal_back->GetMean();
  float rms_xb=X_h1_HS1_Cal_back->GetRMS();
  float mean_yb=Y_h1_HS1_Cal_back->GetMean();
  float rms_yb=Y_h1_HS1_Cal_back->GetRMS();



  TF1* fgaus_xf = new TF1("fgaus_xf","gaus(0)+pol2(3)",mean_xf-3*rms_xf,mean_xf+3*rms_xf);
  fgaus_xf->SetParameter(0,10);
  fgaus_xf->SetParameter(1,0.0);
  fgaus_xf->SetParLimits(1,-10,10);
  fgaus_xf->SetParameter(2,0.5);
  fgaus_xf->SetParLimits(2,0.1,3.0);
  
  fgaus_xf->SetParameter(4,0.0);
  fgaus_xf->SetParameter(5,0.0);
 

  TF1* fgaus_xb = new TF1("fgaus_xb","gaus(0)+pol2(3)",mean_xb-3*rms_xb,mean_xb+3*rms_xb);
  fgaus_xb->SetParameter(0,10);
  fgaus_xb->SetParameter(1,0.0);
  fgaus_xb->SetParLimits(1,-10,10);
  fgaus_xb->SetParameter(2,0.5);
  fgaus_xb->SetParLimits(2,0.1,3.0);

  fgaus_xb->SetParameter(4,0.0);
  fgaus_xb->SetParameter(5,0.0);
 

  

  TF1* fgaus_yf = new TF1("fgaus_yf","gaus(0)+pol2(3)",mean_yf-2*rms_yf,mean_yf+2*rms_yf);
  fgaus_yf->SetParameter(0,10);
  fgaus_yf->SetParameter(1,0.0);
  fgaus_yf->SetParLimits(1,-10,10);
  fgaus_yf->SetParameter(2,0.5);
  fgaus_yf->SetParLimits(2,0.1,3.0);

  fgaus_yf->SetParameter(4,0.0);
  fgaus_yf->SetParameter(5,0.0);
 

  TF1* fgaus_yb = new TF1("fgaus_yb","gaus(0)+pol2(3)",mean_yb-3*rms_yb,mean_yb+3*rms_yb);
  fgaus_yb->SetParameter(0,10);
  fgaus_yb->SetParameter(1,0.0);
  fgaus_yb->SetParLimits(1,-10,10);
  fgaus_yb->SetParameter(2,0.5);
  fgaus_yb->SetParLimits(2,0.1,3.0);


  
  fgaus_yb->SetParameter(4,0.0);
  fgaus_yb->SetParameter(5,0.0);
  
  TF1* fxy = new TF1 ("fxy","x",-20,20);

  
  //scatter plots


  cc_hodo_scat->Divide(2,2);

  cc_hodo_scat->cd(1)->SetGrid();
  Tg->SetMarkerColor(5);
  Tg->Draw("ap");
  Tg->GetXaxis()->SetTitle("hodo1_x");  
  Tg->GetYaxis()->SetTitle("calo_x");  
  Tg->SetTitle("HODO1X VS CALOX");
  Tg->Fit("lin","R");

  cc_hodo_scat->cd(2)->SetGrid();
  Tg->SetMarkerStyle(6);
  Tg2->SetMarkerColor(6);
  Tg2->Draw("ap");
  Tg2->GetXaxis()->SetTitle("hodo2_x");  
  Tg2->GetYaxis()->SetTitle("calo_x");  
  Tg2->SetTitle("HODO2X VS CALOX");
  Tg2->Fit("lin","R");

  cc_hodo_scat->cd(3)->SetGrid();
  Tg3->Draw("ap");
  Tg3->GetXaxis()->SetTitle("hodo1_x");  
  Tg3->GetYaxis()->SetTitle("calo_x_back");  
  Tg3->SetTitle("HODO1X VS CALOX_BACK");
  Tg3->Fit("lin","R");
   
  cc_hodo_scat->cd(4)->SetGrid();
  Tg4->Draw("ap");
  Tg4->GetXaxis()->SetTitle("hodo2_x");  
  Tg4->GetYaxis()->SetTitle("calo_x_back");  
  Tg4->SetTitle("HODO2X VS CALOX_BACK");
  Tg4->Fit("lin","R");

  cc_hodo_scat->SaveAs("plots_evan_calib/scat_x.pdf");

  cc_hodo_scat_y->Divide(2,2);

  cc_hodo_scat_y->cd(1)->SetGrid();
  Tg5->Draw("ap");
  Tg5->GetXaxis()->SetTitle("hodo1_y");  
  Tg5->GetYaxis()->SetTitle("calo_y");  
  Tg5->SetTitle("HODO1Y VS CALOY");
  Tg5->Fit("lin","R");
   
  cc_hodo_scat_y->cd(2)->SetGrid();
  Tg6->Draw("ap");
  Tg6->GetXaxis()->SetTitle("hodo2_y");  
  Tg6->GetYaxis()->SetTitle("calo_y");  
  Tg6->SetTitle("HODO2Y VS CALOY");
  Tg6->Fit("lin","R");

  cc_hodo_scat_y->cd(3)->SetGrid();
  Tg7->Draw("ap");
  Tg7->GetXaxis()->SetTitle("hodo1_y");  
  Tg7->GetYaxis()->SetTitle("calo_y_back");  
  Tg7->SetTitle("HODO1Y VS CALOY_BACK");
  Tg7->Fit("lin","R");
   
  cc_hodo_scat_y->cd(4)->SetGrid();
  Tg8->Draw("ap");
  Tg8->GetXaxis()->SetTitle("hodo2_y");  
  Tg8->GetYaxis()->SetTitle("calo_y_back");  
  Tg8->SetTitle("HODO2Y VS CALOY_BACK");
  Tg8->Fit("lin","R");


  cc_hodo_scat_y->SaveAs("plots_evan_calib/Scat_y.pdf");
  //diff vs hodox

  
  cc_diff_scat->Divide(2,2);

  cc_diff_scat->cd(1)->SetGrid();
  Mg->SetMarkerColor(5);
  Mg->Draw("ap");
  Mg->GetXaxis()->SetTitle("hodo1_x");  
  Mg->GetYaxis()->SetTitle("calo-hodo_x");  
  // Mg->SetTitle("HODO1X VS CALO-HODOX");
  Mg->Fit("lin","R");
   
  cc_diff_scat->cd(2)->SetGrid();
  Mg->SetMarkerStyle(6);
  Mg2->SetMarkerColor(6);
  Mg2->Draw("ap");
  Mg2->GetXaxis()->SetTitle("hodo2_x");  
  Mg2->GetYaxis()->SetTitle("calo-hodo_x");  
  // Mg2->SetTitle("HODO2X VS CALO-HODOX");
  Mg2->Fit("lin","R");

  cc_diff_scat->cd(3)->SetGrid();
  Mg3->Draw("ap");
  Mg3->GetXaxis()->SetTitle("hodo1_x");  
  Mg3->GetYaxis()->SetTitle("calo-hodo_x_back");  
  //Mg3->SetTitle("HODO1X VS CALO-HODOX_BACK");
  Mg3->Fit("lin","R");
   
  cc_diff_scat->cd(4)->SetGrid();
  Mg4->Draw("ap");
  Mg4->GetXaxis()->SetTitle("hodo2_x");  
  Mg4->GetYaxis()->SetTitle("calo-hodo_x_back");  
  // Mg4->SetTitle("HODO2X VS CALO-HODOX_BACK");
  Mg4->Fit("lin","R");


  cc_diff_scat->SaveAs("plots_evan_calib/diff_x.pdf");

  cc_diff_scat_y->Divide(2,2);

  cc_diff_scat_y->cd(1)->SetGrid();
  Mg5->Draw("ap");
  Mg5->GetXaxis()->SetTitle("hodo1_y");  
  Mg5->GetYaxis()->SetTitle("calo-hodo_y");  
  // Mg5->SetTitle("HODO1Y VS CALO-HODOY");
  Mg5->Fit("lin","R");
   
  cc_diff_scat_y->cd(2)->SetGrid();
  Mg6->Draw("ap");
  Mg6->GetXaxis()->SetTitle("hodo2_y");  
  Mg6->GetYaxis()->SetTitle("calo-hodo_y");  
  //Mg6->SetTitle("HODO2Y VS CALO-HODOY");
  Mg6->Fit("lin","R");

  cc_diff_scat_y->cd(3)->SetGrid();
  Mg7->Draw("ap");
  Mg7->GetXaxis()->SetTitle("hodo1_y");  
  Mg7->GetYaxis()->SetTitle("calo-hodo_y_back");  
  //Mg7->SetTitle("HODO1Y VS CALO-HODOY_BACK");
  Mg7->Fit("lin","R");
   
  cc_diff_scat_y->cd(4)->SetGrid();
  Mg8->Draw("ap");
  Mg8->GetXaxis()->SetTitle("hodo2_y");  
  Mg8->GetYaxis()->SetTitle("calo-hodo_y_back");  
  //Mg8->SetTitle("HODO2Y VS CALO-HODOY_BACK");
  Mg8->Fit("lin","R");

  cc_diff_scat_y->SaveAs("plots_evan_calib/diff_y.pdf");
  //-----------


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
    
  cc_X->SaveAs("plots_evan_calib/calx.pdf");
  
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
  

  cc_DX->cd(3)->SetGrid();
  X_h1_HS1_Cal_front->Draw();
  X_h1_HS1_Cal_front->GetXaxis()->SetTitle("calo - hodoscope");
  X_h1_HS1_Cal_front->Fit("fgaus_xf","R");
 
  //  X_h1_HS1_Cal_front->GetFunction("gaus")->SetLineColor(3);


  cc_DX->cd(4)->SetGrid();
  X_h1_HS2_Cal_front->Draw();
  X_h1_HS2_Cal_front->GetXaxis()->SetTitle("calo - hodoscope");
  X_h1_HS2_Cal_front->Fit("fgaus_xf","R");
  //X_h1_HS2_Cal_front->GetFunction("gaus")->SetLineColor(3);

  cc_DX->cd(7)->SetGrid();
  X_h1_HS1_Cal_back->Draw();
  X_h1_HS1_Cal_back->GetXaxis()->SetTitle("calo - hodoscope");
  X_h1_HS1_Cal_back->Fit("fgaus_xb","R");
  //X_h1_HS1_Cal_back->GetFunction("gaus")->SetLineColor(3);
  
  cc_DX->cd(8)->SetGrid();
  X_h1_HS2_Cal_back->Draw();
  X_h1_HS2_Cal_back->GetXaxis()->SetTitle("calo - hodoscope");
  X_h1_HS2_Cal_back->Fit("fgaus_xb","R");
  //X_h1_HS1_Cal_back->GetFunction("gaus")->SetLineColor(3);
  
  cc_DX->SaveAs("plots_evan_calib/x.pdf");
  
  
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

  cc_Y->SaveAs("plots_evan_calib/caly.pdf");


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
  Y_h1_HS1_Cal_front->Fit("fgaus_yf","R");
  //  Y_h1_HS1_Cal_front->GetFunction("gaus")->SetLineColor(3);

  cc_DY->cd(4)->SetGrid();
  Y_h1_HS2_Cal_front->Draw();
  Y_h1_HS2_Cal_front->GetXaxis()->SetTitle("calo - hodoscope");
  Y_h1_HS2_Cal_front->Fit("fgaus_yf","R");
  // Y_h1_HS2_Cal_front->GetFunction("gaus")->SetLineColor(3);  

  cc_DY->cd(7)->SetGrid();
  Y_h1_HS1_Cal_back->Draw();
  Y_h1_HS1_Cal_back->GetXaxis()->SetTitle("calo - hodoscope");
  Y_h1_HS1_Cal_back->Fit("fgaus_yb","R"); 
  //  Y_h1_HS1_Cal_back->GetFunction("gaus")->SetLineColor(3);
  
  cc_DY->cd(8)->SetGrid();
  Y_h1_HS2_Cal_back->Draw();
  Y_h1_HS2_Cal_back->GetXaxis()->SetTitle("calo - hodoscope");
  Y_h1_HS2_Cal_back->Fit("fgaus_yb","R");
  //Y_h1_HS2_Cal_back->GetFunction("gaus")->SetLineColor(3);  
  cc_DY->SaveAs("plots_evan_calib/y.pdf");
  
  gMyRootApp->Run(); 

}

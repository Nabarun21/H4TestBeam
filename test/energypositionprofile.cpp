#include <string>
#include "TStyle.h"
#include "TTree.h"
#include "TFile.h"
#include "TApplication.h"
#include "TChain.h"
#include "TTreeFormula.h"
#include "TCanvas.h"
#include "TProfile2D.h"
#include "TLine.h"

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
 void doHodoReconstruction( std::vector<int> input_values, std::vector<int>& nFibres, std::vector<float>& cluster_position, float shift) {
 std::vector<HodoCluster*> clusters = getHodoClusters( input_values );
 for( unsigned i=0; i<clusters.size(); ++i ) {
  nFibres.push_back( clusters[i]->getSize() );
  cluster_position.push_back( clusters[i]->getPosition() + shift );
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
 
 float w0 = 4.9;
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
  
  TH1F *Energy_Cal_back   = new TH1F("Energy_Cal_back",  "Energy Cal back " , 200, 0, 200);
  TH1F *Energy_Cal_front  = new TH1F("Energy_Cal_front", "Energy Cal front" , 200, 0, 200);

  TProfile2D * energyvxy_hodo1=new TProfile2D("energypositionprofile_hodo1","enery v xy:hodo1",64, -32, 32, 64, -32, 32,0,55);
  TProfile2D * energyvxy_hodo2=new TProfile2D("energypositionprofile_hodo2","enery v xy:hodo2",64, -32, 32, 64, -32, 32,0,55);
  


  
  //  TApplication* gMyRootApp = new TApplication("My ROOT Application", &argc, argv);
  
  
  

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
  for (int num_runs=0;num_runs<input_files_vector.size();num_runs++) 
    {
      // cout<<"chupa "<<num_runs<<endl;
      TString filename="../../../public/new_merged_files/calo_hodo_beam_"+input_files_vector.at(num_runs)+".root"; 
      TFile *f= new TFile(filename); 
      // cout<<"chupa "<<num_runs<<endl;
      // //The tree 
 
       TTree* H4tree = (TTree*)f->Get("t1041");
  
       //  TChain* H4tree = new TChain("t1041");
      
      cout<<"Running on run :"<<input_files_vector.at(num_runs)<<endl;
      //  for (unsigned int i=0; i<input_files_vector.size(); i++) {
      //H4tree->Add(input_files_vector.at(num_runs).c_str());
      // }
  
  
      //---- read file
      int nEntries = H4tree->GetEntries(); 
      std::cout << " nEntries = " << nEntries << std::endl;
      if (maxEvents != -1) nEntries = maxEvents>nEntries ? nEntries : maxEvents ;
      std::cout << " new nEntries = " << nEntries << std::endl;
  
  
      TBSpill* tbspill = new TBSpill();
      TBEvent* tbevent = new TBEvent();
  
      H4tree->SetBranchAddress("tbevent", &tbevent);
      H4tree->SetBranchAddress("tbspill", &tbspill);

  
  
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
      caloCluster->setInterCalibrationConstants("data/InterCalibrationConstants3.txt");

      //      cout<<"chupi1 "<<endl;

  
      ///---- loop ----
      for (int i=0; i<nEntries; i++) {
   
	if ((i%1000)==0) {
	  // std::cout <<  " entry: " << i << "::" << nEntries << std::endl;
	}
   	
	H4tree->GetEntry(i);
	

	float table_x_shift = tbspill->GetShiftX();
	float table_y_shift = tbspill->GetShiftY();
   
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
   
   
	//	Energy_Beam->Fill(tbspill->GetMomentum());

   	   
	//---- calorimeter data
	if (i==0) caloCluster->setMapperEpoch(tbevent->GetTimeStamp());
   
   
	// std::vector<float> caloCluster_position_X_front;
	// std::vector<float> caloCluster_position_Y_front;
	std::vector<float> caloCluster_Energy_front;
	//   std::vector<float> caloCluster_Energies_front;
   
	caloCluster->doCalorimeterReconstruction( rechits, 1, 30, doFiber);
   
	// caloCluster_position_X_front.push_back( caloCluster->getPositionX() );
	// caloCluster_position_Y_front.push_back( caloCluster->getPositionY() );
	caloCluster_Energy_front.push_back( caloCluster->getEnergy() );
	// caloCluster_Energies_front = caloCluster->getCaloClusterComponents();
   
	Energy_Cal_front->Fill(caloCluster_Energy_front.at(0));
   
	// std::vector<float> caloCluster_position_X_back;
	// std::vector<float> caloCluster_position_Y_back;
	std::vector<float> caloCluster_Energy_back;
	//   std::vector<float> caloCluster_Energies_back;;
   
     
	caloCluster->doCalorimeterReconstruction( rechits, -1, 30, doFiber);
   
	// caloCluster_position_X_back.push_back( caloCluster->getPositionX() );
	// caloCluster_position_Y_back.push_back( caloCluster->getPositionY() );
	caloCluster_Energy_back.push_back( caloCluster->getEnergy() );
	//   caloCluster_Energies_back = caloCluster->getCaloClusterComponents();
   
	Energy_Cal_back->Fill(caloCluster_Energy_back.at(0));
   
   
	// //---- modular level DR = 5 mm
	// std::vector<float> caloCluster_position_X_front_module;
	// std::vector<float> caloCluster_position_Y_front_module;
	// std::vector<float> caloCluster_Energy_front_module;
   
	// caloCluster->doCalorimeterReconstruction( rechits, 1, 5, doFiber);

	// caloCluster_position_X_front_module.push_back( caloCluster->getPositionX() );
	// caloCluster_position_Y_front_module.push_back( caloCluster->getPositionY() );
	// caloCluster_Energy_front_module.push_back( caloCluster->getEnergy() );
      
   
	// std::vector<float> caloCluster_position_X_back_module;
	// std::vector<float> caloCluster_position_Y_back_module;
	// std::vector<float> caloCluster_Energy_back_module;
   
	// caloCluster->doCalorimeterReconstruction( rechits, -1, 5, doFiber);
   
	// caloCluster_position_X_back_module.push_back( caloCluster->getPositionX() );
	// caloCluster_position_Y_back_module.push_back( caloCluster->getPositionY() );
	// caloCluster_Energy_back_module.push_back( caloCluster->getEnergy() );


   
   
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
   
   
	doHodoReconstruction( fibers_X1, n_fibers_X1, pos_fibers_X1, table_x_shift); //---- change of coordinate system using numbers from googledoc
	doHodoReconstruction( fibers_X2, n_fibers_X2, pos_fibers_X2, table_x_shift); //---- change of coordinate system using numbers from googledoc
	doHodoReconstruction( fibers_Y1, n_fibers_Y1, pos_fibers_Y1, table_y_shift); //---- change of coordinate system using numbers from googledoc
	doHodoReconstruction( fibers_Y2, n_fibers_Y2, pos_fibers_Y2, table_y_shift); //---- change of coordinate system using numbers from googledoc



	float av_energy=(caloCluster_Energy_front.at(0)+caloCluster_Energy_back.at(0))/200.;
	if (caloCluster_Energy_front.at(0)/100>5.)
	  {
	    //---- now merge and compare
	    if ((pos_fibers_X1.size() == 1) && (pos_fibers_Y1.size() == 1)) 
	      {
		energyvxy_hodo1->Fill(pos_fibers_X1.at(0),pos_fibers_Y1.at(0),caloCluster_Energy_front.at(0)/100);
	      }
	    
	    if ((pos_fibers_X2.size() == 1) && (pos_fibers_Y2.size() == 1)) 
	      {
		energyvxy_hodo2->Fill(pos_fibers_X2.at(0),pos_fibers_Y2.at(0),caloCluster_Energy_front.at(0)/100);
	      }
	  }
   
      }
      delete H4tree; 
      //  }
    }
  //gRoot
      //TStyle * myStyle=new TStyle();
  //  myStyle->cd();
  // myStyle->SetPalette(); 
  //  myStyle->SetOptStat(000000000); 
  gStyle->SetOptStat(0);
  TPad* tempPad;
  
  //---- plot ----

  TCanvas* cc_Cal = new TCanvas ("cc_Cal","",1200,1200);
  cc_Cal->cd();
  energyvxy_hodo1->Draw("colz");
  energyvxy_hodo1->GetXaxis()->SetTitle("X1");
  energyvxy_hodo1->GetYaxis()->SetTitle("Y1");
  tempPad = (TPad*) gPad;
  DrawShashlikModule(tempPad);
  cc_Cal->SaveAs("plots_evan_calib/energypositionprofile_hodo1back.pdf");


  cc_Cal->Clear();

  cc_Cal->cd();
  energyvxy_hodo2->Draw("colz");
  energyvxy_hodo2->GetXaxis()->SetTitle("X2");
  energyvxy_hodo2->GetYaxis()->SetTitle("Y2");
  tempPad = (TPad*) gPad;
  DrawShashlikModule(tempPad);
  cc_Cal->SaveAs("plots_evan_calib/energypositionprofile_hodo2back.pdf");

  //  gMyRootApp->Run(); 
}




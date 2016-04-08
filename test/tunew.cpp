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
  

  

  
  
  bool haverechits = false;
  std::vector<TBRecHit> *rechits=0;
  if(H4tree->GetListOfBranches()->FindObject("tbrechits")) {
    std::cout << " found rechits " << std::endl;
    H4tree->SetBranchAddress("tbrechits",&rechits);
    haverechits = true;
  }
  

  std::vector<vector<float> > X1;
  std::vector<vector<float> > X2;
  std::vector<vector<float> > Y1;
  std::vector<vector<float> > Y2;


  
  
  ///---- loop ----
  for (int i=0; i<nEntries; i++) {


    H4tree->GetEntry(i);

   
    float table_x_shift = tbspill->GetTableX();
    float table_y_shift = tbspill->GetTableY();
   
    if ( ((table_x_reference - table_x_shift) != table_x) || (i == 0) ) {
      //  std::cout << " Table: " << std::endl;
      std::cout << "   x = " << table_x_reference<<"::" <<table_x_shift << " mm " << std::endl;
      std::cout << "   y = " << table_y_reference<<"::" << table_y_shift << " mm " << std::endl;
    }
   
    table_x = table_x_reference - table_x_shift;
    table_y = table_y_reference - table_y_shift;
   
   
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

    X1.push_back(pos_fibers_X1);
    X2.push_back(pos_fibers_X2);
    Y1.push_back(pos_fibers_Y1);
    Y2.push_back(pos_fibers_Y2);

     
  }

  float x1[100];float x2[100];float x3[100];float x4[100];
  float y1[100];float y2[100];float y3[100];float y4[100];



  float w0s[100];
  
  for (int i=0;i<100;i++)
    {w0s[i]=4.5+0.06*(50-i);}
  
  





  for (int i=0;i<100;i++)
    {

  
      TH1F *X_h1_HS1_Cal_front  = new TH1F("X_h1_HS1_Cal_front", "X Hodoscope 1 vs Cal front ", 256, -32, 32);
      TH1F *X_h1_HS2_Cal_front  = new TH1F("X_h1_HS2_Cal_front", "X Hodoscope 2 vs Cal front ", 256, -32, 32);
      TH1F *X_h1_HS1_Cal_back   = new TH1F("X_h1_HS1_Cal_back",  "X Hodoscope 1 vs Cal back " , 256, -32, 32);
      TH1F *X_h1_HS2_Cal_back   = new TH1F("X_h1_HS2_Cal_back",  "X Hodoscope 2 vs Cal back " , 256, -32, 32);
  
  
      TH1F *Y_h1_HS1_Cal_front  = new TH1F("Y_h1_HS1_Cal_front", "Y Hodoscope 1 vs Cal front ", 256, -32, 32);
      TH1F *Y_h1_HS2_Cal_front  = new TH1F("Y_h1_HS2_Cal_front", "Y Hodoscope 2 vs Cal front ", 256, -32, 32);
      TH1F *Y_h1_HS1_Cal_back   = new TH1F("Y_h1_HS1_Cal_back",  "Y Hodoscope 1 vs Cal back " , 256, -32, 32);
      TH1F *Y_h1_HS2_Cal_back   = new TH1F("Y_h1_HS2_Cal_back",  "Y Hodoscope 2 vs Cal back " , 256, -32, 32);
  

      cout<<i<<"th iteration"<<endl;
      cout<<endl;      cout<<endl;      
      CaloCluster* caloCluster = new CaloCluster();
      caloCluster->setW0(w0s[i]);
      caloCluster->setInterCalibrationConstants("data/InterCalibrationConstants.txt");

      for (int j=0; j<nEntries; j++) 
	{
	  H4tree->GetEntry(j);
	  
  
	  //---- calorimeter data
	  if (j==0) caloCluster->setMapperEpoch(tbevent->GetTimeStamp());
	  
	  
	  std::vector<float> caloCluster_position_X_front;
	  std::vector<float> caloCluster_position_Y_front;
	  std::vector<float> caloCluster_Energy_front;
	  std::vector<float> caloCluster_Energies_front;
	  
	  caloCluster->doCalorimeterReconstruction( rechits, 1, 60, doFiber);
	  
	  caloCluster_position_X_front.push_back( caloCluster->getPositionX() );
	  caloCluster_position_Y_front.push_back( caloCluster->getPositionY() );
	  caloCluster_Energy_front.push_back( caloCluster->getEnergy() );
	  caloCluster_Energies_front = caloCluster->getCaloClusterComponents();
	  
	  
	  
	  
  
	  std::vector<float> caloCluster_position_X_back;
	  std::vector<float> caloCluster_position_Y_back;
	  std::vector<float> caloCluster_Energy_back;
	  std::vector<float> caloCluster_Energies_back;;
  
	  
	  caloCluster->doCalorimeterReconstruction( rechits, -1, 60, doFiber);
  
	  caloCluster_position_X_back.push_back( caloCluster->getPositionX() );
	  caloCluster_position_Y_back.push_back( caloCluster->getPositionY() );
	  caloCluster_Energy_back.push_back( caloCluster->getEnergy() );
	  caloCluster_Energies_back = caloCluster->getCaloClusterComponents();
   
	  
	  
     
	  if (X1.at(j).size() ==1 && X2.at(j).size() ==1)
	    
	    {     
	      double sd1=0;
	      double ave1=0;
	      mean_stdev(sd1,ave1,X1.at(j));
	      
	      double sd2=0;
	      double ave2=0;
	      mean_stdev(sd2,ave2,X2.at(j));
	      
	      float c_ave1=mean(caloCluster_position_X_front);    
	      float c_ave2=mean(caloCluster_position_X_back);    
      
	      
	      if( (ave1-ave2) >=-3.92 && (ave1-ave2) <=-1.92  )  {
		//---- X
		
		for (int iCalo = 0; iCalo < caloCluster_position_X_front.size(); iCalo++) {
		  
		  X_h1_HS1_Cal_front->Fill(caloCluster_position_X_front.at(iCalo) - ave1);
	  
		}   
		for (int iCalo = 0; iCalo < caloCluster_position_X_back.size(); iCalo++) {
		  
		  X_h1_HS1_Cal_back->Fill(caloCluster_position_X_back.at(iCalo) - ave1);
		  
		  
		}
		
		
		for (int iCalo = 0; iCalo < caloCluster_position_X_front.size(); iCalo++) {
		  
	
		  X_h1_HS2_Cal_front->Fill(caloCluster_position_X_front.at(iCalo) - ave2);
		  
		}   
		for (int iCalo = 0; iCalo < caloCluster_position_X_back.size(); iCalo++) {
		  
	
		  X_h1_HS2_Cal_back->Fill(caloCluster_position_X_back.at(iCalo) - ave2);
	  
	  
		}
		
	      }
	    }
	  
	  
	  
	  
	  if (Y1.at(j).size() ==1 && Y2.at(j).size() ==1) {
	    
	    
	    double sd1=0;
	    double ave1=0;
	    mean_stdev(sd1,ave1,Y1.at(j));
	    
	    
	    double sd2=0;    
	    double ave2=0;
	    mean_stdev(sd2,ave2,Y2.at(j));
	    
	    float c_ave1=mean(caloCluster_position_Y_front);    
	    float c_ave2=mean(caloCluster_position_Y_back);    
	    
	    
	    if( (ave1-ave2) >=-0.6035 && (ave1-ave2) <=1.3965  )  {	
	      
	      
	      for (int iCalo = 0; iCalo < caloCluster_position_Y_front.size(); iCalo++) {
		
	
		Y_h1_HS1_Cal_front->Fill(caloCluster_position_Y_front.at(iCalo) - ave1);
	
		
	      }
	      for (int iCalo = 0; iCalo < caloCluster_position_Y_back.size(); iCalo++) {
		
	
		Y_h1_HS1_Cal_back->Fill(caloCluster_position_Y_back.at(iCalo) - ave1);
		
	      }
	      
	      
	      for (int iCalo = 0; iCalo < caloCluster_position_Y_front.size(); iCalo++) {
		
	
		Y_h1_HS2_Cal_front->Fill(caloCluster_position_Y_front.at(iCalo) - ave2);
		
		
	      }
	      for (int iCalo = 0; iCalo < caloCluster_position_Y_back.size(); iCalo++) {
		
		
		Y_h1_HS2_Cal_back->Fill(caloCluster_position_Y_back.at(iCalo) - ave2);
		
	      }
	    }
	  }

	}

               
      float mean_xf=X_h1_HS1_Cal_front->GetMean();
      float rms_xf=X_h1_HS1_Cal_front->GetRMS();

  
      //      TF1* fgaus_xf = new TF1("fgaus_xf","gaus(0)+pol2(3)+pol2(3)",0,10);
      TF1* fgaus_xf = new TF1("fgaus_xf","gaus(0)+pol2(3)",mean_xf-3*rms_xf,mean_xf+3*rms_xf);
      fgaus_xf->SetParameter(0,10);
      fgaus_xf->SetParameter(1,0.0);
      fgaus_xf->SetParLimits(1,-10,10);
      fgaus_xf->SetParameter(2,0.3);
      fgaus_xf->SetParLimits(2,0.2,2.0);
  
      fgaus_xf->SetParameter(4,0.0);
      fgaus_xf->SetParameter(5,0.0);
      float mean_xb=(X_h1_HS1_Cal_back->GetMean()+X_h1_HS2_Cal_back->GetMean())/2;
      float rms_xb=(X_h1_HS1_Cal_back->GetRMS()+X_h1_HS2_Cal_back->GetRMS())/2;

 

      //TF1* fgaus_xb = new TF1("fgaus_xb","gaus(0)+pol2(3)+pol2(3)",-8,2);
      TF1* fgaus_xb = new TF1("fgaus_xb","gaus(0)+pol2(3)",mean_xb-4*rms_xb,mean_xb+3*rms_xb);
      fgaus_xb->SetParameter(0,10);
      fgaus_xb->SetParameter(1,0.0);
      fgaus_xb->SetParLimits(1,-10,10);
      fgaus_xb->SetParameter(2,0.3);
      fgaus_xb->SetParLimits(2,0.2,2.0);

      fgaus_xb->SetParameter(4,0.0);
      fgaus_xb->SetParameter(5,0.0);

      float mean_yf=Y_h1_HS1_Cal_front->GetMean();
      float rms_yf=Y_h1_HS1_Cal_front->GetRMS();

 
      //      TF1* fgaus_yf = new TF1("fgaus_yf","gaus(0)+pol2(3)+pol2(3)",10,10);
      TF1* fgaus_yf = new TF1("fgaus_yf","gaus(0)+pol2(3)",mean_yf-3*rms_yf,mean_yf+3*rms_yf);
      fgaus_yf->SetParameter(0,10);
      fgaus_yf->SetParameter(1,0.0);
      fgaus_yf->SetParLimits(1,-10,10);
      fgaus_yf->SetParameter(2,0.3);
      fgaus_yf->SetParLimits(2,0.2,2.0);

      fgaus_yf->SetParameter(4,0.0);
      fgaus_yf->SetParameter(5,0.0);

      float mean_yb=(Y_h1_HS1_Cal_back->GetMean()+Y_h1_HS2_Cal_back->GetMean())/2;
      float rms_yb=(Y_h1_HS1_Cal_back->GetRMS()+Y_h1_HS2_Cal_back->GetRMS())/2;

      //      TF1* fgaus_yb = new TF1("fgaus_yb","gaus(0)+pol2(3)+pol2(3)",2,11);
      TF1* fgaus_yb = new TF1("fgaus_yb","gaus(0)+pol2(3)",mean_yb-4*rms_yb,mean_yb+3*rms_yb);
      fgaus_yb->SetParameter(0,10);
      fgaus_yb->SetParameter(1,0.0);
      fgaus_yb->SetParLimits(1,-10,10);
      fgaus_yb->SetParameter(2,0.3);
      fgaus_yb->SetParLimits(2,0.2,2.0);


  
      fgaus_yb->SetParameter(4,0.0);
      fgaus_yb->SetParameter(5,0.0);
  

      //      X_h1_HS1_Cal_front->Fit("fgaus_xf","QR0");
 

      // X_h1_HS2_Cal_front->Fit("fgaus_xf","QR0");



      X_h1_HS1_Cal_back->Fit("fgaus_xb","QR0");

      X_h1_HS2_Cal_back->Fit("fgaus_xb","QR0");

  
      Y_h1_HS1_Cal_back->Fit("fgaus_yb","QR0"); 

  
      Y_h1_HS2_Cal_back->Fit("fgaus_yb","QR0");
  
      // TF1* funcx1=X_h1_HS1_Cal_front->GetFunction("fgaus_xf");
      // TF1* funcx2=X_h1_HS2_Cal_front->GetFunction("fgaus_xf");
      TF1* funcx3=X_h1_HS1_Cal_back->GetFunction("fgaus_xb");
      TF1* funcx4=X_h1_HS2_Cal_back->GetFunction("fgaus_xb");

      // TF1* funcy1=Y_h1_HS1_Cal_front->GetFunction("fgaus_yf");
      // TF1* funcy2=Y_h1_HS2_Cal_front->GetFunction("fgaus_yf");
      TF1* funcy3=Y_h1_HS1_Cal_back->GetFunction("fgaus_yb");
      TF1* funcy4=Y_h1_HS2_Cal_back->GetFunction("fgaus_yb");


      // x1[i]=funcx1->GetParameter(2);
      // x2[i]=funcx2->GetParameter(2);
      x3[i]=funcx3->GetParameter(2);
      x4[i]=funcx4->GetParameter(2);

      // y1[i]=funcy1->GetParameter(2);
      // y2[i]=funcy2->GetParameter(2);
      y3[i]=funcy3->GetParameter(2);
      y4[i]=funcy4->GetParameter(2);
  
  
      delete X_h1_HS1_Cal_front; 
      delete X_h1_HS2_Cal_front; 
      delete X_h1_HS1_Cal_back;  
      delete X_h1_HS2_Cal_back;  

  
      delete Y_h1_HS1_Cal_front; 
      delete Y_h1_HS2_Cal_front; 
      delete Y_h1_HS1_Cal_back;  
      delete Y_h1_HS2_Cal_back;  
  
        



    }
  






  // TGraph* Tgx=new TGraph(10,w0s,x1);
  // TGraph* Tgx2=new TGraph(10,w0s,x2);
  TGraph* Tgx3=new TGraph(100,w0s,x3);
  TGraph* Tgx4=new TGraph(100,w0s,x4);


  // TGraph* Tgy=new TGraph(4,w0s,y5);
  //TGraph* Tgy2=new TGraph(4,w0s,y2);
  TGraph* Tgy3=new TGraph(100,w0s,y3);
  TGraph* Tgy4=new TGraph(100,w0s,y4);


  TCanvas* cc_X = new TCanvas ("cc_X","",1000,450);
  TCanvas* cc_Y = new TCanvas ("cc_Y","",1000,450);




  cc_X->Divide(2,1);

  /*  cc_X->cd(1);  
      Tgx->Draw("ap");
      Tgx->GetXaxis()->SetTitle("w0");
      Tgx->GetYaxis()->SetTitle("resolution");
      Tgx->SetTitle("X:caloFront-Hodoscope1");  

      cc_X->cd(2);  
      Tgx2->Draw("ap");
      Tgx2->GetXaxis()->SetTitle("w0");
      Tgx2->GetYaxis()->SetTitle("resolution");
      Tgx2->SetTitle("X:caloFront-Hodoscope2");
  */
  cc_X->cd(1);  
  Tgx3->SetLineWidth(3);
  Tgx3->SetMarkerSize(1);
  Tgx3->SetMarkerStyle(21);
  Tgx3->Draw("ap");
  Tgx3->GetXaxis()->SetTitle("w0");
  Tgx3->GetYaxis()->SetTitle("resolution");
  Tgx3->SetTitle("X:caloBack-Hodoscope1");  
  
  cc_X->cd(2);  
  Tgx4->SetLineWidth(3);
  Tgx4->SetMarkerSize(1);
  Tgx4->SetMarkerStyle(21);
  Tgx4->Draw("ap");
  Tgx4->GetXaxis()->SetTitle("w0");
  Tgx4->GetYaxis()->SetTitle("resolution");
  Tgx4->SetTitle("X:caloBack-Hodoscope2");  

  cc_X->SaveAs("/afs/cern.ch/work/n/ndev/CMSSW_7_4_4_patch4/src/H4TestBeam/plot1.pdf");


  cc_Y->Divide(2,1);
  /*
    cc_Y->cd(1);  
    Tgy->Draw("ap");
    Tgy->GetXaxis()->SetTitle("w0");
    Tgy->GetYaxis()->SetTitle("resolution");
    Tgy->SetTitle("Y:caloFront-Hodoscope1");

    cc_Y->cd(2);  
    Tgy2->Draw("ap");
    Tgy2->GetXaxis()->SetTitle("w0");
    Tgy2->GetYaxis()->SetTitle("resolution");
    Tgy2->SetTitle("Y:caloFront-Hodoscope2");  */
    
  cc_Y->cd(1);  
  Tgy3->SetLineWidth(3);
  Tgy3->SetMarkerSize(1);
  Tgy3->SetMarkerStyle(21);
  Tgy3->Draw("ap");
  Tgy3->GetXaxis()->SetTitle("w0");
  Tgy3->GetYaxis()->SetTitle("resolution");
  Tgy3->SetTitle("Y:caloBack-Hodoscope1");
    
  cc_Y->cd(2);  
  Tgy4->SetMarkerSize(1);
  Tgy4->SetMarkerStyle(21);
  Tgy4->SetLineWidth(3);
  Tgy4->Draw("ap");
  Tgy4->GetXaxis()->SetTitle("w0");
  Tgy4->GetYaxis()->SetTitle("resolution");
  Tgy4->SetTitle("Y:caloBack-Hodoscope2");

  cc_Y->SaveAs("/afs/cern.ch/work/n/ndev/CMSSW_7_4_4_patch4/src/H4TestBeam/plot2.pdf");
      
  gMyRootApp->Run(); 
  
}

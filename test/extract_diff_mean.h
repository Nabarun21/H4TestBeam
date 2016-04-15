#include "aux_functions.h"  
#include "TString.h"
#include "TFile.h"
#include "TTree.h"

void  extract_diff_mean(TString fdat, vector <float>& calo_Xf,vector <float>& calo_Xb,vector <float>& calo_Yf,vector <float>& calo_Yb,vector <float>& pos_X1,vector <float>& pos_X2,vector <float>& pos_Y1,vector <float>& pos_Y2, double& x_diff, double& y_diff){

  // extracts the mean of the difference between hodoscope 1 and hodoscope 2

  std::string input_file;
  int maxEvents = -1;
  int doFiber = 0;
  float table_x_reference = 200; //---- mm
  float table_y_reference = 350; //---- mm
  float table_x = 200; //---- mm
  float table_y = 350; //---- mm
  float w0 = 5.0;
 
  //---- configuration
 

  
  //---- read file
  TFile *f= new TFile(fdat); 

  //The tree 
 
  TTree* H4tree = (TTree*)f->Get("t1041");
  
  //---- read file
  int nEntries = H4tree->GetEntries(); 
  std::cout << " nEntries = " << nEntries << std::endl;
  if (maxEvents != -1) nEntries = maxEvents>nEntries ? nEntries : maxEvents ;
  std::cout << " new nEntries = " << nEntries << std::endl;
 
 
  TBSpill* tbspill = new TBSpill();
  TBEvent* tbevent = new TBEvent();
 

  H4tree->SetBranchAddress("tbevent", &tbevent);
  H4tree->SetBranchAddress("tbspill", &tbspill);
 
   
  TH1F *h_X  = new TH1F("h_X", "X Hodoscope 1 vs Hodoscope 2 ", 80, -8, 8);
  TH1F *h_Y  = new TH1F("h_Y", "Y Hodoscope 1 vs Hodoscope 2 ", 80, -8, 8);

  CaloCluster* caloCluster = new CaloCluster();
  caloCluster->setW0(w0);
  caloCluster->setInterCalibrationConstants("data/InterCalibrationConstants.txt");


  bool haverechits = false;
  std::vector<TBRecHit> *rechits=0;
  if(H4tree->GetListOfBranches()->FindObject("tbrechits")) {
    std::cout << " found rechits " << std::endl;
    H4tree->SetBranchAddress("tbrechits",&rechits);
    haverechits = true;}



  for (int i=0;i<nEntries;i++){  
  
    H4tree->GetEntry(i);

    
    if (i==0) {caloCluster->setMapperEpoch(tbevent->GetTimeStamp());
      float table_x_shift = tbspill->GetTableX();
      float table_y_shift = tbspill->GetTableY();
      cout<<"corrosive"<<endl;
      cout<<"Get table X ="<< table_x_shift<<endl;
      cout<<"Get table Y ="<< table_y_shift<<endl;

      if ( ((table_x_reference - table_x_shift) != table_x) || (i == 0) ) {
	std::cout << " Table: " << std::endl;
	std::cout << "   x = " << table_x_reference - table_x_shift << " mm " << std::endl;
	std::cout << "   y = " << table_y_reference - table_y_shift << " mm " << std::endl;
      }
      
      table_x = table_x_reference - table_x_shift;
      table_y = table_y_reference - table_y_shift;
      
     
      std::cout << " Table: " << std::endl;
      std::cout << "   x = " << table_x << " mm " << std::endl;
      std::cout << "   y = " << table_y << " mm " << std::endl;
     
    }
   
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
   

    if (i==8083){
      cout<<"2"<<endl;}   
      
   
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
    if (i==8083){
      cout<<"4"<<endl;}   
   

    //----                                                  table position in mm
    doHodoReconstruction( fibers_X1, n_fibers_X1, pos_fibers_X1, (table_x - 200)); //---- change of coordinate system using numbers from googledoc
    doHodoReconstruction( fibers_X2, n_fibers_X2, pos_fibers_X2, (table_x - 200)); //---- change of coordinate system using numbers from googledoc
    doHodoReconstruction( fibers_Y1, n_fibers_Y1, pos_fibers_Y1, (table_y - 350)); //---- change of coordinate system using numbers from googledoc
    doHodoReconstruction( fibers_Y2, n_fibers_Y2, pos_fibers_Y2, (table_y - 350)); //---- change of coordinate system using numbers from googledoc

 
   
    //---- just hodoscope information

    if (pos_fibers_X1.size()==1 && pos_fibers_X2.size()==1 && caloCluster_position_X_front.size()>0 && caloCluster_position_X_back.size()>0){
      h_X->Fill(pos_fibers_X1.at(0)-pos_fibers_X2.at(0));
      pos_X1.push_back(pos_fibers_X1.at(0));
      pos_X2.push_back(pos_fibers_X2.at(0));
      calo_Xf.push_back(mean(caloCluster_position_X_front));
      calo_Xb.push_back(mean(caloCluster_position_X_back));
    
    }
    if (pos_fibers_Y1.size()==1 && pos_fibers_Y2.size()==1 && caloCluster_position_Y_front.size()>0 && caloCluster_position_Y_back.size()>0){
      h_Y->Fill(pos_fibers_Y1.at(0)-pos_fibers_Y2.at(0));
      pos_Y1.push_back(pos_fibers_Y1.at(0));
      pos_Y2.push_back(pos_fibers_Y2.at(0));
      calo_Yf.push_back(mean(caloCluster_position_Y_front));
      calo_Yb.push_back(mean(caloCluster_position_Y_back));
 
    }
  }

  cout<<"9"<<endl;
  x_diff=h_X->GetMean();
  y_diff=h_Y->GetMean();
}

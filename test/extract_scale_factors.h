#include "include/CaloCluster.h"
#include "include/TBRecHit.h"

void  extract_scale_factors( vector <float> calo_f,vector <float> calo_b,vector<float> pos_1,vector <float> pos_2,vector <float>& n_calo_f,vector <float>& n_calo_b,vector<float>& n_pos_1,vector <float>&n_pos_2,double& front,double& back,double diff){
// extracts the mean of the difference between hodoscope 1 and hodoscope 2


    
  std::vector<float> h1;
  std::vector<float> h2;
  

  std::vector<float> calo;
  std::vector<float> calo_bk;
 
  




  for (int i =0;i<pos_1.size();i++){     
   

    if( (pos_1.at(i)-pos_2.at(i)) >=(diff-1) && (pos_1.at(i)-pos_2.at(i)) <=(diff+1) )  {
     //---- 

      h1.push_back(pos_1.at(i));
      h2.push_back(pos_2.at(i));
      
      calo.push_back(calo_f.at(i));
      calo_bk.push_back(calo_b.at(i));

      n_pos_1.push_back(pos_1.at(i));
      n_pos_2.push_back(pos_2.at(i));
      
      n_calo_f.push_back(calo_f.at(i));
      n_calo_b.push_back(calo_b.at(i));
            
    }
      
  }


  float *ah1=&(h1[0]);
  float *ah2=&(h2[0]);

 
  float *acalo=&(calo[0]);
  
  float *acalo_bk=&(calo_bk[0]);


  TGraph* Tg=new TGraph(h1.size(),ah1,acalo);

  TGraph* Tg2=new TGraph(h2.size(),ah2,acalo);
  
  TGraph* Tg3=new TGraph(h1.size(),ah1,acalo_bk);

  TGraph* Tg4=new TGraph(h2.size(),ah2,acalo_bk);
  
  TF1* lin = new TF1("lin","[0]*x+[1]",-30,30);
  lin->SetParameter(0,1);
  lin->SetParameter(1,0);


  Tg->Fit("lin","R");
  Tg2->Fit("lin","R");
  Tg3->Fit("lin","R");
  Tg4->Fit("lin","R");


  TF1 *f1=Tg->GetFunction("lin");
  TF1 *f2=Tg2->GetFunction("lin");
  TF1 *f3=Tg3->GetFunction("lin");
  TF1 *f4=Tg4->GetFunction("lin");
  
  back = (f4->GetParameter(0)+f3->GetParameter(0))/2;
  front = (f1->GetParameter(0)+f2->GetParameter(0))/2;

    
  }
   

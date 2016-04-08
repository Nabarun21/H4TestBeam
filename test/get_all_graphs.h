#include "extract_diff_mean.h"
#include "extract_scale_factors.h"
#include "extract_TGraph.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include<string>

//void get_all_graphs(string fdat,double shift_x,double  shift_y,TMultiGraph* Mg[8],int s){
void get_all_graphs(string fdat,double shift_x,double  shift_y,TH2F* Mg[8],int s){

   std::string run_num = std::regex_replace( fdat,
				       std::regex("[^0-9]*([0-9]+).*"),
				       string("$1")
				       ); 


     
  double x_diff=0;
  double y_diff=0;

  vector<float> calo_Xf;
  vector<float> calo_Xb;
  vector<float> calo_Yf;
  vector<float> calo_Yb;
  
  vector<float> pos_X1;
  vector<float> pos_X2;
  vector<float> pos_Y1;
  vector<float> pos_Y2;

  cout<<"runnin ext_diff_mean"<<endl;
  
  extract_diff_mean( fdat,  calo_Xf, calo_Xb, calo_Yf, calo_Yb, pos_X1, pos_X2, pos_Y1, pos_Y2, x_diff,y_diff); //form the vectors of position and get the shift between mean hodoscope positions for both coordinates

  cout<<"diffx "<<x_diff<<"diffy "<<y_diff<<endl;


  vector<float> n_calo_Xf;
  vector<float> n_calo_Xb;
  vector<float> n_calo_Yf;
  vector<float> n_calo_Yb;
  
  vector<float> n_pos_X1;
  vector<float> n_pos_X2;
  vector<float> n_pos_Y1;
  vector<float> n_pos_Y2;

  double frontx=1;
  double backx=1;
  double fronty=1;
  double backy=1;

cout<<"runnin ext_scale_fact"<<endl;
  
  extract_scale_factors( calo_Xf, calo_Xb, pos_X1, pos_X2,n_calo_Xf, n_calo_Xb, n_pos_X1,n_pos_X2,frontx, backx,x_diff);
  //X :re form vectors of position after making cut based on relative mean difference and get scaling factor  
  extract_scale_factors( calo_Yf, calo_Yb, pos_Y1, pos_Y2,n_calo_Yf, n_calo_Yb, n_pos_Y1,n_pos_Y2,fronty, backy,y_diff);
  //Y :re form vectors of position after making cut based on relative mean difference and get scaling factor  
  

cout<<"runnin ext_shft_graph"<<endl;
  
cout<<"scales"<<frontx<<" "<<backx<<" "<<fronty<<" "<<backy<<endl;
  
  //X
  extract_shift_Tgraph(n_pos_X1,n_calo_Xf,frontx,shift_x,shift_y,1,Mg,s,run_num);
  extract_shift_Tgraph(n_pos_X2,n_calo_Xf,frontx,shift_x,shift_y,2,Mg,s,run_num);
  extract_shift_Tgraph(n_pos_X1,n_calo_Xb,backx,shift_x,shift_y,3,Mg,s,run_num);
  extract_shift_Tgraph(n_pos_X2,n_calo_Xb,backx,shift_x,shift_y,4,Mg,s,run_num);
  

  //Y
  extract_shift_Tgraph(n_pos_Y1,n_calo_Yf,fronty,shift_x,shift_y,5,Mg,s,run_num);
  extract_shift_Tgraph(n_pos_Y2,n_calo_Yf,fronty,shift_x,shift_y,6,Mg,s,run_num);
  extract_shift_Tgraph(n_pos_Y1,n_calo_Yb,backy,shift_x,shift_y,7,Mg,s,run_num);
  extract_shift_Tgraph(n_pos_Y2,n_calo_Yb,backy,shift_x,shift_y,8,Mg,s,run_num);
    

}
  
  

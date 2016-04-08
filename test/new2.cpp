#include "TH2F.h"
#include <regex> 
#include "get_all_graphs2.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TDR.h"


#define  GLOBAL_X 210.25
#define  GLOBAL_Y 348.00


int main(int argc, char**argv){

  std::string input_file;
 std::string x_positions;
 std::string y_positions;
 int doFiber = 0;
 float table_x_reference = 200; //---- mm
 float table_y_reference = 350; //---- mm
 float table_x = 200; //---- mm
 float table_y = 350; //---- mm
 
 float w0 = 5.0;
 //---- configuration

 // setTDRStyle1(0,0,0,0);
 
//---- configuration

 int c;
 while ((c = getopt (argc, argv, "i:x:y:f")) != -1)
   switch (c)
     {
     case 'i': //---- input
       input_file = string(optarg);
       break;
     case 'x':
       x_positions = string(optarg);
       break;
     case 'y':
       y_positions = string(optarg);
       break;
     case 'f':
       doFiber = atoi(optarg);
       
     case '?':
       if (optopt == 'i' || optopt == 'x' || optopt == 'y'||optopt == 'f')
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
 
 //  std::cout << " Table: " << std::endl;
 //  std::cout << "   x = " << table_x << " mm " << std::endl;
 //std::cout << "   y = " << table_y << " mm " << std::endl;
 
 
 //---- get vector of files
 
 std::vector<std::string> input_files_vector;
 std::vector<std::string> x_positions_vector;
 std::vector<std::string> y_positions_vector;
 
 std::stringstream ss(input_file);
 std::stringstream ss2(x_positions);
 std::stringstream ss3(y_positions); 
 
 cout<<"input files :"<<endl;
 
 std::string token_string;
  while(std::getline(ss, token_string, ',')) {
    std::cout << token_string << '\n';
    input_files_vector.push_back(token_string);
  }
 
  cout<<"x positions"<<endl;
 
  std::string token_string2;
  while(std::getline(ss2, token_string2, ',')) {
    std::cout << token_string2 << '\n';
    x_positions_vector.push_back(token_string2);
  }
  
  cout<<"y positions"<<endl;
  
  std::string token_string3;
  while(std::getline(ss3, token_string3, ',')) {
    std::cout << token_string3 << '\n';
    y_positions_vector.push_back(token_string3);
  }
  

  float x_tempval =GLOBAL_X;
  vector<float> x_pos;
 

 for (int i=0; i<x_positions_vector.size(); i++) 
   {
     char tab[1024];    //temp variable to use in sscanf                                                                                                                    
     strcpy(tab, x_positions_vector.at(i).c_str()); //copy the string into the array of charactars for use in sscanf       
     if (sscanf(tab,"%f",&x_tempval)!=1){
	 cout<<"can't read x_position of the "<<i<<"th input file"<<endl;
	 exit(0);}
       x_pos.push_back(x_tempval);
   }
     

  float y_tempval=GLOBAL_Y;
  vector<float> y_pos;
 

 for (int i=0; i<y_positions_vector.size(); i++) 
   {
     char tab[1024];    //temp variable to use in sscanf                                                                                                                    
     strcpy(tab, y_positions_vector.at(i).c_str()); //copy the string into the array of charactars for use in sscanf       
     if (sscanf(tab,"%f",&y_tempval)!=1){
	 cout<<"can't read y_position of the "<<i<<"th input file"<<endl;
	 exit(0);}
       y_pos.push_back(y_tempval);
   }

  TApplication* gMyRootApp = new TApplication("My ROOT Application", &argc, argv); 



  /*  TMultiGraph* Mg[8];
  char buffer[30];
  for (int k ; k<8;k++)
    {
      sprintf(buffer,"Mg%d",k);
      Mg[k]=new TMultiGraph (buffer,buffer);
    }
  */
  TH2F * histarray[8];

 
 
  char buffer[30];
  for (int k ; k<8;k++)
    {
      sprintf(buffer,"TG%d",k);
      histarray[k]=new TH2F ("buffer","buffer",76, -38, 38, 76, -38, 38);
    }
  


  for (int j=0;j<input_files_vector.size();j++){
    float shift_x=GLOBAL_X-x_pos.at(j);
    float shift_y=GLOBAL_Y-y_pos.at(j);
    cout<<"running_get_all_graphs"<<j<<"th time"<<endl;
    get_all_graphs2(input_files_vector.at(j),shift_x,shift_y,histarray,j);
  }



  TCanvas* X_g = new TCanvas ("X_g","",800,600);
  TCanvas* Y_g = new TCanvas ("Y_g","",800,600);

  X_g->Divide(2,2);
  Y_g->Divide(2,2);
  
  X_g->cd(1)->SetGrid();
  histarray[0]->Draw("Colz");
  histarray[0]->GetXaxis()->SetTitle("Hodo1_x");
  histarray[0]->GetYaxis()->SetTitle("Calox_front");
  gPad->Modified();
  
  X_g->cd(2)->SetGrid();
  histarray[1]->Draw("Colz");
  histarray[1]->GetYaxis()->SetTitle("Calox_front");
  histarray[1]->GetXaxis()->SetTitle("Hodo2_x");
  gPad->Modified();
  
  X_g->cd(3)->SetGrid();
  histarray[2]->Draw("Colz");  
  histarray[2]->GetYaxis()->SetTitle("Hodo1_x");
  histarray[2]->GetYaxis()->SetTitle("Calox_back");
  gPad->Modified();
  
  X_g->cd(4)->SetGrid();
  histarray[3]->Draw("Colz");
  histarray[3]->GetYaxis()->SetTitle("Calox_back");
  histarray[3]->GetXaxis()->SetTitle("Hodo2_x");
  gPad->Modified();
    
  Y_g->cd(1)->SetGrid();   
  histarray[4]->Draw("Colz");
  histarray[4]->GetXaxis()->SetTitle("Hodo1_y");
  histarray[4]->GetYaxis()->SetTitle("Caloy_front"); 
  gPad->Modified();
  
  Y_g->cd(2)->SetGrid();   
  histarray[5]->Draw("Colz");
  histarray[5]->GetXaxis()->SetTitle("Hodo2_y"); 
  histarray[5]->GetYaxis()->SetTitle("Caloy_front"); 
  gPad->Modified();
 
  Y_g->cd(3)->SetGrid();   
  histarray[6]->Draw("Colz");
  histarray[6]->GetXaxis()->SetTitle("Hodo1_y");
  histarray[6]->GetYaxis()->SetTitle("Caloy_back");
  gPad->Modified();


  Y_g->cd(4)->SetGrid();   
  histarray[7]->Draw("Colz");
  histarray[7]->GetYaxis()->SetTitle("Caloy_back");
  histarray[7]->GetXaxis()->SetTitle("Hodo2_y");
  gPad->Modified();



  /*
  X_g->Divide(2,2);
  Y_g->Divide(2,2);

  X_g->cd(1)->SetGrid();
  Mg[0]->Draw("AP");
  Mg[0]->GetXaxis()->SetTitle("Hodo1_x");
  Mg[0]->GetYaxis()->SetTitle("Calox_front");
  gPad->Modified();
  X_g->cd(1)->BuildLegend();

  X_g->cd(2)->SetGrid();
  Mg[1]->Draw("AP");
  Mg[1]->GetYaxis()->SetTitle("Calox_front");
  Mg[1]->GetXaxis()->SetTitle("Hodo2_x");
  gPad->Modified();
  X_g->cd(2)->BuildLegend();
  
  X_g->cd(3)->SetGrid();
  Mg[2]->Draw("AP");  
  Mg[2]->GetYaxis()->SetTitle("Hodo1_x");
  Mg[2]->GetYaxis()->SetTitle("Calox_back");
  gPad->Modified();
  X_g->cd(3)->BuildLegend();
  
  X_g->cd(4)->SetGrid();
  Mg[3]->Draw("AP");
  Mg[3]->GetYaxis()->SetTitle("Calox_back");
  Mg[3]->GetXaxis()->SetTitle("Hodo2_x");
  gPad->Modified();
  X_g->cd(4)->BuildLegend();
    
  Y_g->cd(1)->SetGrid();   
  Mg[4]->Draw("AP");
  Mg[4]->GetXaxis()->SetTitle("Hodo1_y");
  Mg[4]->GetYaxis()->SetTitle("Caloy_front"); 
  gPad->Modified();
  Y_g->cd(1)->BuildLegend();
  
  Y_g->cd(2)->SetGrid();   
  Mg[5]->Draw("AP");
  Mg[5]->GetXaxis()->SetTitle("Hodo2_y"); 
  Mg[5]->GetYaxis()->SetTitle("Caloy_front"); 
  gPad->Modified();
  Y_g->cd(2)->BuildLegend();
 
  Y_g->cd(3)->SetGrid();   
  Mg[6]->Draw("AP");
  Mg[6]->GetXaxis()->SetTitle("Hodo1_y");
  Mg[6]->GetYaxis()->SetTitle("Caloy_back");
  gPad->Modified();
  Y_g->cd(3)->BuildLegend();

  Y_g->cd(4)->SetGrid();   
  Mg[7]->Draw("AP");
  Mg[7]->GetYaxis()->SetTitle("Caloy_back");
  Mg[7]->GetXaxis()->SetTitle("Hodo2_y");
  gPad->Modified();
  Y_g->cd(4)->BuildLegend();

  */

  gMyRootApp->Run(); 
  
}

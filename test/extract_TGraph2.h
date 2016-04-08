#include "include/CaloCluster.h"
#include "include/TBRecHit.h"
#include "TGraph.h"
#include "TMultiGraph.h"
//void  extract_shift_Tgraph2( vector<float> hodo, vector<float> calo, double scale ,double shiftx,double shifty,int k,TMultiGraph* Mg[8],int s ,TString run_num){
void  extract_shift_Tgraph2( vector<float> hodo, vector<float> calo, double scale ,double shiftx,double shifty,int k,TH2F* Mg[2],int s ,TString run_num){

// extracts the  of the difference between hodoscope 1 and hodoscope 2



//   float y_array[hodo.size()];
 // float x_array[hodo.size()];


if((4<k) &&(k<9))
{
 for (int j=0 ;j<hodo.size();j++)
    {
	Mg[k-1]->Fill( (hodo.at(j))+shifty,(calo.at(j)/scale));

	}	
}
else
{	
for (int j=0 ;j<hodo.size();j++)
    {
	Mg[k-1]->Fill( (hodo.at(j)+shiftx),(calo.at(j)/scale));

	}	
}
/*  for (int j=0 ;j<hodo.size();j++)
    {
      y_array[j]=(calo.at(j)/scale)+shifty;
      x_array[j]=hodo.at(j)+shiftx;
      
    }
  


  TGraph* Tg=new TGraph(hodo.size(),x_array,y_array);
  Tg->SetMarkerColor(s+1);
  Tg->SetMarkerStyle(s+2);
  Tg->SetTitle(run_num);
   cout<<"Mojo"<<endl;
  Mg[k-1]->Add(Tg);*/
}

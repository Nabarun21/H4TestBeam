#include "CaloCluster.h"
#include <iostream>

//---- Author: Andrea Massironi

//---- position units here is "mm"
//---- energy   units here is "GeV"





//---- distance definition in calorimeter position
void CaloCluster::Reset () {
 _pos_x = 0.;
 _pos_y = 0.;
 _energy = 0.; 
}



//---- distance definition in calorimeter position
float CaloCluster::DR (float x1, float x2, float y1, float y2) {
 return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}





//---- Reconstruct Calorimeter clusters
void CaloCluster::doCalorimeterReconstruction( std::vector<TBRecHit>* rechits, int face, float maxDR, int fiberLevel) { 
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
  _mapper->ChannelID2ModuleFiber(channelID,moduleID,fiberID);  // get module and fiber IDs
  
  double x,y;
  if (fiberLevel == 0) _mapper->ModuleXY(moduleID,x,y);
  else                 _mapper->FiberXY(fiberID, x, y);
  
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
 
 
 //---- do clustering ----
 //  int num_clusters = 0;
 float x_cluster_logE = 0;
 float y_cluster_logE = 0;
 float weight_cluster = 0;
 float energy_cluster = 0;
 
 float x_cluster_max = -1000;
 float y_cluster_max = -1000;
 
 //---- first calculate cluster energy
 for( std::map < float, std::pair<float, float> >::iterator ii=map_of_calo_clusters.begin(); ii!=map_of_calo_clusters.end(); ii++) {
  //   std::cout << " energy = " << - ii->first << std::endl;
  if (x_cluster_max == -1000) {
   x_cluster_max = (ii->second.first); //---- seed position X
   y_cluster_max = (ii->second.second); //---- seed position X   
   float energy_temp = - ii->first;
   energy_cluster = energy_cluster + energy_temp;
  }
  else {
   if (DR(ii->second.first, x_cluster_max, ii->second.second, y_cluster_max) < maxDR) {
    float energy_temp = - ii->first;
    energy_cluster = energy_cluster + energy_temp;
    //   num_clusters++;
   }
  }
 }
 
 //---- then calculate position 
 //  std::cout<< " new cluster " << std::endl;
 if (energy_cluster > 0) {
  for( std::map < float, std::pair<float, float> >::iterator ii = map_of_calo_clusters.begin(); ii != map_of_calo_clusters.end(); ii++) {
   //   std::cout << " energy = " << - ii->first << std::endl;
   if (DR(ii->second.first, x_cluster_max, ii->second.second, y_cluster_max) < maxDR) {
    float energy_temp = - ii->first;
    float wi = (_w0 + log(energy_temp/energy_cluster));
    if (wi > 0) {
     x_cluster_logE = x_cluster_logE + (ii->second.first)  * wi ;
     y_cluster_logE = y_cluster_logE + (ii->second.second) * wi ;
     weight_cluster = weight_cluster + wi;
    }
   }
  }
 }
 
 float x_cluster_final = 0;
 float y_cluster_final = 0;
 
 if (weight_cluster != 0) {
  x_cluster_final = x_cluster_logE / weight_cluster;
  y_cluster_final = y_cluster_logE / weight_cluster;
  
  _pos_x = x_cluster_final;
  _pos_y = y_cluster_final;
  _energy = energy_cluster;
 }
 
 //  std::cout << " num_clusters = " << num_clusters << std::endl;
} 










#!/bin/env python                                                                                                                                                                    
import os
PATH='/afs/cern.ch/work/n/ndev/public/merged_files'


path=PATH
job_num=2170
def add_shift():
    for file in os.listdir(path):
        print file
        job_num=file.split('_')[1]

adage= "AddInformation -i "+path+"/merged_"+str(job_num)+"_reco.root"+" -o "+"new"+str(job_num)+".root"

os.popen(adage)                                                                         

add_shift()

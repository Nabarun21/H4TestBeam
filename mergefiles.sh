#!/bin/bash

n=2051;
max=2060;
while [ "$n" -le "$max" ]; do
    if [ -d /afs/cern.ch/user/e/ewolfe/public/testbeam4/hodoscope/"$n" ] && [ -f reco/rec_capture_"$n"_reco.root ] ; then
	AddHodoInfo -s reco/rec_capture_"$n"_reco.root -b $(echo /afs/cern.ch/user/e/ewolfe/public/testbeam4/hodoscope/"$n"/*.root |sed 's/\ /,/g') -o /afs/cern.ch/work/n/ndev/public/merged_files/merged_"$n"_reco.root
    else
	echo "Either the reco or hodo files does not exist"
    fi
 #xrdcp root://eoscms//eos/cms/store/user/ewolfe/TESTBEAM/reco/rec_capture_"$n"_reco.root .
#echo $root://eoscms//eos/cms/store/user/ewolfe/TESTBEAM/reco/rec_capture_"$n"_reco.root
  n=`expr "$n" + 1`;

done

#
# Put this in the folder you want your files in.
# Then run the commands below:
# > chmod +x CopyFilesFromEOS.sh
# > ./CopyFilesFromEOS.sh
#






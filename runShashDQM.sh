#!/bin/bash                                                                                                                                                                          
#set defaukt options
inputDirPade="../../DAQ_2015"
inputDirHodo="../../DAQ_2015/hodoscope"
outputDir="."
outputPlotsDir="offlineDQMplots"
run="2422"


TEMP=`getopt -o i:h:o:p:r: --long inputDirPade:,inputDirHodo:,outputDir:,outputPlotsDir:,run: -n 'runShashDQM.sh' -- "$@"`
if [ $? != 0 ] ; then echo "Options are wrong..." >&2 ; exit 1 ; fi

eval set -- "$TEMP"


while true; do
case "$1" in
-i | --inputDirPade ) inputDirPade="$2"; shift 2 ;;
-h | --inputDirHodo ) inputDirHodo="$2"; shift 2 ;;
-o | --outputDir ) outputDir="$2"; shift 2 ;;
-p | --outputPlotsDir ) outputPlotsDir="$2"; shift 2 ;;
-r | --run ) run="$2"; shift 2;;
-- ) shift; break ;;
* ) break ;;
esac
done
echo $inputDirPade
echo $inputDirHodo
echo $outputDir
echo $outputPlotsDir
TBTreeMaker.py -P $inputDirPade"/rec_capture_"$run".txt.bz2" -o $outputDir
echo " "
echo " "
echo " "
echo " "
echo "running merging script"
echo " "
echo " "
echo " "
echo " "
AddHodoInfo -s $outputDir/rec_capture_$run.root -b $(echo $inputDirHodo/$run/?.root |sed 's/\ /,/g'),$(echo $inputDirHodo/$run/??.root |sed 's/\ /,/g'),$(echo $inputDirHodo/$run/???.root |sed 's/\ /,/g') -o $outputDir/merged_$run.root
echo " "
echo " "
echo " "
echo " "
echo "running plotter"
DQMPlots -i   $outputDir/merged_$run.root -o $outputPlotsDir
echo " "
echo " "


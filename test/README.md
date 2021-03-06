Trees location
====

Hodoscope and Shashlik merged 

    /eos/cms/store/group/comm_ecal/upgrade/testbeam/H4_Fall2014/LYSO/Merged
    /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/H4_Fall2014/LYSO/Merged
    
Hodoscope and Shashlik merged, table position and energy added by automatic script

    /eos/cms/store/group/comm_ecal/upgrade/testbeam/H4_Fall2014/LYSO/Merged_and_table_and_energy
    /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/H4_Fall2014/LYSO/Merged_and_table_and_energy

    
    
Add Hodoscope information
====

Automatic script to add Hodoscope information

        
    ls ../DAQ/rec_capture_*_reco.root      | sed s:../DAQ/rec_capture_:: | sed s/_reco.root// | awk '{print "./bin/AddHodoInfo -s ../DAQ/rec_capture_"$1"_reco.root  -b ../DAQ/hodoscope/"$1"/?.root,../DAQ/hodoscope/"$1"/??.root,../DAQ/hodoscope/"$1"/???.root  -o  ../DAQ/rec_capture_"$1"_reco_beam.root"}'
    ls ../DAQ/rec_capture_*_reco_beam.root | sed s:../DAQ/rec_capture_:: | sed s/_reco_beam.root// | awk '{print "./bin/AddInformation -i ../DAQ/rec_capture_"$1"_reco_beam.root  -o rec_capture_"$1"_reco_beam_table_energy.root"}'
    ls ../DAQ/rec_capture_*_reco_beam.root | sed s:../DAQ/rec_capture_:: | sed s/_reco_beam.root// | awk '{print "./bin/AddInformation -i ../DAQ/rec_capture_"$1"_reco_beam.root  -o rec_capture_"$1"_reco_beam_table_energy.root  -f 1"}'
    ls ../DAQ/rec_capture_*_reco_beam.root | sed s:../DAQ/rec_capture_:: | sed s/_reco_beam.root// | awk '{print "./bin/AddInformation -i ../DAQ/rec_capture_"$1"_reco_beam.root  -o rec_capture_"$1"_reco_beam_table_energy.root  -f 1  -w 4.8"}'

    ls /tmp/amassiro/eos/cms/store/group/comm_ecal/upgrade/testbeam/H4_Fall2014/LYSO/Merged/rec_capture_*_reco_beam.root \
        | sed s:/tmp/amassiro/eos/cms/store/group/comm_ecal/upgrade/testbeam/H4_Fall2014/LYSO/Merged/rec_capture_:: | sed s/_reco_beam.root// | awk '{print "./bin/AddInformation -i /tmp/amassiro/eos/cms/store/group/comm_ecal/upgrade/testbeam/H4_Fall2014/LYSO/Merged/rec_capture_"$1"_reco_beam.root  -o /tmp/amassiro/eos/cms/store/group/comm_ecal/upgrade/testbeam/H4_Fall2014/LYSO/Merged_and_table_and_energy/rec_capture_"$1"_reco_beam_table_energy.root  -f 1  -w 4.8"}'
    ls /tmp/amassiro/eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/H4_Fall2014/LYSO/Merged/rec_capture_*_reco_beam.root \
        | sed s:/tmp/amassiro/eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/H4_Fall2014/LYSO/Merged/rec_capture_:: | sed s/_reco_beam.root// | awk '{print "./bin/AddInformation -i /tmp/amassiro/eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/H4_Fall2014/LYSO/Merged/rec_capture_"$1"_reco_beam.root  -o /tmp/amassiro/eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/H4_Fall2014/LYSO/Merged_and_table_and_energy/rec_capture_"$1"_reco_beam_table_energy.root  -f 1  -w 4.8"}'

    
    
Add Table position and energy information
====

    ./bin/AddInformation -i ../DAQ/rec_capture_1098_reco_beam.root -o rec_capture_1098_reco_beam_table_energy.root    
    ./bin/AddInformation -i ../DAQ/rec_capture_1163_reco_beam.root -o rec_capture_1163_reco_beam_table_energy.root
    ./bin/AddInformation -i ../DAQ/rec_capture_1164_reco_beam.root -o rec_capture_1164_reco_beam_table_energy.root
    ./bin/AddInformation -i ../DAQ/rec_capture_1100_reco_beam.root -o rec_capture_1100_reco_beam_table_energy.root
    ./bin/AddInformation -i ../DAQ/rec_capture_1101_reco_beam.root -o rec_capture_1101_reco_beam_table_energy.root
    ./bin/AddInformation -i ../DAQ/rec_capture_1102_reco_beam.root -o rec_capture_1102_reco_beam_table_energy.root
    ./bin/AddInformation -i ../DAQ/rec_capture_1081_reco_beam.root -o rec_capture_1081_reco_beam_table_energy.root

    ./bin/AddInformation -i ../DAQ/rec_capture_1104_reco_beam.root -o rec_capture_1104_reco_beam_table_energy.root
    ./bin/AddInformation -i ../DAQ/rec_capture_1110_reco_beam.root -o rec_capture_1110_reco_beam_table_energy.root
    ./bin/AddInformation -i ../DAQ/rec_capture_1111_reco_beam.root -o rec_capture_1111_reco_beam_table_energy.root
    ./bin/AddInformation -i ../DAQ/rec_capture_1112_reco_beam.root -o rec_capture_1112_reco_beam_table_energy.root

      
Plot beam and calorimeter information
====

    
    ./bin/PlotHodoAndShashlik -i rec_capture_1098_reco_beam_table_energy.root
    ./bin/PlotHodoAndShashlik -i rec_capture_1102_reco_beam_table_energy.root
    ./bin/PlotHodoAndShashlik -i rec_capture_1163_reco_beam_table_energy.root
    ./bin/PlotHodoAndShashlik -i rec_capture_1164_reco_beam_table_energy.root
    ./bin/PlotHodoAndShashlik -i rec_capture_1100_reco_beam_table_energy.root
    ./bin/PlotHodoAndShashlik -i rec_capture_1101_reco_beam_table_energy.root
    ./bin/PlotHodoAndShashlik -i rec_capture_1102_reco_beam_table_energy.root
    ./bin/PlotHodoAndShashlik -i rec_capture_1081_reco_beam_table_energy.root

    ./bin/PlotHodoAndShashlik -i rec_capture_1104_reco_beam_table_energy.root
    ./bin/PlotHodoAndShashlik -i rec_capture_1110_reco_beam_table_energy.root
    ./bin/PlotHodoAndShashlik -i rec_capture_1111_reco_beam_table_energy.root
    ./bin/PlotHodoAndShashlik -i rec_capture_1112_reco_beam_table_energy.root

    ./bin/PlotHodoAndShashlik -i rec_capture_1110_reco_beam_table_energy.root,rec_capture_1111_reco_beam_table_energy.root
    
    ./bin/PlotHodoAndShashlik -i rec_capture_1026_reco_beam_table_energy.root,rec_capture_1012_reco_beam_table_energy.root,rec_capture_1013_reco_beam_table_energy.root,rec_capture_1081_reco_beam_table_energy.root   -f 1   -w 4.8
    ./bin/PlotHodoAndShashlik -i rec_capture_1026_reco_beam_table_energy.root,rec_capture_1012_reco_beam_table_energy.root,rec_capture_1013_reco_beam_table_energy.root,rec_capture_1081_reco_beam_table_energy.root   -f 0   -w 3.35

    50 GeV
    ./bin/PlotHodoAndShashlik -i rec_capture_1081_reco_beam_table_energy.root   -f 0   -w 3.35
    ./bin/PlotHodoAndShashlik -i Merged_and_table_and_energy/rec_capture_1081_reco_beam_table_energy.root   -f 0   -w 3.35
    
    100 GeV
    ./bin/PlotHodoAndShashlik -i rec_capture_1026_reco_beam_table_energy.root   -f 0   -w 3.35
    ./bin/PlotHodoAndShashlik -i Merged_and_table_and_energy/rec_capture_1026_reco_beam_table_energy.root   -f 0   -w 3.35
    
    100 GeV
    ./bin/PlotHodoAndShashlik -i rec_capture_1012_reco_beam_table_energy.root   -f 0   -w 3.35
    ./bin/PlotHodoAndShashlik -i Merged_and_table_and_energy/rec_capture_1012_reco_beam_table_energy.root   -f 0   -w 3.35
    
    100 GeV
    ./bin/PlotHodoAndShashlik -i rec_capture_1013_reco_beam_table_energy.root   -f 0   -w 3.35

    ./bin/PlotHodoAndShashlik -i rec_capture_1020_reco_beam_table_energy.root
    ./bin/PlotHodoAndShashlik -i rec_capture_1026_reco_beam_table_energy.root
    ./bin/PlotHodoAndShashlik -i rec_capture_1012_reco_beam_table_energy.root
    ./bin/PlotHodoAndShashlik -i rec_capture_1013_reco_beam_table_energy.root

    ./bin/PlotHodoAndShashlik -i rec_capture_1024_reco_beam_table_energy.root
    ./bin/PlotHodoAndShashlik -i rec_capture_1025_reco_beam_table_energy.root
    ./bin/PlotHodoAndShashlik -i rec_capture_1017_reco_beam_table_energy.root
         

Calculate position resolution
====

    ./bin/PositionResolution -i Merged_and_table_and_energy/rec_capture_1012_reco_beam_table_energy.root   -f 0   -w 3.35
    
    
Plot Hodoscope information
====

    ./bin/PlotHodo  -i  ../DAQ/rec_capture_9999.root
    ./bin/PlotHodo  -i  ../DAQ/rec_capture_1163_reco_beam.root
    
    
Plot Hodoscope and Shashlik information
====
    
    ./bin/PlotHodoAndShashlik -i mix_1140_reco.root
    ./bin/PlotHodoAndShashlik -i ../DAQ/rec_capture_1003_reco_beam.root
    ./bin/PlotHodoAndShashlik -i ../DAQ/rec_capture_1224_reco_beam.root
    ./bin/PlotHodoAndShashlik -i ../DAQ/rec_capture_1253_reco_beam.root
    ./bin/PlotHodoAndShashlik -i ../DAQ/rec_capture_1244_reco_beam.root
    ./bin/PlotHodoAndShashlik -i ../DAQ/rec_capture_1140_reco_beam.root
    ./bin/PlotHodoAndShashlik -i ../DAQ/rec_capture_1211_reco_beam.root
    ./bin/PlotHodoAndShashlik -i ../DAQ/rec_capture_1163_reco_beam.root   -f  1
    ./bin/PlotHodoAndShashlik -i ../DAQ/rec_capture_1163_reco_beam.root   -f  0
    ./bin/PlotHodoAndShashlik -i ../DAQ/rec_capture_1221_reco_beam.root
    ./bin/PlotHodoAndShashlik -i ../DAQ/rec_capture_1121_reco_beam.root   -f  1
    
    ./bin/PlotHodoAndShashlik -i ../DAQ/rec_capture_1163_reco_beam.root   -f  0   -x 192.75  -y   357.25
    ./bin/PlotHodoAndShashlik -i ../DAQ/rec_capture_1142_reco_beam.root   -f  0   -x 207.25  -y   342.75
    
    ./bin/PlotHodoAndShashlik -i ../DAQ/rec_capture_1096_reco_beam.root   -f  0   -x 163.75  -y   342.75
    
    ./bin/PlotHodoAndShashlik -i ../DAQ/rec_capture_1101_reco_beam.root   -f  0   -x 192.75  -y   357.25    
    ./bin/PlotHodoAndShashlik -i ../DAQ/rec_capture_1101_reco_beam.root   -f  0   -x 207.25  -y   357.25
        
    ./bin/PlotHodoAndShashlik -i ../DAQ/rec_capture_1098_reco_beam.root   -f  0   -x 178.25  -y   357.25

    
    ./bin/PlotHodoAndShashlik -i rec_capture_1098_reco_beam_table_energy.root   -f  0

    
    
    
test
====

    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1140_reco.root  -b ../DAQ/hodoscope/1140/1.root   -o mix_1140.root


    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1224.root  -b ../DAQ/hodoscope/1224/?.root,../DAQ/hodoscope/1224/??.root,../DAQ/hodoscope/1224/???.root  -o  ../DAQ/rec_capture_9999.root
    
    
    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1003_reco.root  -b ../DAQ/hodoscope/1003/?.root,../DAQ/hodoscope/1003/??.root,../DAQ/hodoscope/1003/???.root  -o  ../DAQ/rec_capture_1003_reco_beam.root
    
    
    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1224_reco.root  -b ../DAQ/hodoscope/1224/?.root,../DAQ/hodoscope/1224/??.root,../DAQ/hodoscope/1224/???.root  -o  ../DAQ/rec_capture_1224_reco_beam.root
    
    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1253_reco.root  -b ../DAQ/hodoscope/1253/?.root,../DAQ/hodoscope/1253/??.root,../DAQ/hodoscope/1253/???.root  -o  ../DAQ/rec_capture_1253_reco_beam.root
 
    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1244_reco.root  -b ../DAQ/hodoscope/1244/?.root,../DAQ/hodoscope/1244/??.root,../DAQ/hodoscope/1244/???.root  -o  ../DAQ/rec_capture_1244_reco_beam.root

    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1140_reco.root  -b ../DAQ/hodoscope/1140/?.root,../DAQ/hodoscope/1140/??.root,../DAQ/hodoscope/1140/???.root  -o  ../DAQ/rec_capture_1140_reco_beam.root

    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1211_reco.root  -b ../DAQ/hodoscope/1211/?.root,../DAQ/hodoscope/1211/??.root,../DAQ/hodoscope/1211/???.root  -o  ../DAQ/rec_capture_1211_reco_beam.root

    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1163_reco.root  -b ../DAQ/hodoscope/1163/?.root,../DAQ/hodoscope/1163/??.root,../DAQ/hodoscope/1163/???.root  -o  ../DAQ/rec_capture_1163_reco_beam.root

    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1221_reco.root  -b ../DAQ/hodoscope/1221/?.root,../DAQ/hodoscope/1221/??.root,../DAQ/hodoscope/1221/???.root  -o  ../DAQ/rec_capture_1221_reco_beam.root

    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1121_reco.root  -b ../DAQ/hodoscope/1121/?.root,../DAQ/hodoscope/1121/??.root,../DAQ/hodoscope/1121/???.root  -o  ../DAQ/rec_capture_1121_reco_beam.root

    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1142_reco.root  -b ../DAQ/hodoscope/1142/?.root,../DAQ/hodoscope/1142/??.root,../DAQ/hodoscope/1142/???.root  -o  ../DAQ/rec_capture_1142_reco_beam.root

    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1096_reco.root  -b ../DAQ/hodoscope/1096/?.root,../DAQ/hodoscope/1096/??.root,../DAQ/hodoscope/1096/???.root  -o  ../DAQ/rec_capture_1096_reco_beam.root

    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1101_reco.root  -b ../DAQ/hodoscope/1101/?.root,../DAQ/hodoscope/1101/??.root,../DAQ/hodoscope/1101/???.root  -o  ../DAQ/rec_capture_1101_reco_beam.root

    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1102_reco.root  -b ../DAQ/hodoscope/1102/?.root,../DAQ/hodoscope/1102/??.root,../DAQ/hodoscope/1102/???.root  -o  ../DAQ/rec_capture_1102_reco_beam.root

    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1098_reco.root  -b ../DAQ/hodoscope/1098/?.root,../DAQ/hodoscope/1098/??.root,../DAQ/hodoscope/1098/???.root  -o  ../DAQ/rec_capture_1098_reco_beam.root
    
    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1081_reco.root  -b ../DAQ/hodoscope/1081/?.root,../DAQ/hodoscope/1081/??.root,../DAQ/hodoscope/1081/???.root  -o  ../DAQ/rec_capture_1081_reco_beam.root
    
    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1164_reco.root  -b ../DAQ/hodoscope/1164/?.root,../DAQ/hodoscope/1164/??.root,../DAQ/hodoscope/1164/???.root  -o  ../DAQ/rec_capture_1164_reco_beam.root

    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1100_reco.root  -b ../DAQ/hodoscope/1100/?.root,../DAQ/hodoscope/1100/??.root,../DAQ/hodoscope/1100/???.root  -o  ../DAQ/rec_capture_1100_reco_beam.root

    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1101_reco.root  -b ../DAQ/hodoscope/1101/?.root,../DAQ/hodoscope/1101/??.root,../DAQ/hodoscope/1101/???.root  -o  ../DAQ/rec_capture_1101_reco_beam.root
    
    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1102_reco.root  -b ../DAQ/hodoscope/1102/?.root,../DAQ/hodoscope/1102/??.root,../DAQ/hodoscope/1102/???.root  -o  ../DAQ/rec_capture_1102_reco_beam.root

    
    
    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1104_reco.root  -b ../DAQ/hodoscope/1104/?.root,../DAQ/hodoscope/1104/??.root,../DAQ/hodoscope/1104/???.root  -o  ../DAQ/rec_capture_1104_reco_beam.root
    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1110_reco.root  -b ../DAQ/hodoscope/1110/?.root,../DAQ/hodoscope/1110/??.root,../DAQ/hodoscope/1110/???.root  -o  ../DAQ/rec_capture_1110_reco_beam.root
    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1111_reco.root  -b ../DAQ/hodoscope/1111/?.root,../DAQ/hodoscope/1111/??.root,../DAQ/hodoscope/1111/???.root  -o  ../DAQ/rec_capture_1111_reco_beam.root
    ./bin/AddHodoInfo -s ../DAQ/rec_capture_1112_reco.root  -b ../DAQ/hodoscope/1112/?.root,../DAQ/hodoscope/1112/??.root,../DAQ/hodoscope/1112/???.root  -o  ../DAQ/rec_capture_1112_reco_beam.root
    
 
 
Transform intercalibration constants
====

    for (int i=0; i<32; i++) {if (i<16) {std::cout << i+1 << "   " << cc2[i] << std::endl;} else {std::cout << -(i+1-16) << "   " << cc2[i] << std::endl;}}


Check beam energy
====
    
    r99t test/rootLogon.C
    TChain* tree = new TChain("t1041");
    tree->Add("Merged_and_table_and_energy/rec_capture_*_reco_beam_table_energy.root");
    tree->Draw("tbspill.GetMomentum()")
    
    
    ls Merged_and_table_and_energy/rec_capture_*_reco_beam_table_energy.root | tr "_" " " |  awk '{print "r99t -q test/rootLogon.C Merged_and_table_and_energy/rec_capture_"$7"_reco_beam_table_energy.root  test/GetBeamInformation.C\\("$7"\\) "}'
    ls Merged_and_table_and_energy/rec_capture_*_reco_beam_table_energy.root | tr "_" " " |  awk '{print "root -l -q test/rootLogon.C Merged_and_table_and_energy/rec_capture_"$7"_reco_beam_table_energy.root  test/GetBeamInformation.C\\("$7"\\) "}' | /bin/sh
    


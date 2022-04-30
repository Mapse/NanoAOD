
# Crab Report
crab report ../../CharmoniumRun2016UL/crab_CharmoniumRun2016B_21Feb2020_ver1_UL2016_HIPM-v1_AOD
crab report ../../CharmoniumRun2016UL/crab_CharmoniumRun2016B_21Feb2020_ver2_UL2016_HIPM-v1_AOD
crab report ../../CharmoniumRun2016UL/crab_CharmoniumRun2016C_21Feb2020_UL2016_HIPM-v1_AOD
crab report ../../CharmoniumRun2016UL/crab_CharmoniumRun2016D_21Feb2020_UL2016_HIPM-v1_AOD
crab report ../../CharmoniumRun2016UL/crab_CharmoniumRun2016E_21Feb2020_UL2016_HIPM-v1_AOD
crab report ../../CharmoniumRun2016UL/crab_CharmoniumRun2016F_21Feb2020_UL2016-v1_AOD
crab report ../../CharmoniumRun2016UL/crab_CharmoniumRun2016F_21Feb2020_UL2016_HIPM-v1_AOD
crab report ../../CharmoniumRun2016UL/crab_CharmoniumRun2016G_21Feb2020_UL2016-v1_AOD
crab report ../../CharmoniumRun2016UL/crab_CharmoniumRun2016H_21Feb2020_UL2016-v1_AOD


cd ../../

brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i CharmoniumRun2016UL/crab_CharmoniumRun2016B_21Feb2020_ver1_UL2016_HIPM-v1_AOD/results/processedLumis.json > lumi_charmonium_2016B_ver1_HIPM-v1.txt
brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i CharmoniumRun2016UL/crab_CharmoniumRun2016B_21Feb2020_ver2_UL2016_HIPM-v1_AOD/results/processedLumis.json > lumi_charmonium_2016B_ver2_HIPM-v1.txt
brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i CharmoniumRun2016UL/crab_CharmoniumRun2016C_21Feb2020_UL2016_HIPM-v1_AOD/results/processedLumis.json > lumi_charmonium_2016C_HIPM.txt
brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i CharmoniumRun2016UL/crab_CharmoniumRun2016D_21Feb2020_UL2016_HIPM-v1_AOD/results/processedLumis.json > lumi_charmonium_2016D_HIPM.txt
brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i CharmoniumRun2016UL/crab_CharmoniumRun2016E_21Feb2020_UL2016_HIPM-v1_AOD/results/processedLumis.json > lumi_charmonium_2016E_HIPM.txt
brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i CharmoniumRun2016UL/crab_CharmoniumRun2016F_21Feb2020_UL2016-v1_AOD/results/processedLumis.json > lumi_charmonium_2016F.txt
brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i CharmoniumRun2016UL/crab_CharmoniumRun2016F_21Feb2020_UL2016_HIPM-v1_AOD/results/processedLumis.json > lumi_charmonium_2016F_HIPM.txt
brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i CharmoniumRun2016UL/crab_CharmoniumRun2016G_21Feb2020_UL2016-v1_AOD/results/processedLumis.json > lumi_charmonium_2016G.txt
brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i CharmoniumRun2016UL/crab_CharmoniumRun2016H_21Feb2020_UL2016-v1_AOD/results/processedLumis.json > lumi_charmonium_2016H.txt

mv lumi_charmonium_2016* luminosity_brilcal/2016













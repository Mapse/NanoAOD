i=100

while [ $i -gt 0 ]
do

  crab status -d crab_projects_monte_carlo/crab_NANO_MC_jpsi_ccbar_3FS_4FS_SPS_2017_13TeV_13-08-2022 
  crab resubmit crab_projects_monte_carlo/crab_NANO_MC_jpsi_ccbar_3FS_4FS_SPS_2017_13TeV_13-08-2022 

  sleep 3600
  echo Hours left: $i
  ((i--))
done

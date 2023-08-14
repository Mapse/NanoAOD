i=1000

while [ $i -gt 0 ]
do
  #2017
  #crab resubmit CharmoniumRun2017UL/crab_CharmoniumRun2017F_09Aug2019_UL2017-v1_AOD  
  #crab resubmit CharmoniumRun2017UL/crab_CharmoniumRun2017B_09Aug2019_UL2017-v1_AOD 
  #2018
  #crab resubmit CharmoniumRun2018UL/crab_CharmoniumRun2018B_12Nov2019_UL2018-v1_AOD
  crab resubmit CharmoniumRun2018UL/crab_CharmoniumRun2018D_12Nov2019_UL2018-v1_AOD 

 sleep 400
  echo Hours left: $i
  ((i--))
done

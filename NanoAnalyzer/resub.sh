i=24

while [ $i -gt 0 ]
do
  crab resubmit CharmoniumRun2018UL/crab_CharmoniumRun2018A_AOD/
  crab resubmit CharmoniumRun2018UL/crab_CharmoniumRun2018B_AOD/
  crab resubmit CharmoniumRun2018UL/crab_CharmoniumRun2018C_AOD/
  crab resubmit CharmoniumRun2018UL/crab_CharmoniumRun2018D_AOD/

  sleep 3600
  echo Hours left: $i
  ((i--))
done

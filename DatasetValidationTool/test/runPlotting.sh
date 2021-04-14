#!/bin/bash

if [ "$#" == 0 ] || [ $1 == "-h" ] || [ $1 == "--help" ]
then
  printf "\nUsage: "
  printf "Plotter.sh [InputFileName1] [InputFileName2] \n\n"
  exit 1;
fi

INPUTFILE1=$1
INPUTFILE2=$2
CWD=`pwd`
FILE1=\"$CWD/$INPUTFILE1\"
FILE2=\"$CWD/$INPUTFILE2\"

if [ "$#" == 2 ]
then
   printf "Calling TestingGetKeys.C......"
   root -l -b -q "$CMSSW_BASE/src/DatasetValidation/DatasetValidationTool/test/TestingGetKeys.C(${FILE1},${FILE2})"
   printf "Calling complete.!"
fi

if [ "$#" -gt 2 ]
then
   root -l -b -q "$CMSSW_BASE/src/CosmicRateTool/TrackAnalyzer/macros/CosmicRateTool_CosmicRates.C(${FILE},$runMin,$runMax)"
fi

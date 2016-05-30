#!/usr/bin/env bash

################
# Instructions #
################

# 1. Copy GEN-Fragment.py in Configuration/Generator/python/
# 2. If not in repository copy name.dec in GeneratorInterface/ExternalDecays/data/
# 3. scram b
# 4. source GenMyMC.sh

if [ "$1" != "" ]; then
    targetPy="$1"
else
    targetPy="PYTHIA6_Bd2Psi2SKpi_TuneZ2star_8TeV_cff.py"
fi
globalTag="START53_V19F::All"
puScenario="NoPileUp" # NoPileUp, AVE_40_BX_25ns, CSA14_50ns_PoissonOOT, CAS14_inTime
numEvents=1000

echo targetPy="${targetPy}"

### Full steps: GEN,SIM,DIGI,L1,DIGI2RAW,HLT,RAW2DIGI,RECO
### First consalt PREP( http://cms.cern.ch/iCMS/prep/campaignmanagement?id=any ) to get official configuration

### Special step1: GEN+SIM for offline regional tracking ###
#cmsDriver.py Configuration/GenProduction/python/${targetPy} --step GEN,SIM --beamspot Realistic8TeVCollision --conditions START72_V1::All --pileup NoPileUp --datamix NODATAMIXER --eventcontent RAWSIM --datatier GEN-SIM --no_exec -n ${numEvents}

### Special step2: DIGI+L1+DIGI2RAW+PU mixing for offline regional tracking ###
#cmsDriver.py Configuration/GenProduction/python/${targetPy} --step DIGI,L1,DIGI2RAW,HLT:2014 --beamspot Realistic8TeVCollision --conditions START72_V1::All --pileup ${puScenario} --datamix NODATAMIXER --eventcontent RAWSIM --datatier GEN-SIM-RAW --no_exec -n ${numEvents}

### STEP0: GEN for acceptance measurement ###
#cmsDriver.py Configuration/GenProduction/python/${targetPy} --step GEN --beamspot Realistic8TeVCollision --conditions START72_V1::All --pileup NoPileUp --datamix NODATAMIXER --eventcontent RAWSIM --datatier GEN --no_exec -n ${numEvents}

### STEP1: GEN-SIM + DIGI + PU mixing + HLT + RECO (or AOD, --eventcontent AODSIM)###
cmsDriver.py Configuration/GenProduction/python/${targetPy} --step GEN,SIM,DIGI,L1,DIGI2RAW,HLT:7E33v2,RAW2DIGI,L1Reco,RECO --beamspot Realistic8TeVCollision --conditions ${globalTag} --pileup NoPileUp --datamix NODATAMIXER --eventcontent RECOSIM --datatier GEN-SIM-RECO --no_exec -n ${numEvents}

#cmsDriver.py Configuration/GenProduction/python/${targetPy} --step RAW2DIGI,L1Reco,RECO --conditions ${globalTag} --datamix NODATAMIXER --eventcontent RECOSIM --datatier GEN-SIM-RECO --no_exec -n ${numEvents}


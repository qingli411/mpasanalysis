#!/bin/bash

varname=$1
component="ocn"

case "${HOSTNAME}" in
    theta* )
        caseid=theta.20180906.branch_noCNT.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG
        drc_out=/lus/theta-fs0/projects/ClimateEnergy_3/qingli/e3sm_ts
        drc_in=/projects/ClimateEnergy_3/azamatm/E3SM_simulations/theta.20180906.branch_noCNT.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG/run
        drc_in2=/projects/ClimateEnergy_3/azamatm/E3SM_simulations/theta.20180906.branch_noCNT.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG/run-0006-01-01-180907--0046-01-01-190111
        e3sm_config=/lus/theta-fs0/projects/ccsm/acme/tools/e3sm-unified/load_latest_e3sm_unified_x.sh
        ts_ys=6
        ts_ye=55
        ;;
    edison* )
        caseid=edison.20181204.noCNT.A_WCYCL1950S_CMIP6_LRtunedHR.ne30_oECv3_ICG
        drc_out=/global/project/projectdirs/acme/qingli/e3sm_ts
        drc_in=/global/cscratch1/sd/tang30/ACME_simulations/edison.20181204.noCNT.A_WCYCL1950S_CMIP6_LRtunedHR.ne30_oECv3_ICG/archive/ocn/hist
        e3sm_config=/global/project/projectdirs/acme/software/anaconda_envs/load_latest_e3sm_unified_x.sh
        ts_ys=1
        ts_ye=50
        ;;
    blogin* )
        caseid=20190212.A_WCYCL1950S_CMIP6_LRtunedHR-noCNT.ne30_oECv3_ICG.anvil
        drc_in=/lcrc/group/acme/jwolfe/acme_scratch/anvil/20190212.A_WCYCL1950S_CMIP6_LRtunedHR-noCNT.ne30_oECv3_ICG.anvil/run
        drc_out=/lcrc/group/acme/qingli/e3sm_ts
        e3sm_config=/lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_x.sh
        ts_ys=1
        ts_ye=25
        ;;
    * )
        echo "This script should be executed on edison, theta or blues."
        exit 1
esac

# load nco
if [[ ! $(which ncclimo) ]]; then
    echo "Using e3sm_unified environment..."
    source ${e3sm_config}
fi

# variables
case ${varname} in
    "heatFlux" )
        varlist="timeMonthly_avg_sensibleHeatFlux,timeMonthly_avg_latentHeatFlux,timeMonthly_avg_shortWaveHeatFlux,timeMonthly_avg_longWaveHeatFluxUp,timeMonthly_avg_longWaveHeatFluxDown"
        ;;
    "freshWaterFlux" )
        varlist="timeMonthly_avg_evaporationFlux,timeMonthly_avg_rainFlux,timeMonthly_avg_snowFlux,timeMonthly_avg_seaIceSalinityFlux,timeMonthly_avg_seaIceFreshWaterFlux,timeMonthly_avg_riverRunoffFlux,timeMonthly_avg_iceRunoffFlux"
        ;;
    "windStress" )
        varlist="timeMonthly_avg_windStressZonal,timeMonthly_avg_windStressMeridional"
        ;;
    "mixedLayerDepth" )
        varlist="timeMonthly_avg_dThreshMLD,timeMonthly_avg_tThreshMLD"
        ;;
    "temperature" )
        varlist="timeMonthly_avg_activeTracers_temperature"
        ;;
    "salinity" )
        varlist="timeMonthly_avg_activeTracers_salinity"
        ;;
    "potentialDensity" )
        varlist="timeMonthly_avg_potentialDensity"
        ;;
    "velocity" )
        varlist="timeMonthly_avg_velocityZonal,timeMonthly_avg_velocityMeridional"
esac

ys_str=$(printf "%04d" ${ts_ys})
ye_str=$(printf "%04d" ${ts_ye})
path_out=${drc_out}/${caseid}/${ys_str}-${ye_str}/${component}/${varname}
mkdir -p ${path_out}

model="mpaso"

# Time series
ncclimo -p serial -v ${varlist} -c ${caseid} -m ${model} -s ${ts_ys} -e ${ts_ye} -o ${path_out} ${drc_in}/mpaso.hist.am.timeSeriesStatsMonthly.00??-??-01.nc

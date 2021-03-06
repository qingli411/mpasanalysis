#!/bin/bash

varname=$1
# variables
case ${varname} in
    "heatFlux" )
        varlist="timeMonthly_avg_sensibleHeatFlux,timeMonthly_avg_latentHeatFlux,timeMonthly_avg_shortWaveHeatFlux,timeMonthly_avg_longWaveHeatFluxUp,timeMonthly_avg_longWaveHeatFluxDown"
        component="ocn"
        ;;
    "freshWaterFlux" )
        varlist="timeMonthly_avg_evaporationFlux,timeMonthly_avg_rainFlux,timeMonthly_avg_snowFlux,timeMonthly_avg_seaIceSalinityFlux,timeMonthly_avg_seaIceFreshWaterFlux,timeMonthly_avg_riverRunoffFlux,timeMonthly_avg_iceRunoffFlux"
        component="ocn"
        ;;
    "windStress" )
        varlist="timeMonthly_avg_windStressZonal,timeMonthly_avg_windStressMeridional"
        component="ocn"
        ;;
    "mixedLayerDepth" )
        varlist="timeMonthly_avg_dThreshMLD,timeMonthly_avg_tThreshMLD"
        component="ocn"
        ;;
    "temperature" )
        varlist="timeMonthly_avg_activeTracers_temperature"
        component="ocn"
        ;;
    "salinity" )
        varlist="timeMonthly_avg_activeTracers_salinity"
        component="ocn"
        ;;
    "potentialDensity" )
        varlist="timeMonthly_avg_potentialDensity"
        component="ocn"
        ;;
    "velocity" )
        varlist="timeMonthly_avg_velocityZonal,timeMonthly_avg_velocityMeridional"
        component="ocn"
        ;;
    "iceFraction" )
        varlist="timeMonthly_avg_iceAreaCell"
        component="ice"
        ;;
    "iceThickness" )
        varlist="timeMonthly_avg_iceVolumeCell"
        component="ice"
        ;;
esac


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
        drc_in=/global/cscratch1/sd/tang30/ACME_simulations/edison.20181204.noCNT.A_WCYCL1950S_CMIP6_LRtunedHR.ne30_oECv3_ICG/archive/${component}/hist
        e3sm_config=/global/project/projectdirs/acme/software/anaconda_envs/load_latest_e3sm_unified_x.sh
        ts_ys=1
        ts_ye=50
        ;;
    blueslogin* )
        # caseid=20190212.A_WCYCL1950S_CMIP6_LRtunedHR-noCNT.ne30_oECv3_ICG.anvil
        # drc_in=/lcrc/group/acme/jwolfe/acme_scratch/anvil/20190212.A_WCYCL1950S_CMIP6_LRtunedHR-noCNT.ne30_oECv3_ICG.anvil/run
        # drc_out=/lcrc/group/acme/qingli/e3sm_ts
        # e3sm_config=/lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_x.sh
        # ts_ys=1
        # ts_ye=25
        caseid=20190419.test.A_WCYCL1850.ne30_oGNLD30to10.anvil
        drc_in=/lcrc/group/acme/jwolfe/acme_scratch/anvil/20190212.A_WCYCL1950S_CMIP6_LRtunedHR-noCNT.ne30_oECv3_ICG.anvil/run
        drc_in=/lcrc/group/acme/jwolfe/acme_scratch/anvil/20190419.test.A_WCYCL1850.ne30_oGNLD30to10.anvil/run
        drc_out=/lcrc/group/acme/qingli/e3sm_ts
        e3sm_config=/lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_x.sh
        ts_ys=1
        ts_ye=20
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

ys_str=$(printf "%04d" ${ts_ys})
ye_str=$(printf "%04d" ${ts_ye})
path_out=${drc_out}/${caseid}/${ys_str}-${ye_str}/${component}/${varname}
mkdir -p ${path_out}

# Time series
if [[ ${component} == "ocn" ]]; then
    model="mpaso"
elif [[ ${component} == "ice" ]]; then
    model="mpascice"
fi
ncclimo -p serial -v ${varlist} -c ${caseid} -m ${model} -s ${ts_ys} -e ${ts_ye} -o ${path_out} ${drc_in}/${model}.hist.am.timeSeriesStatsMonthly.00??-??-01.nc


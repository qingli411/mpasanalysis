#!/bin/bash

varname=$1
component="ocn"

case "${HOSTNAME}" in
    theta* )
        caseid=theta.20180906.branch_noCNT.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG
        drc_out=/lus/theta-fs0/projects/ClimateEnergy_3/qingli/e3sm_climo
        drc_in=/projects/ClimateEnergy_3/azamatm/E3SM_simulations/theta.20180906.branch_noCNT.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG/run
        drc_in2=/projects/ClimateEnergy_3/azamatm/E3SM_simulations/theta.20180906.branch_noCNT.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG/run-0006-01-01-180907--0046-01-01-190111
        e3sm_config=/lus/theta-fs0/projects/ccsm/acme/tools/e3sm-unified/load_latest_e3sm_unified_x.sh
        climo_ys=46
        climo_ye=55
        ;;
    edison* )
        caseid=edison.20181204.noCNT.A_WCYCL1950S_CMIP6_LRtunedHR.ne30_oECv3_ICG
        drc_out=/global/project/projectdirs/acme/qingli/e3sm_climo
        drc_in=/global/cscratch1/sd/tang30/ACME_simulations/edison.20181204.noCNT.A_WCYCL1950S_CMIP6_LRtunedHR.ne30_oECv3_ICG/archive/ocn/hist
        e3sm_config=/global/project/projectdirs/acme/software/anaconda_envs/load_latest_e3sm_unified_x.sh
        climo_ys=41
        climo_ye=50
        ;;
    * )
        echo "This script should be executed on either edison or theta."
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

ys_str=$(printf "%04d" ${climo_ys})
ye_str=$(printf "%04d" ${climo_ye})
path_out=${drc_out}/${caseid}/${ys_str}-${ye_str}/${component}/${varname}
mkdir -p ${path_out}

model="mpaso"

# climatology
ncclimo -p serial -v ${varlist} -c ${caseid} -m ${model} -s ${climo_ys} -e ${climo_ye} --seasons=none --dec_md=sdd -i ${drc_in} -o ${path_out}

# Time series
# ncclimo -p serial -v ${varlist} -c ${caseid} -m ${model} -s 6 -e 46 -o $drc_out $drc_in/mpaso.hist.am.timeSeriesStatsMonthly.00??-0?-01.nc

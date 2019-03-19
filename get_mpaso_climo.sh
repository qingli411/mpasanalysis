#!/bin/bash

varname=$1
climo_ys=$2
climo_ye=$3
component="ocn"

case "${HOSTNAME}" in
    theta* )
        caseid=theta.20180906.branch_noCNT.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG
        drc_out=/lus/theta-fs0/projects/ClimateEnergy_3/qingli/e3sm_climo
        drc_in=/projects/ClimateEnergy_3/azamatm/E3SM_simulations/theta.20180906.branch_noCNT.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG/run
        drc_in2=/projects/ClimateEnergy_3/azamatm/E3SM_simulations/theta.20180906.branch_noCNT.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG/run-0006-01-01-180907--0046-01-01-190111
        e3sm_config=/lus/theta-fs0/projects/ccsm/acme/tools/e3sm-unified/load_latest_e3sm_unified_x.sh
        ;;
    edison* )
        caseid=edison.20181204.noCNT.A_WCYCL1950S_CMIP6_LRtunedHR.ne30_oECv3_ICG
        drc_out=/global/project/projectdirs/acme/qingli/e3sm_climo
        drc_in=/global/cscratch1/sd/tang30/ACME_simulations/edison.20181204.noCNT.A_WCYCL1950S_CMIP6_LRtunedHR.ne30_oECv3_ICG/archive/ocn/hist
        e3sm_config=/global/project/projectdirs/acme/software/anaconda_envs/load_latest_e3sm_unified_x.sh
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
if [[ ${climo_ys} == ${climo_ye} ]]; then
    echo "Single year output of ${varname} for year ${climo_ys}."
    yyyy=$(printf %04d ${climo_ys})
    for i in {1..12}; do
        mm=$(printf %02d ${i})
        ncks -O -v ${varlist} ${drc_in}/mpaso.hist.am.timeSeriesStatsMonthly.${yyyy}-${mm}-01.nc ${path_out}/mpaso_${mm}_${yyyy}${mm}_${yyyy}${mm}_climo.nc
        ln -sf ${path_out}/mpaso_${mm}_${yyyy}${mm}_${yyyy}${mm}_climo.nc ${path_out}/mpaso_${mm}_climo.nc
    done
else
    ncclimo -p serial -v ${varlist} -c ${caseid} -m ${model} -s ${climo_ys} -e ${climo_ye} --seasons=none --dec_md=sdd -i ${drc_in} -o ${path_out}
fi

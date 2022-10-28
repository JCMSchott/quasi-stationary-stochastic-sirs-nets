#!/bin/bash

##########################################################

amostra=${1}

num_base=1

tam=${num_base}000; grau_min=3; gama=5.5;
##########################################################
# dlambda = 0.0125/divisor

qtidade=1.0

#n18g35a0p2_7,7,10000000,3,3.5,0.0,20.0,1.5,0.2,200000,600000,false,0,0,0.05

divisor=10.0; # para n = 1e5, 15.0 para rodar rapido, mas sera 17.5 na criticalidade

lambda0=0.0;

lambdaf=1.5;

interpola=false;

ind_interp=0
                 # -1: aceita o lamb0 dado aqui
                 # 0: apenas continua do ultimo lambda lido
                 # 1: dlamb/4; 2: dlamb/3;
                 # 3: dlamb/2; 4: 2 * dlamb/3
                 # 5: 3 * dlamb/4
                 # 6: lamb0 = lamb0

ind_soCalculaIPR=0

##########################################################
alp=1.0;
##########################################################
tRelaxacao=500000; tempoTotal=1500000
##########################################################

ind_no=( 0 1 2 3 5 8 12 14 15 16 )

##########################################################
if ${interpola}; then
   lambda0=$( echo "scale=20; ${lambda0} + 0.0125/(2 * ${divisor}) " | bc -l )
   echo "Interpolara" 
else
#   lambda0=$( echo "scale=20; ${lambda0} + 0.0125/${divisor} " | bc -l )
   echo "Nao interpolara"
#   echo "O valor de lambda0 eh  ${lambda0}"
fi
##########################################################

#flags='-check all -traceback'

#flags='-heap-arrays -O3 -g -fp-stack-check'

#flags='-traceback -check all -heap-arrays -O3 -g -fp-stack-check'

flags=''

##########################################################
# Nao mexer daqui pra baixo!
nzeros=$( echo "scale=0; l(${tam}/${num_base})/l(10) " | bc -l)

IFS='.'

arr=( ${gama} )

gama_sp=${arr[0]}${arr[1]}

alp_arr=( ${alp} )

IFS=' '
##########################################################
dependencias='mod_rndgen_multiple.f90 geraRede.f90 mod_tools_redes.f90'
principal=' main_SIRS_Estocastico.f90'
executavel='n'${num_base}${nzeros}g${gama_sp}'a'${alp_arr[0]}p${alp_arr[1]}'_'${amostra}
#############################################################################
rm ${executavel} &> erroRmExe.log

ifort ${dependencias} ${principal} ${flags} -o ${executavel}

rm -r *.mod

#time ./${executavel} ${amostra} ${tam} ${grau_min} ${gama} ${lambda0} ${divisor} ${lambdaf} ${alp} ${tRelaxacao} ${tempoTotal} ${interpola} ${ind_interp} ${ind_soCalculaIPR} ${qtidade}

qsub -N ${executavel} -cwd executa.sh ${executavel} ${amostra} ${tam} ${grau_min} ${gama} ${lambda0} ${divisor} ${lambdaf} ${alp} ${tRelaxacao} ${tempoTotal} ${interpola} ${ind_interp} ${ind_soCalculaIPR} ${qtidade}

#qsub -N ${executavel} -q all.q@compute-0-${ind_no[${amostra}-1]} -cwd executa.sh ${executavel} ${amostra} ${tam} ${grau_min} ${gama} ${lambda0} ${divisor} ${lambdaf} ${alp} ${tRelaxacao} ${tempoTotal} ${interpola} ${ind_interp} ${ind_soCalculaIPR} ${qtidade}

#$ -S /bin/sh


programa=${1}
tam_rede=${2}
grau_rrn=${3}
lamb0=${4}
divisor=${5}
alp=${6}
n_lamb=${7}
tTotal=${8}
tRelax=${9}
interpola=${10}
ind_ams=${11}

echo $@ > argumentos.log

time ./${programa} ${tam_rede} ${grau_rrn} ${lamb0} ${divisor} ${alp} ${n_lamb} ${tTotal} ${tRelax} ${interpola} ${ind_ams}

echo ${tam_rede} ${grau_rrn} ${lamb0} ${divisor} ${alp} ${n_lamb} ${tTotal} ${tRelax} ${interpola} ${ind_ams} > argumentos.log

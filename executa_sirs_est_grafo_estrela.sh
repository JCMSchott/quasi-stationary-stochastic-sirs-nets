#$ -S /bin/sh


programa=${1}
tam_rede=${2}
lamb0=${3}
divisor=${4}
alp=${5}
n_lamb=${6}
tTotal=${7}
tRelax=${8}
interpola=${9}

echo $@ > argumentos.log

time ./${programa} ${tam_rede} ${lamb0} ${divisor} ${alp} ${n_lamb} ${tTotal} ${tRelax} ${interpola}

echo ${tam_rede} ${lamb0} ${divisor} ${alp} ${n_lamb} ${tTotal} ${tRelax} ${interpola} > argumentos.log

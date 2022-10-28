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
grauStar=${11}
ind_lamb=${12}
Eh_Limiar=${13}

echo $@ > argumentos.log

time ./${programa} ${tam_rede} ${grau_rrn} ${lamb0} ${divisor} ${alp} ${n_lamb} ${tTotal} ${tRelax} ${interpola} ${grauStar} ${ind_lamb} ${Eh_Limiar}

echo ${tam_rede} ${grau_rrn} ${lamb0} ${divisor} ${alp} ${n_lamb} ${tTotal} ${tRelax} ${interpola} ${grauStar} ${ind_lamb} ${Eh_Limiar} > argumentos.log


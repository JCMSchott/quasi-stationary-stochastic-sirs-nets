#$ -S /bin/sh


programa=${1};

amostra=${2}; tam=${3}; grau_min=${4}; gama=${5};

lambda0=${6}; divisor=${7}; lambdaf=${8};

alp=${9};

tRelaxacao=${10}; tempoTotal=${11}

interpola=${12}

ind_interp=${13}

ind_soCalculaIPR=${14}

qtidade=${15}

echo $@ > argumentos.log

time ./${programa} ${amostra} ${tam} ${grau_min} ${gama} ${lambda0} ${divisor} ${lambdaf} ${alp} ${tRelaxacao} ${tempoTotal} ${interpola} ${ind_interp} ${ind_soCalculaIPR} ${qtidade}


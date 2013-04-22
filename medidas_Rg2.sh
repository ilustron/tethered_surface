#!/bin/bash

mkdir -p ./MEDIDAS_Rg2

for l in 16 32 46 64
do
    mkdir -p ./MEDIDAS_Rg2/L$l #Crea el directorio correspondiente al tamaño L de la membrana
    echo Datos necesarios para las medidas del radio del giratón para L$l:

    echo Número de archivos:
    read f

    echo Periodo de correlación estimado:
    read tau

    echo Tiempo de correlación estimado:
    read termal

    for i in `seq 5 11`;
    do 
	k=`echo "$i/10" | bc -l | awk '{printf("%2.1f", $1);}'`
	mkdir -p ./MEDIDAS_Rg2/L$l/K$k
	echo Medidas para K=$k:
	gcc -O2 -DK=$k -DL=$l -DNF=$f -DTAU=$tau -DNTERMAL=$termal medida_Rg2.c -lm -o medida_Rg2.out
	./medida_Rg2.out
    done
    
    k=2.0
    mkdir -p ./MEDIDAS_Rg2/L$l/K$k
    echo Medidas para K=$k:
    gcc -O2 -DK=$k -DL=$l -DNF=$f -DTAU=$tau -DNTERMAL=$termal medida_Rg2.c -lm -o medida_Rg2.out
    ./medida_Rg2.out
done

#!/bin/bash

mkdir ./MEDIDAS_Se

for l in 16 32 46 64
do
    mkdir -p ./MEDIDAS_Se/L$l #Crea el directorio correspondiente al tamaño L de la membrana
    echo Datos necesarios para las medidas de la energía de curvatura para L$l:

    echo Número de archivos:
    read f

    for i in `seq 5 11`;
    do 
	k=`echo "$i/10" | bc -l | awk '{printf("%2.1f", $1);}'`
	mkdir -p ./MEDIDAS_Se/L$l/K$k
	echo Medidas para K=$k:
	gcc -O2 -DK=$k -DL=$l -DNF=$f medida_Se.c -lm -o medida_Se.out
	./medida_Se.out
	echo Convirtiendo Phistorico_Se-L$l-K$i.eps a .pdf
	epstopdf ./MEDIDAS_Se/L$l/K$k/Phistorico_Se-L$l-K$i.eps
	echo Hecho
	echo Convirtiendo Perror_Se-L$l-K$i.eps a .pdf
	epstopdf ./MEDIDAS_Se/L$l/K$k/Perror_Se-L$l-K$i.eps
	echo Hecho
	echo Convirtiendo Plotermal_Se-L$l-K$i.eps a .pdf
	epstopdf ./MEDIDAS_Se/L$l/K$k/Plogtermal_Se-L$l-K$i.eps
	echo Hecho
    done
    
    k=2.0
    mkdir -p ./MEDIDAS_Se/L$l/K$k
    echo Medidas para K=$k:
    gcc -O2 -DK=$k -DL=$l -DNF=$f medida_Se.c -lm -o medida_Se.out
    ./medida_Se.out
    echo Convirtiendo Phistorico_Se-L$l-K$i.eps a .pdf
    epstopdf ./MEDIDAS_Se/L$l/K$k/Phistorico_Se-L$l-K20.eps
    echo Hecho
    echo Convirtiendo Ploterror_Se-L$l-K20.eps a .pdf
    epstopdf ./MEDIDAS_Se/L$l/K$k/Perror_Se-L$l-K20.eps
    echo Hecho
    echo Convirtiendo Plotermal_Se-L$l-K20.eps a .pdf
    epstopdf ./MEDIDAS_Se/L$l/K$k/Plogtermal_Se-L$l-K20.eps
    echo Hecho
done

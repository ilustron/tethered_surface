#!/bin/bash

mkdir ./MEDIDAS_Drg2 #Crea el directorio donde se almacenarán la medidas y gráficos del Cv 

for l in 16 32 46 64 128 #Tamaños de la membrana donde se realizan las medidas
do
    mkdir -p ./MEDIDAS_Drg2/L$l #Crea el directorio correspondiente al tamaño L de la membrana
    echo Datos necesarios para las medidas del calor específico para L$l:

    cat ./RUNS/L$l/K0.5/regparams.log #Muestra el log de k0.5
    
    echo Número de archivos:
    read f

    for i in `seq 5 11`;
    do 
	k=`echo "$i/10" | bc -l | awk '{printf("%2.1f", $1);}'`
	mkdir -p ./MEDIDAS_Drg2/L$l/K$k
	echo Medidas para K=$k:
	gcc -O2 -DK=$k -DL=$l -DNF=$f medida_Drg2.c -lm -o medida_Drg2.out
	./medida_Drg2.out
        #Pasar a pdf los eps resultantes
	echo Convirtiendo Phistorico_Drg2-L$l-K$i.eps a .pdf
	epstopdf ./MEDIDAS_Drg2/L$l/K$k/Phistorico_Drg2-L$l-K$i.eps
	echo Hecho
	echo Convirtiendo Pterror_Drg2-L$l-K$i.eps a .pdf
	epstopdf ./MEDIDAS_Drg2/L$l/K$k/Perror_Drg2-L$l-K$i.eps
	echo Hecho
	echo Convirtiendo Plotermal_Drg2-L$l-K$i.eps a .pdf
	epstopdf ./MEDIDAS_Drg2/L$l/K$k/Plogtermal_Drg2-L$l-K$i.eps
	echo Hecho
    done
    
    k=2.0
    mkdir -p ./MEDIDAS_Drg2/L$l/K$k
    echo Medidas para K=$k:
    gcc -O2 -DK=$k -DL=$l -DNF=$f medida_Drg2.c -lm -o medida_Drg2.out
    ./medida_Drg2.out
    #Pasar a pdf los eps resultantes
    echo Convirtiendo Phistorico_Drg2-L$l-K$i.eps a .pdf
    epstopdf ./MEDIDAS_Drg2/L$l/K$k/Phistorico_Drg2-L$l-K20.eps
    echo Hecho
    echo Convirtiendo Perror_Drg2-L$l-K20.eps a .pdf
    epstopdf ./MEDIDAS_Drg2/L$l/K$k/Perror_Drg2-L$l-K20.eps
    echo Hecho
    echo Convirtiendo Plotermal_Drg2-L$l-K20.eps a .pdf
    epstopdf ./MEDIDAS_Drg2/L$l/K$k/Plogtermal_Drg2-L$l-K20.eps
    echo Hecho
done

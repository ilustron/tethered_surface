#!/bin/bash

mkdir -p ./MEDIDAS_Rg2

for l in 16 32 46 64
do
    mkdir -p ./MEDIDAS_Rg2/L$l #Crea el directorio correspondiente al tamaño L de la membrana
    echo Datos necesarios para las medidas del radio del giratón para L$l:

    cat ./RUNS/L$l/K0.5/regparams.log #Muestra el log de k0.5

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
       #Pasa a pdf los eps resultantes
	epstopdf ./MEDIDAS_Rg2/L$l/K$k/Ploterror_Rg2-L$l-K$i.eps
	epstopdf ./MEDIDAS_Rg2/L$l/K$k/Plottermal_Rg2-L$l-K$i.eps
    done
    
    k=2.0
    mkdir -p ./MEDIDAS_Rg2/L$l/K$k
    echo Medidas para K=$k:
    gcc -O2 -DK=$k -DL=$l -DNF=$f -DTAU=$tau -DNTERMAL=$termal medida_Rg2.c -lm -o medida_Rg2.out
    ./medida_Rg2.out
    #Pasar a pdf los eps resultantes
    epstopdf ./MEDIDAS_Rg2/L$l/K$k/Ploterror_Rg2-L$l-K20.eps
    epstopdf ./MEDIDAS_Rg2/L$l/K$k/Plottermal_Rg2-L$lK20.eps
done

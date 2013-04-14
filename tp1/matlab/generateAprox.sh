#! /bin/bash

paramfile=""
            if [ $# -eq 0 ]
            then
paramfile="../informe/plots/parameters.txt"
            else
paramfile="../informe/plots/parameters-eps.txt"
            fi



files=(../datos/x_15_9_3.txt ../datos/x_2_62_35.txt ../datos/X1.txt ../datos/X2.txt ../datos/X3.txt ../datos/X4.txt ../datos/X5.txt ../datos/X6.txt ../datos/X7.txt)
if [ $# -eq 0 ] 
then
    variables=(Beta Lambda Sigma)
else
    variables=(Beta Iteraciones)
fi

betas=(9 6.2)
lambdas=(3 3.5)
sigmas=(1.5 2)

index=0
caso=1

metodos=(Newton RegulaFalsi)
nombreposta=(Newton Illinois)
metval=(0 1)

precision=(15 16 17 18 19 20 21 22 23 24 25 26 27)
epsilons=(0.01 0.001 0.0001 0.00001 0.000001 0.0000001 0.00000001 0.000000001)
for file in "${files[@]}"; do
    for met in "${metodos[@]}"; do

        params=""
        for var in "${variables[@]}"; do

            a=$(grep "caso$caso" "$paramfile"| grep $met | grep "$var" | awk '{print $5}' | xargs)
            val=. read -a serie <<< "$a"
            index=0
            if [ $# -eq 0 ]
            then
            for prec in "${precision[@]}"; do
                params+=" ${serie[$index]} "
                index=$[$index +1]
            done
            else
                for prec in "${epsilons[@]}"; do
                params+=" ${serie[$index]} "
                index=$[$index +1]
            done
 
            fi
            if [ $# -eq 0 ]
            then
                filename="Caso$caso-$met"
            else
                filename="Caso$caso-$met-eps"
            fi


        done
            if [ $# -eq 0 ]
            then
                echo "./plotCurve.m ${precision[0]} $params $filename"
               ./plotCurve.m ${precision[0]} $params $filename
            else
                 echo "./plotCurveSoloBeta.m 2 $params $filename"
               ./plotCurveSoloBeta.m 2 $params $filename
            fi
        mv "$filename.png" ../informe/plots/$filename.png
    done
	caso=$[$caso +1]
done



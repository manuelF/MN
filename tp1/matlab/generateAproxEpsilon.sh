#! /bin/bash

paramfile=""
	paramfile="../informe/plots/parameters-eps.txt"



files=(../datos/x_15_9_3.txt ../datos/x_2_62_35.txt ../datos/X1.txt ../datos/X2.txt ../datos/X3.txt ../datos/X4.txt ../datos/X5.txt ../datos/X6.txt ../datos/X7.txt)
variables=(Beta Iteraciones)

betas=(9 6.2)
lambdas=(3 3.5)
sigmas=(1.5 2)

index=0
caso=1

metodos=(Newton RegulaFalsi)
nombreposta=(Newton Illinois)
metval=(0 1)

precision=(21 50)
epsilons=(0.01 0.001 0.0001 0.00001 0.000001 0.0000001 0.00000001 0.000000001)
for file in "${files[@]}"; do
	for met in "${metodos[@]}"; do

		for pre in "${precision[@]}"; do

			params=""
			for var in "${variables[@]}"; do

				a=$(grep "caso$caso" "$paramfile"| grep $met | grep "$var" | grep "\\-$pre"|  awk '{print $5}' | xargs)
				val=. read -a serie <<< "$a"
				index=0
				for prec in "${epsilons[@]}"; do
					params+=" ${serie[$index]} "
					index=$[$index +1]
				done



			done
		    filename="Caso$caso-$met-$pre-eps"
			echo "./plotCurveSoloBeta.m 2 $params $filename"
				./plotCurveSoloBeta.m 2 $params $filename
			mv "$filename.png" ../informe/plots/$filename.png
		done
	done
	caso=$[$caso +1]
done


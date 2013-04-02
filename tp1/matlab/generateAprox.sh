#! /bin/bash


paramfile="../informe/plots/parameters.txt"

files=(../datos/x_15_9_3.txt ../datos/x_2_62_35.txt ../datos/X1.txt ../datos/X2.txt ../datos/X3.txt ../datos/X4.txt ../datos/X5.txt ../datos/X6.txt ../datos/X7.txt)
variables=(Beta Lambda Sigma)

betas=(9 6.2)
lambdas=(3 3.5)
sigmas=(1.5 2)

index=0
caso=1

precision=(12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27)

for file in "${files[@]}"; do
	for var in "${variables[@]}"; do

		a=$(grep "caso$caso" "$paramfile" | grep "$var" | awk '{print $5}' | xargs)
		val=. read -a serie <<< "$a"
		params=""
		index=0
		for prec in "${precision[@]}"; do
			params+=" ${serie[$index]} "
			index=$[$index +1]
		done

		echo "./plot_curve ${precision[0]} $params"
		filename="Caso$caso-$var"
		./plotCurve.m ${precision[0]} $params $filename
		mv "$filename.png" ../informe/plots

	done
	caso=$[$caso +1]
done



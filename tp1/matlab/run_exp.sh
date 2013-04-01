#! /bin/bash
files=(../datos/x_15_9_3.txt ../datos/x_2_62_35.txt ../datos/X1.txt ../datos/X2.txt ../datos/X3.txt ../datos/X4.txt ../datos/X5.txt ../datos/X6.txt ../datos/X7.txt)
betas=(9 6.2)
lambdas=(3 3.5)
sigmas=(1.5 2)

metodos=(Newton PF)
metodosval=(0 1)

precision=(17 35 50)
index=0
caso=1

paramfile="../informe/plots/parameters.txt"
rm $paramfile
touch $paramfile

for file in "${files[@]}"; do
	truebeta=${betas[$index]}
	truelambda=${lambdas[$index]}
	truesigma=${sigmas[$index]}
	echo "$file : " >> $paramfile
	echo "" >> $paramfile
	for met in "${metodosval[@]}"; do

		for prec in "${precision[@]}"; do
			#correr con prec, 100 iteraciones, beta entre 20 y 0
			a=$(../src/Newton  "${prec}" 100 20.0 0 ${metodosval} < "${file}")
			val=. read -a vals <<< "$a"
			beta=${vals[0]}
			lambda=${vals[1]}
			sigma=${vals[2]}
			its=${vals[3]}
			filename="${metodos[$metodosval]}-${prec}-caso$[$caso].png"
			echo "./histogramaYAjuste.m ${file} $sigma $beta $lambda $truesigma $truebeta $truelambda $filename"
			echo "$filename : Beta = $beta" >> $paramfile
			echo "$filename : Lambda = $lambda" >> $paramfile
			echo "$filename : Sigma = $sigma" >> $paramfile
			echo "$filename : Iteraciones = $its" >> $paramfile
			echo "" >> $paramfile
			./histogramaYAjuste.m ${file} $truesigma $truebeta $truelambda $sigma $beta $lambda $filename
			mv $filename ../informe/plots/
		done
	done
	index=$[$index +1]
	caso=$[$caso +1]
done
#! /bin/bash
files=(../datos/x_15_9_3.txt ../datos/x_2_62_35.txt)
betas=(9 6.2)
lambdas=(3 3.5)
sigmas=(1.5 2)

precision=(17 35 50)
index=0
caso=1
for file in "${files[@]}"; do
	truebeta=${betas[$index]}
	truelambda=${lambdas[$index]}
	truesigma=${sigmas[$index]}

	for prec in "${precision[@]}"; do
		#correr con prec, 100 iteraciones, beta entre 20 y 0
		a=$(../src/Newton  "${prec}" 100 20.0 0 < "${file}")
		val=. read -a vals <<< "$a"
		beta=${vals[0]}
		lambda=${vals[1]}
		sigma=${vals[2]}
		filename="Newton-${prec}-caso$[$caso].png"
		echo "./histogramaYAjuste.m ${file} $sigma $beta $lambda $truesigma $truebeta $truelambda $filename"
		./histogramaYAjuste.m ${file} $sigma $beta $lambda $truesigma $truebeta $truelambda $filename
		mv $filename ../informe/plots/
	done
	index=$[$index +1]
	caso=$[$caso +1]
done

#for prec in "${precision[@]}"; do
	#correr con prec, 100 iteraciones, beta entre 20 y 0
#	../src/Newton  "${prec}" 100 20.0 0 < ../datos/x_2_62_35.txt
#done

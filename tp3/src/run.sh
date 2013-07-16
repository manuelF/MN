#! /bin/bash
make entero;
./OCR 100 0 2 5000 > HR_5k_0_2 
./OCR 100 1 2 5000 > HR_5k_1_2 
rm -f V.txt
./OCR 100 0 2 15000 > HR_15k_0_2 
./OCR 100 1 2 15000 > HR_15k_1_2 
rm -f V.txt
./OCR 100 0 2 30000 > HR_30k_0_2  
./OCR 100 1 2 30000 > HR_30k_1_2 



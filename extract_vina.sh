#!/bin/bash
#This command helps extract vina output results from a .txt file previously merged with a cat command

cat vina_output_1.txt vina_output_2.txt vina_output_3.txt vina_output_4.txt vina_output_5.txt vina_output_6.txt vina_output_7.txt vina_output_8.txt vina_output_9.txt vina_output_10.txt > ifit1_docking_results.txt

#or use cat with the *.txt wildcard in a directory which contains all .txt outputs

awk 'NR>3 && ($1==1 || $1==3 || $1==2 || $1==4 || $1==5 || $1==6 || $1==7 || $1==8 || $1==9)' ifit1_docking_results.txt > ifit1_extracted.txt


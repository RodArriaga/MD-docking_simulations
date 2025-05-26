#!/bin/bash

#This script performs 10 vina runs and outputs the energy results in a .txt file, and the pdb results in a .pbqt file. Run this in a dedicated directory for each ligand and receptor changing the config file to the corresponding receptor-ligand docking procedure.
#Copy the receptor, ligand (.pdbqt) and config file (.txt) into the directory where the script will be run

for i in {1..10}
do
OUTFILE="vina_output_${i}.txt"
PDBQTFILE="vina_output_${i}.pdbqt"
/Users/adrianarriaga/Documents/RNA_Genomics/autodock_vina/bin/vina --config /Users/adrianarriaga/Documents/RNA_Genomics/docking_cap0/IFI1_cap0AA_config.txt --out "$PDBQTFILE" > "$OUTFILE"
echo "Run $i complete. Energy output written to $OUTFILE and PDBQT output written to $PDBQTFILE"
done

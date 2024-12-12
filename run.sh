#!/bin/bash

# Directorio de herramientas paralelas
partools="$(pwd)/../parallel"

# Archivo de entrada
ftsv=$1

# Validar archivo de entrada
[[ $ftsv == "" ]] && echo "need .tsv file" && exit 15

# Procesar archivo .tsv
echo "Processing MI calculations for: $ftsv ..."
nom=$(echo $ftsv | cut -d. -f 1)

awk '{print $1}' $ftsv > node.list
cname=$(head -1 node.list)
echo "Column Index Name: $cname"

# Calcular MI con 4 núcleos
SECONDS=0
python ${partools}/aracne-par.py $ftsv node.list $cname 4 &> aracne.log
echo "ARACNe time: $(echo $SECONDS/60 | bc -l) minutes."

# Crear archivo SIF directamente desde los cálculos
SECONDS=0
python ${partools}/adj2sif.py > ${nom}.sif
echo "Creating SIF: $(echo $SECONDS/60 | bc -l) minutes."

# Ordenar y preparar salida
SECONDS=0
sort -r -k3,3 ${nom}.sif > ${nom}.sort
echo "Sorting: $(echo $SECONDS/60 | bc -l) minutes."

# Limpiar archivos temporales
rm -rf adj log node.list

echo "Output files are ready."

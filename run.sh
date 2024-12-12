#!/bin/bash

# Ruta a la carpeta de herramientas
partools="$(pwd)/../parallel"

# Archivo .tsv como argumento
ftsv=$1

# Número de núcleos (por defecto 4)
num_cores=4

# Verificación de argumento
[[ $ftsv == "" ]] && echo "need .tsv file" && exit 15

# Procesamiento del archivo .tsv
echo "Processing MI calculations for: $ftsv ..."
nom=$(echo $ftsv | cut -d. -f 1)

# Creación de la lista de nodos
awk '{print $1}' $ftsv > node.list
cname=$(head -1 node.list)
echo "Column Index Name: $cname"

# Cálculo de MI con ARACNe
SECONDS=0
python ${partools}/aracne-par.py $ftsv node.list $cname $num_cores &> aracne.log 
echo "ARACNe time: $(echo $SECONDS/60 | bc -l) minutes."

# Unión de matrices de adyacencia
SECONDS=0
n=$( (cd adj; ls) | head -1 | cut -d'.' -f 2 )
echo "Parameters to join: $nom $n node.list $cname"
python ${partools}/joinadj.py $nom $n node.list $cname
echo "join ADJ matriz time: $(echo $SECONDS/60 | bc -l) minutes."

# Moviendo matriz de adyacencia
echo "Moving adjacency matrix"
mv adj/mat.adj .

# Creación del archivo SIF
SECONDS=0
python ${partools}/adj2sif.py > ${nom}.sif
echo "Creating SIF: $(echo $SECONDS/60 | bc -l) minutes."

# Ordenamiento del archivo SIF
SECONDS=0
sort -r -k3,3 ${nom}.sif > ${nom}.sort
echo "Sorting: $(echo $SECONDS/60 | bc -l) minutes."

# Limpieza de archivos temporales
rm -rf adj log mat.adj node.list

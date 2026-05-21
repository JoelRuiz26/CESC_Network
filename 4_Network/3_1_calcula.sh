#!/bin/bash
# Uso: bash calcula.sh <NombreDeTuArchivo.tsv>
bash run.sh "$1" &> salida &
echo "Ejecutando: bash run.sh $1 &> salida &"

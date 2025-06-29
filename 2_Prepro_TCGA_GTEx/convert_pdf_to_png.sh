#!/bin/bash

echo "¿En qué carpeta están los PDFs que quieres convertir?"
read -e -p "Ruta de la carpeta: " pdf_dir

# Si el usuario da ruta relativa, la convertimos a absoluta
pdf_dir=$(realpath "$pdf_dir")

if [ ! -d "$pdf_dir" ]; then
    echo "¡Error! La carpeta '$pdf_dir' no existe."
    exit 1
fi

cd "$pdf_dir"

echo "Convirtiendo PDFs en $pdf_dir a PNGs (1200dpi)..."
for pdf in *.pdf; do
    # Verifica que haya archivos PDF
    if [ ! -e "$pdf" ]; then
        echo "No hay PDFs en $pdf_dir"
        exit 0
    fi
    base=$(basename "$pdf" .pdf)
    pdftoppm -png -rx 1200 -ry 1200 "$pdf" "$base"
    # Renombra si solo hay una página
    if [ -f "${base}-1.png" ] && [ ! -f "${base}.png" ]; then
        mv "${base}-1.png" "${base}.png"
    fi
done

echo "¡Listo! PNGs a 1200dpi guardados en $pdf_dir"

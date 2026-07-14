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

# Subcarpetas de salida: antes / después / (otros sin etiqueta)
before_dir="$pdf_dir/png_before"
after_dir="$pdf_dir/png_after"
other_dir="$pdf_dir/png_otros"
mkdir -p "$before_dir" "$after_dir"

echo "Convirtiendo PDFs en $pdf_dir a PNGs (1200dpi)..."
for pdf in *.pdf; do
    # Verifica que haya archivos PDF
    if [ ! -e "$pdf" ]; then
        echo "No hay PDFs en $pdf_dir"
        exit 0
    fi
    base=$(basename "$pdf" .pdf)

    # Elige subcarpeta según el nombre del archivo (before/after)
    case "$base" in
        *_before|*before*) out_dir="$before_dir" ;;
        *_after|*after*)   out_dir="$after_dir" ;;
        *)                 out_dir="$other_dir"; mkdir -p "$other_dir" ;;
    esac

    pdftoppm -png -rx 1200 -ry 1200 "$pdf" "$out_dir/$base"
    # Renombra si solo hay una página
    if [ -f "$out_dir/${base}-1.png" ] && [ ! -f "$out_dir/${base}.png" ]; then
        mv "$out_dir/${base}-1.png" "$out_dir/${base}.png"
    fi
done

echo "¡Listo! PNGs a 1200dpi:"
echo "  - antes:    $before_dir"
echo "  - después:  $after_dir"
[ -d "$other_dir" ] && echo "  - otros:    $other_dir"

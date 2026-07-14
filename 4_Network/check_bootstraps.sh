#!/bin/bash
# check_bootstraps.sh — audita la completitud de los archivos bootstrapNetwork_*.txt de ARACNe.
#
# Metodo (basado en Aracne.java::writeFinal):
#   ARACNe escribe cada bootstrap con un BufferedWriter recorriendo TODOS los reguladores
#   en un orden de HashMap fijo (identico entre bootstraps, misma lista de genes) y cierra al final.
#   Por tanto, un bootstrap COMPLETO cumple 3 invariantes:
#     (1) contiene el numero total de reguladores distintos (= max observado),
#     (2) su ultima linea es el ultimo regulador del orden de escritura (el mismo en todos),
#     (3) el archivo termina en salto de linea.
#   Un archivo truncado (proceso muerto durante writeFinal) falla (1) y/o (2).
#
# Uso: bash check_bootstraps.sh <directorio_con_bootstraps>

set -e
DIR="${1:?Uso: bash check_bootstraps.sh <directorio>}"
shopt -s nullglob
files=("$DIR"/bootstrapNetwork_*.txt)
[ ${#files[@]} -eq 0 ] && { echo "No hay archivos bootstrapNetwork_*.txt en $DIR"; exit 1; }

TMP=$(mktemp)
echo -e "archivo\tlineas\tregs_distintos\tultimo_regulador\tnewline" > "$TMP"
for f in "${files[@]}"; do
  read lines distinct lastreg <<< "$(awk -F'\t' 'NR>1{n++; reg[$1]; last=$1} END{print n, length(reg), last}' "$f")"
  [ -z "$(tail -c1 "$f")" ] && nl="SI" || nl="NO"
  echo -e "$(basename "$f")\t$lines\t$distinct\t$lastreg\t$nl" >> "$TMP"
done

# Referencia = moda (valor mas frecuente) de regs distintos y de ultimo regulador
ref_regs=$(tail -n +2 "$TMP" | cut -f3 | sort | uniq -c | sort -rn | head -1 | awk '{print $2}')
ref_last=$(tail -n +2 "$TMP" | cut -f4 | sort | uniq -c | sort -rn | head -1 | awk '{print $2}')

echo "== Referencia (bootstrap completo): regs_distintos=$ref_regs, ultimo_regulador=$ref_last =="
echo
echo "== Archivos INCOMPLETOS / sospechosos (si aparece alguno) =="
bad=0
while IFS=$'\t' read -r name lines distinct lastreg nl; do
  [ "$name" = "archivo" ] && continue
  motivo=""
  [ "$distinct" != "$ref_regs" ] && motivo="$motivo regs=$distinct(esperado $ref_regs)"
  [ "$lastreg" != "$ref_last" ] && motivo="$motivo ultimo=$lastreg(esperado $ref_last)"
  [ "$nl" != "SI" ] && motivo="$motivo sin_newline"
  if [ -n "$motivo" ]; then echo "  INCOMPLETO: $name ->$motivo"; bad=$((bad+1)); fi
done < "$TMP"

total=$(( ${#files[@]} ))
echo
echo "== RESUMEN: $((total-bad))/$total COMPLETOS, $bad incompletos =="
mv "$TMP" "$DIR/AUDIT_completitud.tsv"
echo "Tabla por archivo guardada en: $DIR/AUDIT_completitud.tsv"

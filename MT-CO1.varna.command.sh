#!/bin/bash
set -euo pipefail

VARNA_JAR="/g/data/lf10/tools/varna/VARNAv3-93.jar"
OUT_PNG="/scratch/lf10/rh1772/MAT2A/analysis/RNAseq/DMS_reactivity/MT-CO1.varna_reactivity.png"
OUT_SVG="/scratch/lf10/rh1772/MAT2A/analysis/RNAseq/DMS_reactivity/MT-CO1.varna_reactivity.svg"

SEQ="$(cat /scratch/lf10/rh1772/MAT2A/analysis/RNAseq/DMS_reactivity/MT-CO1.varna.sequence.txt)"
STRUCT="$(cat /scratch/lf10/rh1772/MAT2A/analysis/RNAseq/DMS_reactivity/MT-CO1.varna.structure.txt)"
VALUES="$(cat /scratch/lf10/rh1772/MAT2A/analysis/RNAseq/DMS_reactivity/MT-CO1.varna.values.txt)"

# PNG
java -cp "$VARNA_JAR" fr.orsay.lri.varna.applications.VARNAcmd \
  -sequenceDBN "$SEQ" \
  -structureDBN "$STRUCT" \
  -o "$OUT_PNG" \
  -algorithm naview \
  -resolution 8.0 \
  -colorMap "$VALUES" \
  -colorMapStyle "0.0:#FFFFFF;0.002:#FFF7BC;0.005:#FEC44F;0.01:#FE9929;0.02:#EC7014;0.05:#CC4C02" \
  -bp "#888888" \
  -title "MT-CO1 DMS reactivity (0.5% - 0%)"

# SVG
java -cp "$VARNA_JAR" fr.orsay.lri.varna.applications.VARNAcmd \
  -sequenceDBN "$SEQ" \
  -structureDBN "$STRUCT" \
  -o "$OUT_SVG" \
  -algorithm naview \
  -resolution 8.0 \
  -colorMap "$VALUES" \
  -colorMapStyle "0.0:#FFFFFF;0.002:#FFF7BC;0.005:#FEC44F;0.01:#FE9929;0.02:#EC7014;0.05:#CC4C02" \
  -bp "#888888" \
  -title "MT-CO1 DMS reactivity (0.5% - 0%)"

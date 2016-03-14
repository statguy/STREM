#!/bin/bash
# Run in figures output directory

function pdf2eps {
  if [ ! -z "$2" ]; then
    pdfcrop --bbox "$2" $1.pdf "$1-crop.pdf"
  else
    pdfcrop $1.pdf "$1-crop.pdf"
  fi

  convert -density 600 -units PixelsPerInch "$1-crop.pdf" -background white \
    -flatten -alpha off -auto-level -depth 8 -compress lzw $1.tif
  chmod a+r $1.tif
}

pdf2eps manuscript-figure-1
pdf2eps manuscript-figure-2
pdf2eps manuscript-figure-3
pdf2eps manuscript-figure-3sample
pdf2eps manuscript-figure-4a
pdf2eps manuscript-figure-5

# 1 2 481 219
pdf2eps manuscript-figure-3a1 "1 2 481 205"
pdf2eps manuscript-figure-3a2 "1 2 481 205"
pdf2eps manuscript-figure-3a3 "1 2 481 205"
pdf2eps manuscript-figure-3a4 "1 2 481 205"

mv manuscript-figure-1.tif Fig1.tif
mv manuscript-figure-2.tif Fig2.tif
mv manuscript-figure-3.tif Fig4.tif
mv manuscript-figure-3sample.tif Fig3.tif
mv manuscript-figure-4a.tif Fig6.tif
mv manuscript-figure-5.tif Fig5.tif

mv manuscript-figure-3a1.tif S1_Fig.tif
mv manuscript-figure-3a2.tif S2_Fig.tif
mv manuscript-figure-3a3.tif S3_Fig.tif
mv manuscript-figure-3a4.tif S4_Fig.tif

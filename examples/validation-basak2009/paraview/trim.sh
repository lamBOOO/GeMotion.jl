#!/bin/bash

# brew install imagemagick
# images come from standard MATLAB jpg export
# PV export with 1000x1000 linewidth 2pt and roughly maximum window size filled vertically on imac

format=png
destinationFolder=.

for file in *.$format; do
    filename=${file%.$format}
    convert -trim $filename.$format ${filename}.$format
    echo $filename converted.
done

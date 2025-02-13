#!/bin/bash

for p in 2 3 4 5 6
do
  for name in co-annulus_unstructured co-annulus_structured
  do
    gmsh -2 -o "$name"_"$p".msh -setnumber p "$p" "$name".geo
  done
done

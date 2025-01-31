#!/bin/bash

for p in 1 2 3 4 5 6
do
  for name in co-annulus_structured co-annulus_unstructured
  do
    gmsh -2 -o "$name"_"$p".msh -setnumber p "$p" "$name".geo
  done
done

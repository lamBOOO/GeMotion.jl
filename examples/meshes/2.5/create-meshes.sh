#!/bin/bash

for p in 1 2 3 4 5 6 7
do
  for name in co-annulus_structured co-annulus_unstructured
  do
    gmsh -2 -o "$name"_"$p".msh -setnumber p "$p" "$name".geo
  done
done

for p in 1 2 3 4 5 6
do
  for name in co-annulus_unstructured_anisotrop
  do
    gmsh -2 -o "$name"_"$p".msh -setnumber p "$p" "$name".geo
  done
done

for phi in 0 45 95
do
  gmsh -2 \
    -o co-annulus_unstructured_eccentric_e0p8_phi"$phi"_4.msh \
    -setnumber p 4 \
    -setnumber eccentricity 0.8 \
    -setnumber phi_deg "$phi" \
    co-annulus_unstructured_eccentric.geo
done

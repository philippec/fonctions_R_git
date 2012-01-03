#!/bin/sh

# Sources
FONCTIONSR_SRC="/Users/philippe/Site ftp/Site ftp.umontreal.ca/labo/fonctions_r/"
FONCTIONSR_DST="/Volumes/bio$/Casgrain/prog/labo/fonctions_r/"
LEGENDRE_SRC="/Users/philippe/Sites/legendre/"
LEGENDRE_DST="/Volumes/bio$/legendre/"
PDF_SRC="/Users/philippe/Site ftp/pdf legendre/"
PDF_DST="/Volumes/bio$/legendre/reprints/"


# Sync
echo "Fonctions R..."
rsync -avz "$FONCTIONSR_SRC" "$FONCTIONSR_DST" --exclude ".DS_Store"
echo "Site Web Legendre..."
rsync -avz "$LEGENDRE_SRC" "$LEGENDRE_DST" --exclude ".DS_Store"
echo "PDF Legendre..."
rsync -avz "$PDF_SRC" "$PDF_DST" --exclude ".DS_Store"

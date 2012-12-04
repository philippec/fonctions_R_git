#!/bin/sh

rsync -avz "/Users/philippe/Site ftp/Site ftp.umontreal.ca/labo/" "numericalecology@adn.biol.umontreal.ca:~/labo/"  --exclude ".DS_Store" --exclude .git

exit 0

# Old script

# Sources
FONCTIONSR_SRC="/Users/philippe/Site ftp/Site ftp.umontreal.ca/labo/fonctions_r/"
FONCTIONSR_DST="/Volumes/bio$/Casgrain/prog/labo/fonctions_r/"
LEGENDRE_SRC="/Users/philippe/Sites/legendre/"
LEGENDRE_DST="/Volumes/bio$/legendre/"
PDF_SRC="/Users/philippe/Site ftp/Reprints/"
PDF_DST="/Volumes/bio$/legendre/reprints/"

# Setup

osascript -e "try" -e 'mount volume "smb://casgrain:F00Bar78@10.112.40.128/bio$"' -e "end try"
if [ ! -d "/Volumes/bio$" ];
then
    echo "Could not mount bio$"
    exit 1
fi

# Sync
echo "Fonctions R..."
rsync -avz "$FONCTIONSR_SRC" "$FONCTIONSR_DST" --exclude ".DS_Store" --exclude .git
echo "Site Web Legendre..."
rsync -avz "$LEGENDRE_SRC" "$LEGENDRE_DST" --exclude ".DS_Store" --exclude .git
echo "PDF Legendre..."
rsync -avz "$PDF_SRC" "$PDF_DST" --exclude ".DS_Store" --exclude .git

# Teardown
osascript -e 'tell application "Finder" to eject disk "bio$"'

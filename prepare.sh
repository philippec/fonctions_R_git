#!/bin/sh

DEST="/Volumes/Louve/Users/philippe/Site ftp/Site ftp.umontreal.ca/labo/fonctions_r"

find . -name ".DS_Store" -delete

zip -f -D "$DEST/AdjustedRsquare.R.zip" "AdjustedRsquare.R"
zip -f -D "$DEST/anova.1way.R.zip" "anova.1way.R"
zip -f -D "$DEST/anova.2way.R.zip" "anova.2way.R"
zip -f -D "$DEST/anova.3way.R.zip" "anova.3way.R"
zip -f -D "$DEST/broken.stick.R.zip" "broken.stick.R"
zip -f -D "$DEST/corPerm.R.zip" "corPerm.R"
zip -f -D "$DEST/multRegress.R.zip" "multRegress.R"
zip -f -D "$DEST/nest.anova.perm.R.zip" "nest.anova.perm.R"
zip -f -D "$DEST/Sidak.R.zip" "Sidak.R"
zip -f -D "$DEST/t.perm.R.zip" "t.perm.R"
zip -f -D "$DEST/t.paired.perm.R.zip" "t.paired.perm.R"

zip -f -D -r "$DEST/CCorA.zip" "CCorA"

cd "PCA-CA/CA"
zip -f "$DEST/CA.zip" *
cd "../PCA"
zip -f "$DEST/PCA.zip" *
cd "../../"
#zip -f -j "$DEST/PCA.zip" "PCA-CA/PCA/*"

R CMD build kendall.W
R CMD build --binary kendall.W
R CMD build CADM
R CMD build --binary CADM
R CMD build lmorigin
R CMD build --binary lmorigin
R CMD build mantel.correlog
R CMD build --binary mantel.correlog
R CMD build parafit
R CMD build --binary parafit
R CMD build rdaTest
R CMD build --binary rdaTest

#cp "STI_1.0_win.zip" "$DEST/"
#cp "STI_1.0_mac.tar.gz" "$DEST/"
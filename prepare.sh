#!/bin/sh

DEST="/Users/philippe/Site ftp/Site ftp.umontreal.ca/labo/fonctions_r"

find . -name ".DS_Store" -delete

zip -f -D "$DEST/AdjustedRsquare.R.zip" "AdjustedRsquare.R"
zip -f -D "$DEST/anova.1way.R.zip" "anova.1way.R"
zip -f -D "$DEST/anova.2way.R.zip" "anova.2way.R"
zip -f -D "$DEST/anova.3way.R.zip" "anova.3way.R"
zip -f -D "$DEST/broken.stick.R.zip" "broken.stick.R"
zip -f -D "$DEST/corPerm.R.zip" "corPerm.R"
zip -f -D "$DEST/multRegress.R.zip" "multRegress.R"
zip -f -D "$DEST/nest.anova.perm.R.zip" "nest.anova.perm.R"
zip -f -D "$DEST/pcnm.all.R.zip" "pcnm.all.R"
zip -f -D "$DEST/Sidak.R.zip" "Sidak.R"

zip -f -D -r "$DEST/CCorA.zip" "CCorA"
zip -f -D -r "$DEST/rdaTest.zip" "rdaTest"

cp "quickPCNM.R" "$DEST/quickPCNM-770.R"
zip -f -D "$DEST/quickPCNM-770.R.zip" "$DEST/quickPCNM-770.R"
#zip -j "$DEST/quickPCNM-770.R.zip" "$DEST/quickPCNM-770.R"
rm "$DEST/quickPCNM-770.R"

cd "PCA-CA/CA"
zip -f "$DEST/CA.zip" *
cd "../PCA"
zip -f "$DEST/PCA.zip" *
cd "../../"
#zip -f -j "$DEST/PCA.zip" "PCA-CA/PCA/*"


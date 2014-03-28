#!/bin/sh

DEST="/Users/philippe/Projects/Legendre/prog/fonctions_r"

find . -name ".DS_Store" -delete
find "$DEST" -name ".DS_Store" -delete

zip -f -D "$DEST/AdjustedRsquare.R.zip" "AdjustedRsquare.R"
zip -f -D "$DEST/anova.1way.R.zip" "anova.1way.R"
zip -f -D "$DEST/anova.2way.R.zip" "anova.2way.R"
zip -f -D "$DEST/anova.2way.unbalanced.R.zip" "anova.2way.unbalanced.R"
zip -f -D "$DEST/anova.3way.R.zip" "anova.3way.R"
zip -f -D "$DEST/broken.stick.R.zip" "broken.stick.R"
zip -f -D "$DEST/corPerm.R.zip" "corPerm.R"
zip -f -D "$DEST/dagnelie.test.R.zip" "dagnelie.test.R"
zip -f -D "$DEST/dbRDA.D.R.zip" "dbRDA.D.R"
zip -f -D "$DEST/hclust.PL.R.zip" "hclust.PL.R"
zip -f -D "$DEST/multRegress.R.zip" "multRegress.R"
zip -f -D "$DEST/manova.2way.unbalanced.R.zip" "manova.2way.unbalanced.R"
zip -f -D "$DEST/nest.anova.perm.R.zip" "nest.anova.perm.R"
zip -f -D "$DEST/Sidak.R.zip" "Sidak.R"
zip -f -D "$DEST/t.perm.R.zip" "t.perm.R"
zip -f -D "$DEST/t.paired.perm.R.zip" "t.paired.perm.R"

zip -f -D -r "$DEST/CCorA.zip" "CCorA"
zip -f -D -r "$DEST/CCA.zip" "CCA"
zip -f -D -r "$DEST/seriation.zip" "seriation"
zip -f -D -r "$DEST/periodograph.zip" "periodograph"
zip -f -D -r "$DEST/Periodogram_W-R.zip" "Periodogram_W-R"
zip -f -D -r "$DEST/beta-diversity.zip" "beta-diversity"

cd "PCA-CA/CA"
zip -f "$DEST/CA.zip" *
cd "../PCA"
zip -f "$DEST/PCA.zip" *
cd "../../"
#zip -f -j "$DEST/PCA.zip" "PCA-CA/PCA/*"

R CMD build const.clust
cp const.clust_1.3.tar.gz "$DEST/"
R CMD INSTALL --build const.clust
cp const.clust_1.3.tgz "$DEST/"
R CMD build rdaTest
cp rdaTest_1.10.tar.gz "$DEST/"
R CMD INSTALL --build rdaTest
cp "rdaTest_1.10.tgz" "$DEST/"
cp "rdaTest_1.10.zip" "$DEST/"

#cp "STI_1.0_win.zip" "$DEST/"
#cp "STI_1.0_mac.tar.gz" "$DEST/"
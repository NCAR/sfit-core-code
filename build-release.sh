# Script to build a sfit4 official release
# This script is only intended to be used for creating repositories.
# Please do not use if you do not exactly know what all this is about

# First review the .gitattributes file and include all files
# or regex which are part of the repository but should not go into the repository

basedirname='SFIT4-Official-Release-1-0'

cd docs
find . -name '*.tex' -exec pdflatex {} ';'
cd ..

git archive -o $basedirname.tar --prefix $basedirname\/ V1.0.8 

# The documentation is created from the tex files


# Append the created pdf files to the release archive

find docs -name '*.pdf' -exec tar --transform 's,^\.,$basedirname,' f $basedirname.tar {} ';'

# compress tar file
gzip  -f SFIT4-Official-Release-1-0.tar



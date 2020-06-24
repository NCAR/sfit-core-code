# Script to build a sfit4 official release
# This script is only intended to be used for creating repositories.
# Please do not use if you do not exactly know what all this is about

# First review the .gitattributes file and include all files
# or regex which are part of the repository but should not go into the repository

basedirname='SFIT4-Official-Release-1-0'

git clean -f

cd docs
find . -name '*.tex' -exec pdflatex {} ';'
cd ..

git archive --worktree-attributes -o $basedirname.tar Official_Release_1.0 #V1.0.8 

# The documentation is created from the tex files


# Append the created pdf files to the release archive

find docs -name '*.pdf' -exec tar -rf $basedirname.tar {} ';'

# compress tar file
#gzip  -f SFIT4-Official-Release-1-0.tar



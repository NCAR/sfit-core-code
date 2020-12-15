# Script to build a sfit4 official release
# This script is only intended to be used for creating repositories.
# Please do not use if you do not exactly know what all this is about

# First review the .gitattributes file and include all files
# or regex which are part of the repository but should not go into the repository

function traverse() {   
    # walk the directory tree and execute comand in each directory
    cd $1
    for file in $(ls .); do
	if [ -d ${file} ] ; then
            traverse $file
	fi
    done
    echo $(pwd)
    find . -maxdepth 1 -name '*.tex' -exec pdflatex {} ';'
    find . -maxdepth 1 -name '*.tex' -exec pdflatex {} ';'
    find . -maxdepth 1 -name '*.doc*' -exec lowriter --convert-to pdf {} ';'
    find . -maxdepth 1 -name '*.txt*' -exec lowriter --convert-to pdf {} ';'
    find . -maxdepth 1 -name '*.xls*' -exec localc --convert-to pdf {} ';'
    find . -maxdepth 1 -name '*.ppt*' -exec loimpress --convert-to pdf {} ';'
    
    cd ..
}



basedirname='SFIT4-Official-Release-1-0'

# remove all files which are not in the repositiry
git clean -f

traverse docs

# Edit the file .fitattribues to exclude files and directories from the Release
git archive --prefix=$basedirname/ --worktree-attributes -o $basedirname.tar Official_Release_1.0 #V1.0.8 

# The documentation is created from the tex files


# Append the created pdf files to the release archive

find docs -name '*.pdf' -exec tar --transform 's,^,'$basedirname'/,' -rf $basedirname.tar {} ';'

# compress tar file
#gzip  -f SFIT4-Official-Release-1-0.tar

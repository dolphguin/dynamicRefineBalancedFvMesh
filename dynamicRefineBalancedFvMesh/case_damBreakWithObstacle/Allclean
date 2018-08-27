#!/bin/bash

WM_PROJECT_DIR=/opt/openfoam6
source $WM_PROJECT_DIR/etc/bashrc
echo "version: " && foamVersion

cd ${0%/*} || exit 1    # Run from this directory

# Source tutorial clean functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

rm system/cellSetDict > /dev/null 2>&1
rm -rf 0 > /dev/null 2>&1

cleanCase

#------------------------------------------------------------------------------

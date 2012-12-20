#!/bin/bash

PROJECT_DIR="$( cd "$( dirname "$0" )" && pwd )"
export PROTEOMICS_WINEPREFIX="/var/wine/$USER-pwiz"
$PROJECT_DIR/setup_proteomics_wine_env.sh
export WINEPREFIX=$PROTEOMICS_WINEPREFIX
export WINEARCH=win32
wine 'C:\Program Files\ProteoWizard\ProteoWizard 3.0.4019\msconvert.exe' $*

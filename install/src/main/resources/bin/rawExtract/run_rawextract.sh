#!/bin/bash

PROJECT_DIR="$( cd "$( dirname "$0" )" && pwd )"
export PROTEOMICS_WINEPREFIX="/var/wine/$USER-pwiz"
PWIZ_ENVIRONMENT=/mnt/mprc/software/public/bumbershoot/pwiz-wine-3_0_4019
$PWIZ_ENVIRONMENT/setup_proteomics_wine_env.sh
export WINEPREFIX=$PROTEOMICS_WINEPREFIX
export WINEARCH=win32
wine $PROJECT_DIR/MprcExtractRaw.exe "$@"

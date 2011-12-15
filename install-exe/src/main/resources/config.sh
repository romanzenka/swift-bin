#!/bin/bash

echo Starting Swift configuration

# Absolute path to this script
SCRIPT=`readlink -f $0`
# Absolute path to directory this script is in
SCRIPTPATH=`dirname $SCRIPT`

java  -Dswift.home="${SCRIPTPATH}" -Dswift.daemon=config -jar bin/swift/launcher.jar --config --war ./bin/swift/swift.war  $*

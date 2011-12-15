#!/bin/bash
source ~/.bashrc

echo Host name: `hostname`

umask 002


if [ $# -lt 1 ]
then
	echo "Usage: unixWrapper.sh <command>"
	echo "Where <command> is the command to be executed"
else
       echo "Executing command: $@"
       "$@"
       echo Done
fi

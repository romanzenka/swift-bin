#!/bin/bash
source ~/.bashrc

echo Host name: `hostname`

umask 002

if [ $# -lt 1 ] 
then
	echo "Usage: wrapper.sh <command>"
	echo "Where <command> is the command to be executed by the SGE"
else
	echo "Executing command: $@"
	exec "$@"
fi

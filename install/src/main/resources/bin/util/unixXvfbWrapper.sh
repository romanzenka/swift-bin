#!/bin/bash
source ~/.bashrc

echo Host name: `hostname`

umask 002

create_xvfb () {
	USERNAME=`whoami`
	DISPLAYNO=1
	while [ -z $xvfb_success ]
		do
		echo Xvfb :${DISPLAYNO} -screen 0 1024x1024x8
		Xvfb :${DISPLAYNO} -screen 0 1024x1024x8 >& /dev/null &
		XVFB_PID=$!
		echo $XVFB_PID
		sleep 1
		if ps --pid $XVFB_PID
               		then
                	echo "Started XVFB on display $DISPLAYNO process $XVFB_PID"
                	xvfb_success=1
                else
		        echo "Failed to use display "${DISPLAYNO}
                        DISPLAYNO=$(($DISPLAYNO + 1))
                        XVFB_PID=""
		fi
 		done
	export XVFB_PID
	export DISPLAY=:${DISPLAYNO}
	}

kill_xvfb () {
	kill $XVFB_PID
	}


if [ $# -lt 1 ]
then
	echo "Usage: unixXvfbWrapper <command>"
	echo "Where <command> is the command to be executed with Xvfb set up"
else
    echo Using Xvfb: `which Xvfb`

    if [ $? != 0 ]
    then
        echo "Xvfb command was not found. Install Xvfb and try again."
        exit 1
    fi

       create_xvfb

       echo "Executing command: $@"
       "$@"

       echo Killing Xvfb at $XVFB_PID
       kill_xvfb
       echo Done
fi

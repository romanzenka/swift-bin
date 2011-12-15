#!/bin/bash
echo Running Swift daemon
umask 002
java -Dswift.home=`pwd` -Dlog4j.configuration=file:conf/log4j.properties -jar bin/swift/swift.jar $*

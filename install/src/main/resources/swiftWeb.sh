#!/bin/bash

echo Starting Swift web server
umask 002
java -Dswift.home=`pwd` -Dlog4j.configuration=file:conf/log4j.properties -jar bin/swift/launcher.jar --war ./bin/swift/swift.war $*

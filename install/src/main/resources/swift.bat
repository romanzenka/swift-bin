@echo off
echo Running Swift daemon
java -Xdebug -Xrunjdwp:transport=dt_socket,server=y,suspend=n,address=5122 -jar bin\swift\swift.jar %*

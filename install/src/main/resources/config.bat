@echo off
echo Starting Swift configuration
java  -Dlog4j.configuration=file:conf/log4j.properties  -Xdebug -Xrunjdwp:transport=dt_socket,server=y,suspend=y,address=5200   -jar bin\swift\launcher.jar --config --war .\bin\swift\swift.war %*

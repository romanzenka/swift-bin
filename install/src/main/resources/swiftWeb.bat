@echo off
echo Starting Swift web server
java -Xdebug -Xrunjdwp:transport=dt_socket,server=y,suspend=n,address=5121 -jar bin\swift\launcher.jar --war .\bin\swift\swift.war %*

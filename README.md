Swift binaries
==============

This repository contains supplementary files that are used for Swift install creation and for integration testing.

We use a separate repository so the developers do not have to download a huge
amount of data if they want to do Swift development.

The binaries are deployed to a Maven repository at http://informatics.mayo.edu/maven and automatically downloaded during the build process.

### Executable files

A problem with the way maven creates assemblies forced us to split all the executable installation files into a separate install-exe package.

This should be eventually resolved in a better way - this is hard to manage.

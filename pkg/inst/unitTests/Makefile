TOP=../..
PKG=pkg
COMPLETE_NAME=nonlinearTseries_0.2.tar.gz
SUITE=doRUnit.R
R=R
all: install test

install: # Install package
	cd ${TOP}/..;\
	${R} CMD build ${PKG}
	cd ${TOP}/..;\
	${R} CMD INSTALL ${COMPLETE_NAME}
test: # Run unit tests
	export RCMDCHECK=FALSE;\
	cd ${TOP}/tests;\
	${R} --vanilla --slave < ${SUITE}

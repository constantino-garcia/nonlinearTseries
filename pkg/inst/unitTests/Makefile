TOP=../..
PKG=nonlinearAnalysis
COMPLETE_NAME=nonlinearAnalysis_0.1.tar.gz
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

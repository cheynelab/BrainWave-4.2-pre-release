# new makefile for compiling BW mex functions locally
# requires current copies of ctflib.o and bwlib.o in megcode directory
# VERSION 4.2 April 24, 23, moved libraries to same directory as brainwave folder

CTF_LIB = ../../bw_libraries/ctflib
BW_LIB = ../../bw_libraries/bwlib
# add different path to point to non-threaded bwlib library
BW_LIB_NO_THREADS = ../../bw_libraries/bwlib_nothreads

# updated Nov, 2016 for compilation on VM

# Matlab compiler version on Mac 10.14
# using XCode 10 Apple LLVM version 10.0.0 (clang-1000.11.45.2)
#
mex_imac64 = /Applications/MATLAB_R2021b.app/bin/mex

# added compile statements for Mac ARM chip - requires recompiled libraries! 
mex_amac64 = /Applications/MATLAB_R2023b.app/bin/mex

# Matlab compiler for CentOS 6.0 on VMware 
mex_linux = /usr/local/MATLAB/R2014b/bin/mex

# compile for Windows 10 on msys2 mingw64 from Cecilia's notes
# added static linking of stdlibcc and libstdc++ to eliminate need to have local copy of stdlibg++-6.dll
# inclusion of lpthread will require winpthread-1.dll

MATLABROOT=c:/Progra~1/MATLAB/R2018b
LIBS= -L$(MATLABROOT)/bin/win64 -L$(MATLABROOT)/extern/lib/win64/microsoft -lmex -lmx -lmwlapack -lmwblas -leng -ladvapi32 -luser32 -lgdi32 -lkernel32 -lmingwex

mex_win64= g++ -static-libgcc -static-libstdc++
MEXFLAG= -m64 -shared -DMATLAB_MEX_FILE -I$(MATLABROOT)/extern/include -Wl,--export-all-symbols $(LIBS)

win64:
	$(mex_win64) $(MEXFLAG) bw_geom2res4.cc -o bw_geom2res4.mexw64 $(CTF_LIB)/ctflib_win64.o -lws2_32
	$(mex_win64) $(MEXFLAG) bw_CTFGetParams.cc -o bw_CTFGetParams.mexw64 $(CTF_LIB)/ctflib_win64.o -lws2_32
	$(mex_win64) $(MEXFLAG) -DMX_COMPAT_32 bw_CTFGetHeader.cc -o bw_CTFGetHeader.mexw64 $(CTF_LIB)/ctflib_win64.o -lws2_32
	$(mex_win64) $(MEXFLAG) bw_CTFGetSensors.cc -o bw_CTFGetSensors.mexw64 $(CTF_LIB)/ctflib_win64.o -lws2_32
	$(mex_win64) $(MEXFLAG) bw_CTFGetChannelData.cc -o bw_CTFGetChannelData.mexw64 $(CTF_LIB)/ctflib_win64.o -lws2_32
	$(mex_win64) $(MEXFLAG) bw_CTFGetChannelLabels.cc -o bw_CTFGetChannelLabels.mexw64 $(CTF_LIB)/ctflib_win64.o -lws2_32
	$(mex_win64) $(MEXFLAG) bw_CTFEpochDs.cc -o bw_CTFEpochDs.mexw64 $(CTF_LIB)/ctflib_win64.o -lws2_32
	$(mex_win64) $(MEXFLAG) bw_CTFChangeHeadPos.cc -o bw_CTFChangeHeadPos.mexw64 $(CTF_LIB)/ctflib_win64.o -lws2_32
	$(mex_win64) $(MEXFLAG) bw_getCTFData.cc -o bw_getCTFData.mexw64 $(CTF_LIB)/ctflib_win64.o -lws2_32
	$(mex_win64) $(MEXFLAG) bw_filter.cc -o bw_filter.mexw64 $(CTF_LIB)/ctflib_win64.o -lws2_32
	$(mex_win64) $(MEXFLAG) bw_combineDs.cc -o bw_combineDs.mexw64 $(CTF_LIB)/ctflib_win64.o -lws2_32
	$(mex_win64) $(MEXFLAG) bw_kit2res4.cc -o bw_kit2res4.mexw64 $(CTF_LIB)/ctflib_win64.o -lws2_32
	$(mex_win64) $(MEXFLAG) bw_concatenateDs.cc -o bw_concatenateDs.mexw64 $(CTF_LIB)/ctflib_win64.o -lws2_32
	$(mex_win64) $(MEXFLAG) bw_fitDipole.cc -o bw_fitDipole.mexw64 $(CTF_LIB)/ctflib_win64.o -lws2_32
	$(mex_win64) $(MEXFLAG) bw_CTFNewDs.cc -o bw_CTFNewDs.mexw64 $(CTF_LIB)/ctflib_win64.o -lws2_32

	$(mex_win64) $(MEXFLAG) bw_CTFGetAverage.cc -o bw_CTFGetAverage.mexw64 $(CTF_LIB)/ctflib_win64.o $(BW_LIB_NO_THREADS)/bwlib_win64.o -lws2_32
	$(mex_win64) $(MEXFLAG) bw_makeVS.cc -o bw_makeVS.mexw64 $(CTF_LIB)/ctflib_win64.o $(BW_LIB_NO_THREADS)/bwlib_win64.o -lws2_32
	$(mex_win64) $(MEXFLAG) bw_makeEventRelated.cc -o bw_makeEventRelated.mexw64 $(CTF_LIB)/ctflib_win64.o $(BW_LIB_NO_THREADS)/bwlib_win64.o -lws2_32
	$(mex_win64) $(MEXFLAG) bw_makeDifferential.cc -o bw_makeDifferential.mexw64 $(CTF_LIB)/ctflib_win64.o $(BW_LIB_NO_THREADS)/bwlib_win64.o -lws2_32

	$(mex_win64) $(MEXFLAG) bw_computeFaceNormals.cc -o bw_computeFaceNormals.mexw64 -lws2_32
	$(mex_win64) $(MEXFLAG) trilinear.cpp -o trilinear.mexw64 -lws2_32

# move to parent directory on windows - assumes running on linux-like msys2 environment
	mv *.mexw64 ../

linux64:
	$(mex_linux) bw_geom2res4.cc $(CTF_LIB)/ctflib_glx64.o
	$(mex_linux) bw_CTFGetParams.cc $(CTF_LIB)/ctflib_glx64.o
	$(mex_linux) bw_CTFGetHeader.cc $(CTF_LIB)/ctflib_glx64.o
	$(mex_linux) bw_CTFGetSensors.cc $(CTF_LIB)/ctflib_glx64.o
	$(mex_linux) bw_CTFGetChannelData.cc $(CTF_LIB)/ctflib_glx64.o
	$(mex_linux) bw_CTFGetChannelLabels.cc $(CTF_LIB)/ctflib_glx64.o
	$(mex_linux) bw_CTFEpochDs.cc $(CTF_LIB)/ctflib_glx64.o
	$(mex_linux) bw_CTFChangeHeadPos.cc $(CTF_LIB)/ctflib_glx64.o
	$(mex_linux) bw_getCTFData.cc $(CTF_LIB)/ctflib_glx64.o
	$(mex_linux) bw_filter.cc $(CTF_LIB)/ctflib_glx64.o
	$(mex_linux) bw_combineDs.cc $(CTF_LIB)/ctflib_glx64.o
	$(mex_linux) bw_kit2res4.cc $(CTF_LIB)/ctflib_glx64.o
	$(mex_linux) bw_concatenateDs.cc $(CTF_LIB)/ctflib_glx64.o
	$(mex_linux) bw_fitDipole.cc $(CTF_LIB)/ctflib_glx64.o
	$(mex_linux) bw_CTFNewDs.cc $(CTF_LIB)/ctflib_glx64.o

	$(mex_linux) bw_CTFGetAverage.cc $(CTF_LIB)/ctflib_glx64.o $(BW_LIB)/bwlib_glx64.o
	$(mex_linux) bw_makeVS.cc $(CTF_LIB)/ctflib_glx64.o $(BW_LIB)/bwlib_glx64.o
	$(mex_linux) bw_makeEventRelated.cc $(CTF_LIB)/ctflib_glx64.o $(BW_LIB)/bwlib_glx64.o
	$(mex_linux) bw_makeDifferential.cc $(CTF_LIB)/ctflib_glx64.o $(BW_LIB)/bwlib_glx64.o

	$(mex_linux) bw_computeFaceNormals.cc
	$(mex_linux) trilinear.cpp
	mv *.mexa64 ../

test:
	$(mex_amac64) CXXFLAGS='-target arm64-apple-macos11 -Wno-deprecated-declarations' bw_getCTFData.cc $(CTF_LIB)/ctflib_maca64.o
	mv bw_getCTFData.mexmaca64 ../
	
imac64:
	$(mex_imac64) CXXFLAGS='-target x86_64-apple-macos10.9 -Wno-deprecated-declarations' bw_geom2res4.cc $(CTF_LIB)/ctflib_maci64.o
	$(mex_imac64) CXXFLAGS='-target x86_64-apple-macos10.9 -Wno-deprecated-declarations' bw_CTFGetParams.cc $(CTF_LIB)/ctflib_maci64.o
	$(mex_imac64) CXXFLAGS='-target x86_64-apple-macos10.9 -Wno-deprecated-declarations' -DMX_COMPAT_32 bw_CTFGetHeader.cc $(CTF_LIB)/ctflib_maci64.o
	$(mex_imac64) CXXFLAGS='-target x86_64-apple-macos10.9 -Wno-deprecated-declarations' bw_CTFGetSensors.cc $(CTF_LIB)/ctflib_maci64.o
	$(mex_imac64) CXXFLAGS='-target x86_64-apple-macos10.9 -Wno-deprecated-declarations' bw_CTFGetChannelData.cc $(CTF_LIB)/ctflib_maci64.o
	$(mex_imac64) CXXFLAGS='-target x86_64-apple-macos10.9 -Wno-deprecated-declarations' bw_CTFGetChannelLabels.cc $(CTF_LIB)/ctflib_maci64.o
	$(mex_imac64) CXXFLAGS='-target x86_64-apple-macos10.9 -Wno-deprecated-declarations' bw_CTFEpochDs.cc $(CTF_LIB)/ctflib_maci64.o
	$(mex_imac64) CXXFLAGS='-target x86_64-apple-macos10.9 -Wno-deprecated-declarations' bw_CTFChangeHeadPos.cc $(CTF_LIB)/ctflib_maci64.o
	$(mex_imac64) CXXFLAGS='-target x86_64-apple-macos10.9 -Wno-deprecated-declarations' bw_getCTFData.cc $(CTF_LIB)/ctflib_maci64.o
	$(mex_imac64) CXXFLAGS='-target x86_64-apple-macos10.9 -Wno-deprecated-declarations' bw_filter.cc $(CTF_LIB)/ctflib_maci64.o
	$(mex_imac64) CXXFLAGS='-target x86_64-apple-macos10.9 -Wno-deprecated-declarations' bw_combineDs.cc $(CTF_LIB)/ctflib_maci64.o
	$(mex_imac64) CXXFLAGS='-target x86_64-apple-macos10.9 -Wno-deprecated-declarations' bw_kit2res4.cc $(CTF_LIB)/ctflib_maci64.o
	$(mex_imac64) CXXFLAGS='-target x86_64-apple-macos10.9 -Wno-deprecated-declarations' bw_concatenateDs.cc $(CTF_LIB)/ctflib_maci64.o
	$(mex_imac64) CXXFLAGS='-target x86_64-apple-macos10.9 -Wno-deprecated-declarations' bw_fitDipole.cc $(CTF_LIB)/ctflib_maci64.o
	$(mex_imac64) CXXFLAGS='-target x86_64-apple-macos10.9 -Wno-deprecated-declarations' bw_CTFNewDs.cc $(CTF_LIB)/ctflib_maci64.o
	
	$(mex_imac64) CXXFLAGS='-target x86_64-apple-macos10.9 -Wno-deprecated-declarations' bw_CTFGetAverage.cc $(CTF_LIB)/ctflib_maci64.o $(BW_LIB)/bwlib_maci64.o
	$(mex_imac64) CXXFLAGS='-target x86_64-apple-macos10.9 -Wno-deprecated-declarations' bw_makeVS.cc $(CTF_LIB)/ctflib_maci64.o $(BW_LIB)/bwlib_maci64.o
	$(mex_imac64) CXXFLAGS='-target x86_64-apple-macos10.9 -Wno-deprecated-declarations' bw_makeEventRelated.cc $(CTF_LIB)/ctflib_maci64.o $(BW_LIB)/bwlib_maci64.o
	$(mex_imac64) CXXFLAGS='-target x86_64-apple-macos10.9 -Wno-deprecated-declarations' bw_makeDifferential.cc $(CTF_LIB)/ctflib_maci64.o $(BW_LIB)/bwlib_maci64.o

	$(mex_imac64) CXXFLAGS='-target x86_64-apple-macos10.9 -Wno-deprecated-declarations' bw_computeFaceNormals.cc
	$(mex_imac64) CXXFLAGS='-target x86_64-apple-macos10.9 -Wno-deprecated-declarations' trilinear.cpp

	mv *.mexmaci64 ../

amac64:
	$(mex_amac64) CXXFLAGS='-target arm64-apple-macos11 -Wno-deprecated-declarations' bw_geom2res4.cc $(CTF_LIB)/ctflib_maca64.o
	$(mex_amac64) CXXFLAGS='-target arm64-apple-macos11 -Wno-deprecated-declarations' bw_CTFGetParams.cc $(CTF_LIB)/ctflib_maca64.o
	$(mex_amac64) CXXFLAGS='-target arm64-apple-macos11 -Wno-deprecated-declarations' -DMX_COMPAT_32 bw_CTFGetHeader.cc $(CTF_LIB)/ctflib_maca64.o
	$(mex_amac64) CXXFLAGS='-target arm64-apple-macos11 -Wno-deprecated-declarations' bw_CTFGetSensors.cc $(CTF_LIB)/ctflib_maca64.o
	$(mex_amac64) CXXFLAGS='-target arm64-apple-macos11 -Wno-deprecated-declarations' bw_CTFGetChannelData.cc $(CTF_LIB)/ctflib_maca64.o
	$(mex_amac64) CXXFLAGS='-target arm64-apple-macos11 -Wno-deprecated-declarations' bw_CTFGetChannelLabels.cc $(CTF_LIB)/ctflib_maca64.o
	$(mex_amac64) CXXFLAGS='-target arm64-apple-macos11 -Wno-deprecated-declarations' bw_CTFEpochDs.cc $(CTF_LIB)/ctflib_maca64.o
	$(mex_amac64) CXXFLAGS='-target arm64-apple-macos11 -Wno-deprecated-declarations' bw_CTFChangeHeadPos.cc $(CTF_LIB)/ctflib_maca64.o
	$(mex_amac64) CXXFLAGS='-target arm64-apple-macos11 -Wno-deprecated-declarations' bw_getCTFData.cc $(CTF_LIB)/ctflib_maca64.o
	$(mex_amac64) CXXFLAGS='-target arm64-apple-macos11 -Wno-deprecated-declarations' bw_filter.cc $(CTF_LIB)/ctflib_maca64.o
	$(mex_amac64) CXXFLAGS='-target arm64-apple-macos11 -Wno-deprecated-declarations' bw_combineDs.cc $(CTF_LIB)/ctflib_maca64.o
	$(mex_amac64) CXXFLAGS='-target arm64-apple-macos11 -Wno-deprecated-declarations' bw_kit2res4.cc $(CTF_LIB)/ctflib_maca64.o
	$(mex_amac64) CXXFLAGS='-target arm64-apple-macos11 -Wno-deprecated-declarations' bw_concatenateDs.cc $(CTF_LIB)/ctflib_maca64.o
	$(mex_amac64) CXXFLAGS='-target arm64-apple-macos11 -Wno-deprecated-declarations' bw_fitDipole.cc $(CTF_LIB)/ctflib_maca64.o
	$(mex_amac64) CXXFLAGS='-target arm64-apple-macos11 -Wno-deprecated-declarations' bw_CTFNewDs.cc $(CTF_LIB)/ctflib_maca64.o

	$(mex_amac64) CXXFLAGS='-target arm64-apple-macos11 -Wno-deprecated-declarations' bw_CTFGetAverage.cc $(CTF_LIB)/ctflib_maca64.o $(BW_LIB)/bwlib_maca64.o
	$(mex_amac64) CXXFLAGS='-target arm64-apple-macos11 -Wno-deprecated-declarations' bw_makeVS.cc $(CTF_LIB)/ctflib_maca64.o $(BW_LIB)/bwlib_maca64.o
	$(mex_amac64) CXXFLAGS='-target arm64-apple-macos11 -Wno-deprecated-declarations' bw_makeEventRelated.cc $(CTF_LIB)/ctflib_maca64.o $(BW_LIB)/bwlib_maca64.o
	$(mex_amac64) CXXFLAGS='-target arm64-apple-macos11 -Wno-deprecated-declarations' bw_makeDifferential.cc $(CTF_LIB)/ctflib_maca64.o $(BW_LIB)/bwlib_maca64.o

	$(mex_amac64) CXXFLAGS='-target arm64-apple-macos11 -Wno-deprecated-declarations' bw_computeFaceNormals.cc
	$(mex_amac64) CXXFLAGS='-target arm64-apple-macos11 -Wno-deprecated-declarations' trilinear.cpp

	mv *.mexmaca64 ../

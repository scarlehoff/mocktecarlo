#!/usr/bin/env scons

target  = "mocktecarlo"

## Cuba 
cubaSrc = "Cuba-4.2"
from subprocess import check_output
## Lhapdf 
lhapdflib  = check_output(["lhapdf-config", "--libdir"]).strip()
lhapdfinc  = check_output(["lhapdf-config", "--incdir"]).strip()

include = [lhapdfinc, cubaSrc + "/include"]
includeflags = []
for i in include: includeflags.append('-I' + i)

source  = ["main.cpp", "crossSection.cpp", "MomentumSet.cpp", 
           "phaseSpace.cpp", "FourVector.cpp", "matrixElement.cpp"]
libpath = [cubaSrc + "/lib", lhapdflib]
libs    = ['cuba', 'm', 'LHAPDF']
ccflags = includeflags + ["-std=c++11"]

env = Environment(CCFLAGS = ccflags, LIBS = libs, LIBPATH = libpath)
env.Program(target = target, source = source)

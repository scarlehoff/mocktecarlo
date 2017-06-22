#!/usr/bin/env scons

target  = "mocktecarlo"

## Cuba 
cubaSrc = "Cuba-4.2"
from subprocess import check_output
## Lhapdf 
lhapdflib  = check_output(["lhapdf-config", "--libdir"]).strip()
lhapdfinc  = check_output(["lhapdf-config", "--incdir"]).strip()
## Fastjet
fastjetsrc = check_output(["fastjet-config", "--prefix"]).strip()
fastjetlib = fastjetsrc + "/lib"
fastjetinc = fastjetsrc + "/include"
fastjetlf  = ["fastjettools", "fastjet", "fastjetplugins", "siscone", "siscone_spherical"]

## Generate include flags
include = [lhapdfinc, cubaSrc + "/include", fastjetinc]
includeflags = []
for i in include: 
    includeflags.append('-I' + i)

source_files  = ["main.cpp", "CrossSection.cpp", "MomentumSet.cpp", 
           "PhaseSpace.cpp", "FourVector.cpp", "MatrixElement.cpp", "SubtractionTerm.cpp"]
source=[]
for filename in source_files:
    source.append("src/" + filename)

libpath = [cubaSrc + "/lib", lhapdflib, fastjetlib]
libs    = ['cuba', 'm', 'LHAPDF', 'gsl', 'gslcblas'] + fastjetlf
# debug
dbflags = ["-g", "-O2", 
#            "-Q", # will show which function is causing it to crash
#            "-v", # shows how ccl was invoked
           ]
ccflags = includeflags + ["-std=c++11"] # + dbflags


env = Environment(CCFLAGS = ccflags, LIBS = libs, LIBPATH = libpath, CPPPATH="include")
env.Program(target = target, source = source)

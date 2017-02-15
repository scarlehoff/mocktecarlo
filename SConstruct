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
for i in include: includeflags.append('-I' + i)

source  = ["main.cpp", "crossSection.cpp", "MomentumSet.cpp", 
           "phaseSpace.cpp", "FourVector.cpp", "matrixElement.cpp", "subtractionTerm.cpp"]
libpath = [cubaSrc + "/lib", lhapdflib, fastjetlib]
libs    = ['cuba', 'm', 'LHAPDF'] + fastjetlf
ccflags = includeflags + ["-std=c++11"]


env = Environment(CCFLAGS = ccflags, LIBS = libs, LIBPATH = libpath)
env.Program(target = target, source = source)

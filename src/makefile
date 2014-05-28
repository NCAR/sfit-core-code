# To remove all the objects but leave the executables use the
# command line 	make -f makefile clean

# Note: object channel.o is added by WUJIAN PENG in July, 2002
#       object planck.o added by M. Palm in 2008
#       arguments for g95 anf gfortran added by M. Palm in 2008
.SUFFIXES:

.SUFFIXES: .f .f90 .o .mod .F

# Fixed source, if any:
#.f.o:
#	$(FC) -c $(FFLAGS) $<
#
#.F.o:
#	$(FC) -c $(FFLAGS) $<

#.f90.mod:
#	$(FC) -c $(FFLAGS) $<

# Free form source:
.f90.o:
	$(FC) -c $(FFLAGS) $<

#%.mod %.o: %,f90
#	$(FC) -c $(FFLAGS) $<

FC = gfortran

# error check everything
FFLAGS = -Warray-bounds -Wcharacter-truncation -Wunused-parameter \
         -Wimplicit-interface -Wunderflow -Wconversion-extra \
         -Wunused-dummy-argument -Wcompare-reals \
         -fcheck=all -fimplicit-none -fwhole-file -fall-intrinsics -fno-automatic \
         -fstack-arrays -frepack-arrays \
         -std=gnu -pedantic \
         -ffpe-trap=invalid,zero,overflow,underflow -g -fbacktrace -fsanitize=address

# -Warray-temporaries too many...

#error check a lot
#FFLAGS = -Wall -fno-automatic -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -Warray-bounds \
#         -fbounds-check -fbacktrace

# optomize
FFLAGS = -Wall -fno-automatic -frepack-arrays

# supply compile date:
# -cpp -DBDATE="\"`date +%Y%m%d`"\" -O2 -fno-automatic

# intel
#FC = ifort
#FFLAGS = -fpe-all=0  -O0 -debug full

#FFLAGS = -debug full


# Portland compiler
#FC = pgf90
#FFLAGS = ${PGBUG} ${PGTRAP}
#FFLAGS = -traceback -g -Kieee  ${PGTRAP} ${PGBUG}

PGTRAP = -Ktrap=divz -Ktrap=align -Ktrap=denorm  -Ktrap=inv -Ktrap=ovf -Ktrap=unf
PGBUG = -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -Mpgicoff
PGOPT = -fast

# use intermnal matrix inversion (not lapack)
# no eigenvector output
# mac:
#VLIB = -Wl,-framework -Wl,vecLib

# linux:
#VLIB = -L/usr/lib64 -llapack

LIBS =

VLIBS =

OBJS = \
bandparam.o \
binput_parse_4_0.o \
binput_4_0.o \
channel.o \
datafiles.o \
diagnostic.o \
frwdmdl.o \
galatry.o \
initialize.o \
isotope.o \
lineparam.o \
matrix.o \
molcparam.o \
opt.o \
params.o \
retvparam.o \
raytrace.o\
readin.o \
sfit4.o \
solar.o \
synspec.o \
tips.o \
transmis.o \
vibfcn.o \
voigt.o \
voigt_sdv_lm.o \
writeout.o \
andex.o \
xsections.o

#convert_old_data_files.o

all : sfit4 hbin pspec # convert_rdrv_394

sfit4: $(OBJS) rdrv40.o
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS) $(VLIB)

#convert_rdrv_394 : convert_rdrv_394.o $(OBJS) convert_old_data_files.o read_binput_394.o
#	$(FC) -o $@ $^ $(FFLAGS) $(LIBS) $(VLIB)

hbin: hbin.o $(OBJS)
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS) $(VLIB)

pspec: pspec.o spec.o fitting.o  $(OBJS)
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS) $(VLIB)

#params.o binput_4_0.o binput_parse_4_0.o datafiles.o bandparam.o

clean:
	\rm -f *.o *.mod *.lst
	\rm -rf *.dSYM

cleaner:
	\rm -f *.o *.mod *.lst sfit4 convert_rdrv_394 hbin pspec
	\rm -rf *.dSYM

# following is the list of what files use each included file
#    all names listed use the included file shown after the ':'

#    This is needed so make can check if included files have changed and
#    so .o files need to be recompiled.

#    EXAMPLE:
#     cocalc.o colist.o : bands.dat

#     in the example cocalc.o and colist.o use the included file bands.dat

bandparam.o \
: params.o

binput_4_0 \
: params.o binput_parse_4_0.o retvparam.o bandparam.o lineparam.o solar.o channel.o initialize.o opt.o writeout.o

binput_parse_4_0.o \
: params.o retvparam.o frwdmdl.o bandparam.o lineparam.o solar.o channel.o initialize.o opt.o datafiles.o isotope.o writeout.o

channel.o \
: params.o writeout.o synspec.o datafiles.o

convert_old_data_files.o \
: params.o bandparam.o retvparam.o

convert_rdrv_394.o \
: datafiles.o params.o binput_4_0.o convert_old_data_files.o

datafiles.o \
: params.o

diagnostic.o \
: opt.o retvparam.o

frwdmdl.o \
: params.o transmis.o datafiles.o molcparam.o solar.o synspec.o channel.o bandparam.o retvparam.o writeout.o initialize.o

galatry.o \
: params.o

hbin.o \
: params.o datafiles.o bandparam.o lineparam.o isotope.o binput_4_0.o binput_parse_4_0.o

initialize.o \
: params.o retvparam.o transmis.o molcparam.o xsections.o datafiles.o synspec.o lineparam.o bandparam.o raytrace.o solar.o writeout.o channel.o

isotope.o \
: params.o molcparam.o retvparam.o vibfcn.o datafiles.o

lineparam.o \
: params.o bandparam.o retvparam.o datafiles.o molcparam.o isotope.o raytrace.o

matrix.o \
: params.o

molcparam.o \
: params.o

opt.o \
: params.o datafiles.o frwdmdl.o synspec.o matrix.o writeout.o bandparam.o

pspec.o \
: spec.o params.o datafiles.o bandparam.o lineparam.o isotope.o binput_4_0.o binput_parse_4_0.o fitting.o

spec.o \
: fitting.o

raytrace.o \
: params.o molcparam.o datafiles.o retvparam.o isotope.o writeout.o andex.o

rdrv.o \
: params.o datafiles.o binput_4_0.o

read_binput_394.o \
: params.o binput_parse_4_0.o retvparam.o bandparam.o lineparam.o solar.o channel.o initialize.o opt.o writeout.o

readin.o \
: params.o retvparam.o transmis.o molcparam.o xsections.o datafiles.o synspec.o lineparam.o solar.o bandparam.o initialize.o opt.o channel.o

retvparam.o \
: params.o

sfit4.o \
: params.o vibfcn.o retvparam.o frwdmdl.o bandparam.o molcparam.o datafiles.o solar.o opt.o diagnostic.o synspec.o \
lineparam.o initialize.o readin.o writeout.o binput_4_0.o raytrace.o

solar.o \
: params.o datafiles.o xsections.o writeout.o

synspec.o \
: params.o transmis.o solar.o bandparam.o

tips.o \
: params.o

transmis.o \
: params.o bandparam.o retvparam.o vibfcn.o xsections.o molcparam.o lineparam.o

vibfcn.o \
: params.o retvparam.o molcparam.o #writeout.o

voigt_sdv_lm.o \
: params.o

voigt.o \
: params.o

writeout.o \
: params.o datafiles.o retvparam.o molcparam.o isotope.o bandparam.o

xsections.o \
: params.o retvparam.o vibfcn.o tips.o isotope.o lineparam.o bandparam.o molcparam.o voigt_sdv_lm.o

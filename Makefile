####for mpi/parallel compilation

#MPICH_DIR =/projects/cesm/devtools/mpich-3.0.4-gcc4.8.1
#MPICH_DIR = /opt/openmpi/1.6.4/gcc4
#CXX = ${MPICH_DIR}/bin/mpic++
#CFLAGS = -DPMODE -g3 -O0 -fPIC -fpermissive -fmessage-length=0 -I${MPICH_DIR}/include
#LDFLAGS = -L${MPICH_DIR}/lib
#RM = rm -f

######for non-mpi compilation
CXX = g++
CFLAGS = -O2

CLMOBJ = tclm602.o  

LULCOBJ = tbiome603c.o tlcluc603c.o

TEMOBJ = ttem604.o atms602.o tveg604.o tmcrb604.o \
tsoil604.o qsoiltemp603.o  humnact604.o

ELMNTOBJ = elmnt602.o telm604.o latdat602.o \
tclmdat602.o tco2dat437.o tndepdat602.o \
lulcdat44a.o tmaxcohortdat437.o \
ttemdat603a.o telvdat602.o tsoldat602.o 

OBJ =  tprocessXML602.o

aggrmap603 : tbiome603c.o tprocessXML602.o tclmdat602.o tmaxcohortdat437.o ttemdat603a.o aggrmap603.cxx
	${CXX} ${CFLAGS} -o aggrmap603 tbiome603c.o tprocessXML602.o tclmdat602.o tmaxcohortdat437.o ttemdat603a.o aggrmap603.cxx

mapaubChina10km603 : tbiome603c.o tprocessXML602.o tclmdat602.o tmaxcohortdat437.o ttemdat603a.o mapaubChina10km603.cxx
	${CXX} ${CFLAGS} -o mapaubChina10km603 tbiome603c.o tprocessXML602.o tclmdat602.o tmaxcohortdat437.o ttemdat603a.o mapaubChina10km603.cxx

mapmblChina10km603 : tbiome603c.o tprocessXML602.o tclmdat602.o tmaxcohortdat437.o ttemdat603a.o mapmblChina10km603.cxx
	${CXX} ${CFLAGS} -o mapmblChina10km603 tbiome603c.o tprocessXML602.o tclmdat602.o tmaxcohortdat437.o ttemdat603a.o mapmblChina10km603.cxx

xtem604 : ${CLMOBJ} ${LULCOBJ} ${TEMOBJ} ${ELMNTOBJ} ${OBJ} ptem604.cxx 
	${CXX} ${CFLAGS} ${LDFLAGS} -o xtem604 ${CLMOBJ} ${LULCOBJ} ${TEMOBJ} ${ELMNTOBJ} ${OBJ} ptem604.cxx 

clean:
	$(RM) ${CLMOBJ} ${LULCOBJ} ${TEMOBJ} ${ELMNTOBJ} ${OBJ} xtem604 

atms602.o : atms602.cpp atms602.h
	${CXX} ${CFLAGS} -c atms602.cpp

elmnt602.o : elmnt602.cpp elmnt602.h
	${CXX} ${CFLAGS} -c elmnt602.cpp

humnact604.o : humnact604.cpp humnact604.h
	${CXX} ${CFLAGS} -c humnact604.cpp

latdat602.o : latdat602.cpp latdat602.h
	${CXX} ${CFLAGS} -c latdat602.cpp

lulcdat44a.o : lulcdat44a.cpp lulcdat44a.h
	${CXX} ${CFLAGS} -c lulcdat44a.cpp

qsoiltemp603.o : qsoiltemp603.cpp qsoiltemp603.h
	${CXX} ${CFLAGS} -c qsoiltemp603.cpp

tbiome603c.o : tbiome603c.cpp tbiome603c.h
	${CXX} ${CFLAGS} -c tbiome603c.cpp

tclm602.o : tclm602.cpp tclm602.h
	${CXX} ${CFLAGS} -c tclm602.cpp

tclmdat602.o : tclmdat602.cpp tclmdat602.h
	${CXX} ${CFLAGS} -c tclmdat602.cpp

tco2dat437.o : tco2dat437.cpp tco2dat437.h
	${CXX} ${CFLAGS} -c tco2dat437.cpp

telvdat602.o : telvdat602.cpp telvdat602.h
	${CXX} ${CFLAGS} -c telvdat602.cpp

telm604.o : telm604.cpp telm604.h
	${CXX} ${CFLAGS} -c telm604.cpp

tlcluc603c.o : tlcluc603c.cpp tlcluc603c.h
	${CXX} ${CFLAGS} -c tlcluc603c.cpp

tmaxcohortdat437.o : tmaxcohortdat437.cpp tmaxcohortdat437.h
	${CXX} ${CFLAGS} -c tmaxcohortdat437.cpp

tmcrb604.o : tmcrb604.cpp tmcrb604.h
	${CXX} ${CFLAGS} -c tmcrb604.cpp

tndepdat602.o : tndepdat602.cpp tndepdat602.h
	${CXX} ${CFLAGS} -c tndepdat602.cpp

tprocessXML602.o : tprocessXML602.cpp tprocessXML602.h
	${CXX} ${CFLAGS} -c tprocessXML602.cpp

tsoil604.o : tsoil604.cpp tsoil604.h
	${CXX} ${CFLAGS} -c tsoil604.cpp

tsoldat602.o : tsoldat602.cpp tsoldat602.h
	${CXX} ${CFLAGS} -c tsoldat602.cpp

ttem604.o : ttem604.cpp ttem604.h
	${CXX} ${CFLAGS} -c ttem604.cpp

ttemdat603a.o : ttemdat603a.cpp ttemdat603a.h
	${CXX} ${CFLAGS} -c ttemdat603a.cpp

tveg604.o : tveg604.cpp tveg604.h
	${CXX} ${CFLAGS} -c tveg604.cpp


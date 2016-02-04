MEXSUFFIX  = mexa64
MATLABHOME = /share/Applic.ZIV/Matlab/R2013a
MEX        = pgc++
MEXCXX     = echo
CXX	       = pgc++

#GPUARCH = -acc -ta=nvidia,cc30,keepgpu
GPUARCH = 
GPUFLAGS = -fast $(GPUARCH) -Minfo=accel # pgcc
#GPUFLAGS = -std=c99 # gcc

#CFLAGS     = -fPIC -ansi -pthread -DMX_COMPAT_32 -O2 -DMATLAB_MEX_FILE
#CFLAGS = -fPIC -DMX_COMPAT_32 -O2 -DMATLAB_MEX_FILE
#CFLAGS = -fPIC -DMX_COMPAT_32 -DMATLAB_MEX_FILE
# -pthread -fno-omit-frame-pointer -ansi

# ccff = profile data
#CFLAGS = -Minfo=ccff -fPIC -lpthread -D_GNU_SOURCE  -O2 -DNDEBUG -B -csuffix=cpp $(GPUFLAGS)


CFLAGS = -fPIC -lpthread --c++11 -D_GNU_SOURCE  -O2 -DNDEBUG -B -csuffix=cpp $(GPUFLAGS)


#CFLAGS = -lpthread -DMATLAB_MEX_FILE -D_GNU_SOURCE  -O2 -DNDEBUG -acc -fast -ta=nvidia,cc30
#-Minfo 

LIBS      = 
INCLUDE   = -I$(MATLABHOME)/extern/include -Icommon
#MEXFLAGS =  -shared -Wl,-rpath-link,$(MATLABHOME)/bin/glnxa64 -L$(MATLABHOME)/bin/glnxa64 -lmx -lmex -lmat -lm -lpthread -Minfo
MEXFLAGS =  -shared -DMATLAB_MEX_FILE -Wl,--no-undefined -Wl,-rpath-link,$(MATLABHOME)/bin/glnxa64 -L$(MATLABHOME)/bin/glnxa64 -lmx -lmex -lmat -lm -lpthread $(CFLAGS)
PROJECTS = mf_iterate mf_iterate_standalone
MEXDIR = ..

all: $(PROJECTS)

# mf_iterate/functions.o common/tools.o
mf_iterate: mf_iterate/functions.o common/tools_c.o common/mextools.o common/mf_iterate.o
	$(MEX) $(MEXFLAGS) $(LIBS) -o $(MEXDIR)/$@.$(MEXSUFFIX) $^
	
mf_iterate_standalone: mf_iterate_standalone/functions.o common/tools_c.o common/mf_iterate.o
	$(CXX) $(CFLAGS) $(LIBS) -o $(MEXDIR)/$@ $^
	
.cpp.o:
	$(CXX) -c -o $@ $< $(CFLAGS) $(INCLUDE)
.c.o:
	$(CXX) -c -o $@ $< $(CFLAGS) $(INCLUDE)

clean:
	rm -f common/*.o
	for proj in $(PROJECTS); do \
		rm -f $$proj/*.o; \
		rm -f $(MEXDIR)/$$proj.$(MEXSUFFIX); \
		rm -f $(MEXDIR)/$$proj; \
	done

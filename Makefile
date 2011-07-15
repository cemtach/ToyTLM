

#
# C compiler specification
#
# CC = cc
CC = gcc
#
# C optimization, only one CFLAGS line should be uncommented
# The -Xcpluscomm flags is required for the SGI cc compiler to recognize the C++
# style comments ( // ) used in the code.
#
#CFLAGS = -Xcpluscomm 
#CFLAGS = -Xcpluscomm -O
#CFLAGS = -Xcpluscomm -O2
#CFLAGS = -Xcpluscomm -O3
#CFLAGS = -Xcpluscomm -n32
#CFLAGS = -Xcpluscomm -O -n32 
#CFLAGS = -Xcpluscomm -O2 -n32
#CFLAGS = -Xcpluscomm -O3 -n32
#CFLAGS = -Xcpluscomm -g
CFLAGS = -Wall -g
MAKE   = make
RM     = /bin/rm


all:
	@echo "To make an executable, type one of the following:"
	@echo " "
	@echo " $(MAKE) ToyTLMc"
	@echo " "
	@echo "To test one of the executables, type one of the following:"
	@echo " "
	@echo " $(MAKE) testc"
	@echo " "



.PHONY : clean cleandata all depend testc


#
# Creation of ToyFDTD from C source code
#
ToyTLMc: ToyTLMc.o Makefile
	${CC} ${CFLAGS} -o ToyTLMc ToyTLMc.o -lm


ToyTLMc.o: ToyTLM.c Makefile
	${CC} ${CFLAGS} -c -o ToyTLMc.o ToyTLM.c


# Make this more GNUish
check: testc

testc: ToyTLMc cleandata
	time ./ToyTLMc > c_runLog


#
# Type: make clean
# to remove executables, core files, object files, et cetera
#
clean: cleandata
	-$(RM) -f ToyTLMc *.o *~ core
#
# Type: make cleandata
# to remove only data files created by ToyFDTD
#
cleandata:
	-$(RM) -f c_*0.bob
	-$(RM) -f c_*1.bob
	-$(RM) -f c_*2.bob
	-$(RM) -f c_*3.bob
	-$(RM) -f c_*4.bob
	-$(RM) -f c_*5.bob
	-$(RM) -f c_*6.bob
	-$(RM) -f c_*7.bob
	-$(RM) -f c_*8.bob
	-$(RM) -f c_*9.bob
	-$(RM) -f c_runLog
	-$(RM) -f ToyTLMc.viz
	-$(RM) -f c_*.dat
 


.PHONY : clean cleandata all depend testc check






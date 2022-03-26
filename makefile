OBJ = RACode.o LDPCCode.o TurboCode.o AddNoise.o GF2Xlib.o perm.o indexx.o
NTL = -lntl -lgmp -L/usr/local/lib

example: example.o $(OBJ)
	g++ example.o $(OBJ) $(NTL)
fig1: fig1.o $(OBJ) RATurboCode.o
	g++ fig1.o $(OBJ) RATurboCode.o $(NTL)
table1: table1.o $(OBJ)
	g++ table1.o $(OBJ) $(NTL)
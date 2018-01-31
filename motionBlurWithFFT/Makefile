CC = g++
MAIN = motionblur
OBJS = fft.o motionblur.o
LIBS = -lnetpbm 
INCLUDES = -I /usr/include/netpbm
CFLAGS = -Wall
SRCS = fft.cc

all: $(MAIN)

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) -o $(MAIN) $(OBJS)

fft.o: $(SRCS)
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) -c $(SRCS)

%.o: %.cc
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) -c $< 

clean:
	rm -f *~ *.o 

clean-all: clean
	rm -f $(MAIN)
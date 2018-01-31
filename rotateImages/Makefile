all: image_rotate

image_rotate: image_rotate.cc
	g++ -I /usr/include/netpbm -o image_rotate image_rotate.cc -Wall -lnetpbm

clean:
	rm -f *~

clean-all: clean
	rm -f image_rotate
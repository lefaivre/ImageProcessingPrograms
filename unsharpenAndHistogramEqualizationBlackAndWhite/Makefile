all: unsharp hist_eq

unsharp: unsharp.cc
	g++ -I /usr/include/netpbm -o unsharp unsharp.cc -Wall -lnetpbm

hist_eq: hist_eq.cc
	g++ -I /usr/include/netpbm -o hist_eq hist_eq.cc -Wall -lnetpbm

clean:
	rm -f *~

clean-all: clean
	rm -f unsharp hist_eq

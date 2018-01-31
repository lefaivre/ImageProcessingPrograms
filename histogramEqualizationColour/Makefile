all: hist_eq

hist_eq: hist_eq.cc
	g++ -I /usr/include/netpbm -o hist_eq hist_eq.cc -Wall -lnetpbm

clean:
	rm -f *~

clean-all: clean
	rm -f hist_eq

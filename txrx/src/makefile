.PHONY: clean
all: linreg_mc.so
linreg_mc.so: linreg_mc.cpp
	$(CXX) -I./libeigen linreg_mc.cpp -o linreg_mc.so -shared -fPIC -O3 -march=native
clean: 
	-rm linreg_mc.so


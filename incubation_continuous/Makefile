CXX = g++
CXXFLAGS = -O2 -g -std=gnu++11
LDLIBS = -lm

SRCS = $(shell ls *.cpp)
OBJS = $(subst .cpp,.o,$(SRCS))
APPLS = incubation_continuous.exe

all: $(APPLS)

%.exe: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) $(APPLS)

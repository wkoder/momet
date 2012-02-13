CXX = g++
CXXFLAGS = -O3 -g -Wall -fmessage-length=0 -I/usr/include/python2.7 -fPIC
#CXXFLAGS = -g -Wall -fmessage-length=0
SRC = src/
BIN = bin/
VPATH = $(SRC):
OBJS = $(patsubst $(SRC)%.cpp, $(BIN)%.o, $(wildcard $(SRC)*.cpp)) \
		$(patsubst $(SRC)%.c, $(BIN)%.o, $(wildcard $(SRC)*.c))
DEPS = $(OBJS)
EXE = momet.so

init:
	mkdir -p $(BIN)

$(BIN)%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
$(BIN)%.o: %.c
	$(CXX) $(CXXFLAGS) -c $< -o $@
$(EXE): init $(DEPS)
	$(CXX) $(CXXFLAGS) -shared -Wl,--export-dynamic -lboost_python $(OBJS) -o $(BIN)$(EXE)
	
all:	clean $(EXE)

clean:
	rm -f $(BIN)*.o $(BIN)$(EXE)

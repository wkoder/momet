CXX = g++
PYTHON_VER = 2.7
CXXFLAGS = -O3 -g -Wall -fmessage-length=0 -I/usr/include/python$(PYTHON_VER) -fPIC
SRC = src/
BIN = bin/
VPATH = $(SRC):
OBJS = $(patsubst $(SRC)%.cpp, $(BIN)%.o, $(wildcard $(SRC)*.cpp)) \
		$(patsubst $(SRC)%.c, $(BIN)%.o, $(wildcard $(SRC)*.c))
DEPS = $(OBJS)
EXE = momet.so
INSTALL_DIR = /usr/lib/pymodules/python$(PYTHON_VER)

all:	clean $(EXE)

init:
	mkdir -p $(BIN)

clean:
	rm -f $(BIN)*.o $(BIN)$(EXE)
	
install:
	cp $(BIN)$(EXE) $(INSTALL_DIR)

$(BIN)%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
$(BIN)%.o: %.c
	$(CXX) $(CXXFLAGS) -c $< -o $@
$(EXE): init $(DEPS)
	$(CXX) $(CXXFLAGS) -shared -Wl,-soname,$(BIN)$(EXE) -o $(BIN)$(EXE) $(OBJS) -lpython$(PYTHON_VER) -lboost_python

######################################
#
######################################
#source file
#源文件，自动找所有.c, .cc 和.cpp文件，并将目标定义为同名.o文件
SOURCE  := $(wildcard *.c) $(wildcard *.cc) $(wildcard *.cpp)
OBJS    := $(patsubst %.c,%.o,$(patsubst %.cc,%.o,$(patsubst %.cpp,%.o,$(SOURCE))))

#target you can change test to what you want
#目标文件名，输入任意你想要的执行文件名
TARGET  := Main_Program

#compile and lib parameter
#编译参数
CC      := g++
LDFLAGS :=  -Wl,-rpath /home/shihao/snopt7/lib -L /home/shihao/snopt7/lib
LIBS    := -lsnopt_cpp -lsnopt7 -lsnblas -lsnprint7 -lgfortran -lf2c -lm -lpython2.7
INCLUDE := -I/home/shihao/snopt7/cppsrc/ -I/usr/include/python2.7 -I/home/shihao/dlib-19.2 -I/usr/local/lib/python2.7/dist-packages/numpy/core/include/
CFLAGS  := -g -O3 -std=c++11
CXXFLAGS:= $(CFLAGS) -DHAVE_CONFIG_H

#i think you should do anything here
#下面的基本上不需要做任何改动了
.PHONY : everything objs clean veryclean rebuild

everything : $(TARGET)


%.o: %.c
	$(CC) $(CFLAGS) -c $< $(INCLUDE)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< $(INCLUDE)

%.o: %.cc
	$(CC) $(CFLAGS) -c $< $(INCLUDE)

all : $(TARGET)

objs : $(OBJS)

rebuild: veryclean everything

clean :
	rm -fr *.so
	rm -fr *.o

veryclean : clean
	rm -fr $(TARGET)

$(TARGET) : $(OBJS)
	$(CC) $(CXXFLAGS) -o $@ $(OBJS)  $(LDFLAGS) $(LIBS)

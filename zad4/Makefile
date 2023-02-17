CC = g++

CFLAGS = -Wall -std=c++11

OBJECTS = mn.o

TARGET = mn.out

%.o : %.cpp Makefile
	$(CC) $(CFLAGS) -c $< -I/usr/include/eigen3 -o $@

$(TARGET) : $(OBJECTS)
	$(CC) $(CFLAGS) $(LIBS) $(OBJECTS) -o $(TARGET)

.PHONY : clean

clean :
	$(RM) $(TARGET) *.o core
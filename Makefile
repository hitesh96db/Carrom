CC = g++
CFLAGS = -Wall
PROG = a.out

SRCS = basic.cpp

ifeq ($(shell uname),Darwin)
	LIBS = -framework OpenGL -framework GLUT -g
else
	LIBS = -lglut -lGL -lGLU -lm -g
endif

all: $(PROG) 

$(PROG):	$(SRCS)
	$(CC) $(CFLAGS) -o $(PROG) $(SRCS) $(LIBS)

clean:
	rm -f $(PROG)

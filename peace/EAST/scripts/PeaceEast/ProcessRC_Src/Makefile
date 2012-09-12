
CC = g++
DEBUG = -g
OPTS = ${DEBUG} -Wall -Wextra

OS = $(shell uname -s)
ifeq ($(OS),Darwin)
	OPTS += -DOSX
endif

O_FILES = Prim.o ProcessRC.o

ProcessRC: Main.cpp ${O_FILES}
	${CC} -o ProcessRC Main.cpp ${O_FILES} ${OPTS}

Prim.o: Prim.cpp Prim.h
	${CC} -c Prim.cpp ${OPTS}

ProcessRC.o: Node.h Prim.h ProcessRC.cpp ProcessRC.h Prim.cpp
	${CC} -c ProcessRC.cpp ${OPTS}

clean:
	rm *.o ProcessRC

OBJS=nw.o
CFLAGS=-g -I. -Wall -Wextra -Werror -lpthread
#DEFINES=-DTHINK_TIME
BIN=nw
CC=gcc

%.o:%.c
	$(CC) $(CFLAGS) $(DEFINES) -o $@ -c $<

$(BIN): $(OBJS)
	$(CC) $(CFLAGS) $(DEFINES) -o $(BIN) $^

clean:
	rm $(BIN) $(OBJS)

# to run the program, command line "make run"
SEQ1=ACAATCC
SEQ2=AGCATGC
GP=-1
MB=2
MP=-1
WA=0

run:
	./$(BIN) $(SEQ1) $(SEQ2) $(GP) $(MB) $(MP) $(WA)

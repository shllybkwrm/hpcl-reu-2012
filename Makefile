OBJS=main.o
CFLAGS=-g -I. -Wall -Wextra -Werror -O3 -lpthread
#DEFINES=-DTHINK_TIME
BIN=main
CC=gcc

%.o:%.c
	$(CC) $(CFLAGS) $(DEFINES) -o $@ -c $<

$(BIN): $(OBJS)
	$(CC) $(CFLAGS) $(DEFINES) -o $(BIN) $^

clean:
	rm $(BIN) $(OBJS)

# to run the program, command line "make run"

#WIKIPEDIA SW
#SEQ1=ACACACTA
#SEQ2=AGCACACA
#GP=-1
#MB=2
#MP=-1
#WA=1

#BLUE BOOK
SEQ1=ATAGCT
SEQ2=GATATGCA
GP=-2
MB=1
MP=-1
WA=0

run:
	./$(BIN) $(SEQ1) $(SEQ2) $(GP) $(MB) $(MP) $(WA)
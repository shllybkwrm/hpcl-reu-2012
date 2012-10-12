OBJS=main.o
CFLAGS=-g -I. -Wall -Wextra -Werror -lpthread
#DEFINES=-DTHINK_TIME
BIN=main
CC=gcc

%.o:%.c
	$(CC) $(CFLAGS) $(DEFINES) -o $@ -c $<

$(BIN): $(OBJS)
	$(CC) $(CFLAGS) $(DEFINES) -o $(BIN) $^

clean:
	rm $(BIN) $(OBJS)
# unified_codes: BBGS / AG / Goppa research program

CC       ?= gcc
TARGET   ?= unified_codes
SRC      ?= unified_codes.c
CFLAGS   ?= -O3 -march=native -std=gnu11 -Wall -Wextra
OPENMP   ?= -fopenmp
LDFLAGS  ?=
LDLIBS   ?= -lm

.PHONY: all clean check

# `all` is only a Make target. The executable is named `unified_codes`.
all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(OPENMP) $(LDFLAGS) $< -o $@ $(LDLIBS)

match: 
	gcc -O3 -march=native -mtune=native -std=c11 \
  		-Wall -Wextra -o unified_codes unified_codes_matched_errors.c -lm

suggest:
	gcc -O3 -march=native -mtune=native -std=gnu11 \
	  	-Wall -Wextra -o unified_codes unified_codes_suggest_param.c -lm

# Reproducible short regression test: 100 trials, 1000 injected bit errors.
check: $(TARGET)
	OMP_NUM_THREADS=1 ./$(TARGET) \
	  --curve bbgs \
	  --m 6 --q 4 \
	  --a 1 --b 2 \
	  --basis newton \
	  --concat-inner goppa:6:2:64 \
	  --errors 1000 \
	  -n 100

run:
	OMP_NUM_THREADS=1 ./unified_codes
	
clean:
	rm -f $(TARGET)

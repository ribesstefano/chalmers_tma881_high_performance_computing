CXX := gcc
CFLAGS := -std=c11
LDFLAGS = -lm -lpthread

SOURCES = $(wildcard *.c)
EXECS = $(patsubst %.c, %.exe, $(SOURCES))

.PHONY: all

all: newton

newton: newton.c
	$(CXX) $(CFLAGS) -O3 -march=x86-64 $< -o newton $(LDFLAGS)

images: SHELL:=/bin/bash
images:
	for i in *.ppm ; \
	do \
		name=`echo "$$i" | cut -d'.' -f1` ; \
		echo "# ===========================================================" ; \
		echo "# Converting PPM image '$${name}' to PNG." ; \
		echo "# ===========================================================" ; \
		ffmpeg -i "$$i" "`echo "$$i" | cut -d'.' -f1`.png" ; \
	done

compress:
	tar -czvf submission.tar.gz Makefile newton.c

test: compress
	python3 ../reference_code/newton/check_submission.py submission.tar.gz

clean:
	rm -rf *.exe *.o *.s ./newton *.ppm *.png

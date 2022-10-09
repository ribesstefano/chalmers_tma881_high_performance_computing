CXX := gcc
CFLAGS := -std=c11
LDFLAGS = -lm -lpthread

SOURCES = $(wildcard *.c)
EXECS = $(patsubst %.c, %.exe, $(SOURCES))

.PHONY: all

all: newton multi_stage

newton: newton.c
	$(CXX) $(CFLAGS) -O3 -march=x86-64 $< -o newton $(LDFLAGS)

multi_stage: multi_stage.c
	$(CXX) $(CFLAGS) -O3 -march=x86-64 $< -o multi_stage $(LDFLAGS)

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

test: newton
	./newton -t2 -l1000 5
	$(MAKE) images

compress: Makefile newton.c
	tar -czvf submission.tar.gz Makefile newton.c c11threads.h

check: compress
	python3 ../../newton_iteration/check_submission.py submission.tar.gz

clean:
	rm -rf *.exe *.o *.s ./newton ./multi_stage *.ppm *.png extracted/ pictures/

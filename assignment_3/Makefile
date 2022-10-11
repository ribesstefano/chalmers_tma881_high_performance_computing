CXX := gcc
CFLAGS := -std=c11 -O3
LDFLAGS = -lm -lpthread

SOURCES = $(wildcard *.c)
EXECS = $(patsubst %.c, %.exe, $(SOURCES))

.PHONY: all

all: newton multi_stage

newton: newton.c
	$(CXX) $(CFLAGS) -march=native $< -o newton $(LDFLAGS)

multi_stage: multi_stage.c
	$(CXX) $(CFLAGS) -march=native $< -o multi_stage $(LDFLAGS)

images: SHELL:=/bin/bash
images:
	for i in *.ppm ; \
	do \
		outname="`echo "$$i" | cut -d'.' -f1`.jpg" ; \
		echo "# ===========================================================" ; \
		echo "# Generating JPG image '$${outname}'" ; \
		echo "# ===========================================================" ; \
		convert "$$i" "$$outname" ; \
	done

test: newton
	./newton -t1 -l1000 7
	$(MAKE) images

compress: Makefile newton.c
	tar -czvf submission.tar.gz Makefile newton.c c11threads.h

check: compress
	python3 ../../newton_iteration/check_submission.py submission.tar.gz

clean:
	rm -rf *.exe *.o *.s ./newton ./multi_stage *.ppm *.jpg extracted/ pictures/

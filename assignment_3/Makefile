CXX := gcc
CFLAGS := -std=c11 -O3 -march=native
LDFLAGS = -lm -lpthread

SOURCES = $(wildcard *.c)
EXECS = $(patsubst %.c, %.exe, $(SOURCES))

.PHONY: all

all: newton

newton: newton.c
	$(CXX) $(CFLAGS) $< -o newton $(LDFLAGS)

images: SHELL:=/bin/bash
images:
	for i in *.ppm ; \
	do \
		outname="`echo "$$i" | cut -d'.' -f1`.png" ; \
		echo "# ===========================================================" ; \
		echo "# Generating PNG image '$${outname}'" ; \
		echo "# ===========================================================" ; \
		convert "$$i" "$$outname" ; \
	done

test: newton
	./newton -t1 -l1000 7
	$(MAKE) images

compress: Makefile newton.c
	tar -czvf submission.tar.gz Makefile newton.c

check: compress
	# python3 ../../newton_iteration/check_submission.py submission.tar.gz
	python3 /home/hpc2022/newton_iteration/check_submission.py submission.tar.gz

clean:
	rm -rf *.exe *.o *.s ./newton *.ppm *.png extracted/ pictures/

CXX := gcc
CFLAGS := -std=c11
LDFLAGS = -lm # -lOpenCL

SOURCES = $(wildcard *.c)
EXECS = $(patsubst %.c, %.exe, $(SOURCES))

.PHONY: all

all: heat_diffusion

heat_diffusion: heat_diffusion.c
	$(CXX) $(CFLAGS) -O3 -march=native $< -o heat_diffusion $(LDFLAGS)

test: heat_diffusion
	./heat_diffusion -n2 -d0.02

compress: Makefile heat_diffusion.c
	tar -czvf submission.tar.gz Makefile heat_diffusion.c

check: compress
	python3 ../../heat_diffusion_iteration/check_submission.py submission.tar.gz

clean:
	rm -rf *.exe *.o *.s ./heat_diffusion ./multi_stage *.ppm *.png extracted/ pictures/

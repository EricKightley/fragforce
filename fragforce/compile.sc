#!/bin/bash
gcc -c -O3 -fPIC -Wall -Werror -std=c11 force.c
gcc -shared -fPIC -O3 -o force.so force.o


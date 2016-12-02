#!/bin/bash
gcc -c -fPIC -Wall -Werror -std=c11 force.c -O3
gcc -shared -o force.so force.c -fPIC -O3


#!/bin/bash
gcc -c -fPIC force.c -O3
gcc -shared -o force.so force.c -fPIC -O3
 

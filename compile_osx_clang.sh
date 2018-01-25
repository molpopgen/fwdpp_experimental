#!/bin/bash

CXX=clang++ CC=clang CPPFLAGS="-I$HOME/src/fwdpp -I$HOME/anaconda3/include" LDFLAGS="-L$HOME/anaconda3/lib \
-Wl,-rpath,$HOME/anaconda3/lib" make



#!/bin/bash

CPPFLAGS=-I$HOME/src/fwdpp LDFLAGS=-Wl,-rpath,$HOME/anaconda3/lib make

#!/bin/sh

git clone --recursive git@github.com:libopencm3/libopencm3.git
cd libopencm3
git checkout 6763681c260cf280487d70ca0d1996a8b79fff30
make -j8
cd ../

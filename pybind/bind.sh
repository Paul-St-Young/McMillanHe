#!/bin/bash
for name in mmh.h mmh.cpp rng.h; do
  ln -s ../purecpp/$name .
done

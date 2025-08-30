#!/bin/bash

# Set the number of threads
export JULIA_NUM_THREADS=20

# Measure execution time
echo "Starting execution with JULIA_NUM_THREADS=$JULIA_NUM_THREADS"
time julia rl.jl

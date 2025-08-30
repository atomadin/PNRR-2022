#!/bin/bash

# Set the number of threads 32
export JULIA_NUM_THREADS=32

# Measure execution time
echo "Starting execution with JULIA_NUM_THREADS=$JULIA_NUM_THREADS"
time julia rl_gold.jl

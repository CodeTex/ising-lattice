#!/bin/zsh

DATA_DIR="data"
IMG_DIR="img"

FLAGS=(-Og -Wall -fimplicit-none -fcheck=all -fbacktrace)
MOD=(constants.f95 io.f95 shell_prompts.f95)
EXE=main.out

[ ! -d "$DATA_DIR" ] || mkdir -p "$DATA_DIR"
[ ! -d "$IMG_DIR" ] || mkdir -p "$IMG_DIR"

gfortran "${FLAGS[@]}" "${MOD[@]}" main.f95 -o $EXE &&
time ./$EXE &&
gnuplot plot_magnetization.plot &&
gnuplot plot_energy.plot
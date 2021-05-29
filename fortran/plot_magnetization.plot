set terminal png size 1600,1000 font ',13'

set encoding utf8

set output 'img/ising_lattice_magnetization.png'
set autoscale

# multiplot
if (!exists("MP_LEFT"))   MP_LEFT = .06
if (!exists("MP_RIGHT"))  MP_RIGHT = .95
if (!exists("MP_BOTTOM")) MP_BOTTOM = .1
if (!exists("MP_TOP"))    MP_TOP = .9
if (!exists("MP_GAP"))    MP_GAP = .09

set multiplot layout 2,2 columnsfirst \
   title "Monte Carlo simulated Ising lattices" offset 0,-0.5\
   margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP

# Key and axes
set key box right top
set xlabel 'Temperature [J/k_B]' enhanced offset 0,0.5
set yrange [0:1.1]

# input
file1 = 'data/ising_lattice-004.dat'
file2 = 'data/ising_lattice-008.dat'
file3 = 'data/ising_lattice-016.dat'
file4 = 'data/ising_lattice-032.dat'
file5 = 'data/ising_lattice-064.dat'
file6 = 'data/ising_lattice-128.dat'

# colors
col1 = '#7128a1'
col2 = '#e71d36'
col3 = '#f3722c'
col4 = '#f9c74f'
col5 = '#90be6d'
col6 = '#43aa8b'


# plot <|M|>
vcol = 2
set title "<|M|>" offset 0,-0.7
# set ylabel "<|M|>"
plot file1 using 1:vcol title '04x04' w lp pt 19 ps .5 lc rgb col1 lt 1,\
    file2 using 1:vcol title '08x08' w lp pt 19 ps .5 lc rgb col2 lt 1,\
    file3 using 1:vcol title '16x16' w lp pt 19 ps .5 lc rgb col3 lt 1,\
    file4 using 1:vcol title '32x32' w lp pt 19 ps .5 lc rgb col4 lt 1,\
    file5 using 1:vcol title '64x64' w lp pt 19 ps .5 lc rgb col5 lt 1,\
    file6 using 1:vcol title '128x128' w lp pt 19 ps .5 lc rgb col6 lt 1

# plot <|M^2|>
vcol = 3
set title "<|M^2|>" offset 0,-0.7
# set ylabel "<|M^2|>"
plot file1 using 1:vcol title '04x04' w lp pt 19 ps .5 lc rgb col1 lt 2,\
    file2 using 1:vcol title '08x08' w lp pt 19 ps .5 lc rgb col2 lt 2,\
    file3 using 1:vcol title '16x16' w lp pt 19 ps .5 lc rgb col3 lt 2,\
    file4 using 1:vcol title '32x32' w lp pt 19 ps .5 lc rgb col4 lt 2,\
    file5 using 1:vcol title '64x64' w lp pt 19 ps .5 lc rgb col5 lt 2,\
    file6 using 1:vcol title '128x128' w lp pt 19 ps .5 lc rgb col6 lt 2

# plot <|M^4|>
vcol = 4
set title "<|M^4|>" offset 0,-0.7
# set ylabel "<|M^4|>"
plot file1 using 1:vcol title '04x04' w lp pt 19 ps .5 lc rgb col1 lt 3,\
    file2 using 1:vcol title '08x08' w lp pt 19 ps .5 lc rgb col2 lt 3,\
    file3 using 1:vcol title '16x16' w lp pt 19 ps .5 lc rgb col3 lt 3,\
    file4 using 1:vcol title '32x32' w lp pt 19 ps .5 lc rgb col4 lt 3,\
    file5 using 1:vcol title '64x64' w lp pt 19 ps .5 lc rgb col5 lt 3,\
    file6 using 1:vcol title '128x128' w lp pt 19 ps .5 lc rgb col6 lt 3

# plot Binder parameter
vcol = 5
set title "Binder parameter" offset 0,-0.7
# set ylabel "Binder parameter"
plot file1 using 1:vcol title '04x04' w lp pt 19 ps .5 lc rgb col1 lt 4,\
    file2 using 1:vcol title '08x08' w lp pt 19 ps .5 lc rgb col2 lt 4,\
    file3 using 1:vcol title '16x16' w lp pt 19 ps .5 lc rgb col3 lt 4,\
    file4 using 1:vcol title '32x32' w lp pt 19 ps .5 lc rgb col4 lt 4,\
    file5 using 1:vcol title '64x64' w lp pt 19 ps .5 lc rgb col5 lt 4,\
    file6 using 1:vcol title '128x128' w lp pt 19 ps .5 lc rgb col6 lt 4

unset multiplot
set output

set terminal png size 800,600 font ',13'

set encoding utf8

set output 'img/ising_lattice_energy.png'
set title "Monte Carlo simulated Ising lattices"
set autoscale

# Key and axes
set key box right bottom
set xlabel 'Temperature [J/k_B]' enhanced
set ylabel 'Energy [J]'
set yrange [-2500:100]

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


# plot energy
vcol = 6
plot file1 using 1:vcol title '04x04' w lp pt 19 ps .5 lc rgb col1 lt 1,\
    file2 using 1:vcol title '08x08' w lp pt 19 ps .5 lc rgb col2 lt 1,\
    file3 using 1:vcol title '16x16' w lp pt 19 ps .5 lc rgb col3 lt 1,\
    file4 using 1:vcol title '32x32' w lp pt 19 ps .5 lc rgb col4 lt 1,\
    file5 using 1:vcol title '64x64' w lp pt 19 ps .5 lc rgb col5 lt 1,\
    file6 using 1:vcol title '128x128' w lp pt 19 ps .5 lc rgb col6 lt 1

set output

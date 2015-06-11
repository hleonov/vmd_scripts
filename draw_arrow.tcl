proc vmd_draw_arrow {mol start end} {
    # an arrow is made of a cylinder and a cone
    set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
    graphics $mol cylinder $start $middle radius 0.35
    graphics $mol cone $middle $end radius 0.55
}

# first specify the color with
# graphics 0 color 3
# use the script top generate the axis and extesion to the ends
# helix_axis.pl d133.pdb 130 125 4 2.0 1.0
# draw arrow {0 0 1} {2 2 2}
# use this in vmd
# graphics 0 list
# graphics 0 delete 2
# graphics 0 delete all

#!/usr/bin/env python3
#This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
#Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
"""
Created on Tue Sep 13 10:33:21 2022

@author: Gururaj Wagle

Create a input.wlgrid file for a custom grid: a log-spaced base grid with a
higher-resolution log-spaced insert. Same can be used to create an
input.fluxwl file (via --fname-out).
"""
import argparse
import numpy as np


def main():
    parser = argparse.ArgumentParser(
        description='Create an input.wlgrid file for a custom wavelength'
        ' grid: log-spaced base grid with a higher-resolution insert.')
    parser.add_argument('--fname-out', type=str, default='input.wlgrid',
                        help='Output file name (default: input.wlgrid)')
    parser.add_argument('--show-plots', action='store_true',
                        help='Display a plot of the grid spacing')
    args = parser.parse_args()

    # construct grid
    ## set min, max values and number of bins
    ng = 50
    wlmin = 1e-8
    wlmax = 50e-5
    wlgrid = np.zeros(ng)
    for i in np.arange(ng):
        wlgrid[i] = wlmin*(wlmax/wlmin)**(float(i)/float(ng))

    # split the base grid around the high-resolution window
    wlgrid1 = wlgrid[wlgrid < 3e-5]
    wlgrid2 = wlgrid[wlgrid > 10e-5]

    wlgrid_res = np.logspace(np.log10(3e-5), np.log10(10e-5), 1500)

    wlgrid_new = np.concatenate((wlgrid1, wlgrid_res, wlgrid2))

    if args.show_plots:
        from matplotlib import pyplot as plt
        fig99 = plt.figure(num=99)
        fig99.clear()
        ax99 = fig99.add_subplot(1, 1, 1)
        ax99.plot(wlgrid[:-1], np.diff(wlgrid))
        ax99.plot(wlgrid_new[:-1], np.diff(wlgrid_new))
        fig99.tight_layout()
        plt.show()

    # write grid values
    with open(args.fname_out, 'w') as f:
        print("File name: ", f.name)
        f.write(str(len(wlgrid_new)) + " ")
        for wlval in wlgrid_new:
            f.write("%.6e " % wlval)
        f.write("#")


if __name__ == '__main__':
    main()

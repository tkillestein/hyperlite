#!/usr/bin/env python3
#This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
#Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.

#--------------------------------------------------------------------------
# Post-process flux data (output.flx_luminos) with a Savitzky-Golay filter.
# This yields a flux file with the same format as the input, but smoothed
# over wavelength (the number of wavelength groups is unchanged).
#--------------------------------------------------------------------------

#-- libraries
import numpy as np
import argparse


def main():
    #-- build parse field
    parser = argparse.ArgumentParser(description='Post-process flux data'
                                     ' with a Savitzky-Golay filter.')
    parser.add_argument('flux', type=str, help='File of flux data')
    parser.add_argument('flux_grid', type=str, help='File of flux grid data')
    parser.add_argument('--nfilter', type=int, default=5,
                        help='Number of filter points')
    parser.add_argument('--fname_out', type=str, help='Optional name for'
                        ' output file.')

    #-- parse the arguments
    args = parser.parse_args()

    #-- load the flux data (atleast_2d: a snapshot run writes a single row)
    flux = np.atleast_2d(np.loadtxt(args.flux))

    #-- number of time steps
    nt = np.size(flux, 0)

    #-- get wavelength grid from flux grid
    with open(args.flux_grid) as wlfile:
        wlfile.readline()
        wlline = wlfile.readline().split()
    wl = np.array([float(wlval) for wlval in wlline])

    #-- ensure a sufficient number of groups
    nfilter = args.nfilter
    ng = np.size(wl) - 1
    if ng < nfilter:
        print('WARNING: insufficient number of groups (< nfilter)')

    #-- group centers are used as independent variable
    wlc = .5*(wl[1:] + wl[:ng])
    dwl = wl[1:] - wl[:ng]
    #-- regress over per wl values
    flux_spec = flux/dwl

    #-- initialize new smooth flux matrix
    smooth_flux = flux_spec.copy()

    #-- smooth spectrum for each time step
    nfilr = (nfilter - 1)//2
    for it in range(nt):

        #-- do not smooth on boundary wl values
        for ig in range(nfilr, ng - nfilr):

            #-- approximate inverse matrix for a cubic regression
            wl3mat = np.vander(wlc[ig - nfilr:ig + nfilr], 4)
            wlmat = np.dot(wl3mat.T, wl3mat)
            wlmat = np.linalg.inv(wlmat)
            wlmat = np.dot(wlmat, wl3mat.T)

            #-- calculate SG coefficients
            abar = np.dot(wlmat, flux_spec[it, ig - nfilr:ig + nfilr])

            #-- smooth the flux at this point
            smflux = np.dot(abar, wl3mat[nfilr, :])
            if smflux > 0.0:
                smooth_flux[it, ig] = smflux
            else:
                smooth_flux[it, ig] = flux_spec[it, ig]

        #-- now multiply by dwl get to erg/s and renormalize
        smooth_flux[it, :] *= dwl
        smooth_flux[it, :] *= np.sum(flux[it, :])/np.sum(smooth_flux[it, :])

    #-- save the smoothed flux
    if args.fname_out is None:
        np.savetxt('output.flx_luminos_smoothed', smooth_flux,
                   fmt='% .4e', delimiter=' ')
    else:
        np.savetxt(args.fname_out, smooth_flux, fmt='% .4e', delimiter=' ')


if __name__ == '__main__':
    main()

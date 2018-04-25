#!/usr/local/env python

import netCDF4
import numpy as np

def ap2ep(Au, PHIu, Av, PHIv, plot_demo=False):
    """
    Convert tidal amplitude and phase lag (ap-) parameters into tidal ellipse
    (ep-) parameters. Please refer to ep2ap for its inverse function.

    Usage:

    SEMA, ECC, INC, PHA, w = ap2ep(Au, PHIu, Av, PHIv, plot_demo=False)

    Where:

        Au, PHIu, Av, PHIv are the amplitudes and phase lags (in degrees) of
        u- and v- tidal current components. They can be vectors or
        matrices or multidimensional arrays.

        plot_demo is an optional argument, when it is supplied as an array
        of indices, say [i j k l], the program will plot an ellipse
        corresponding to Au[i, j, k, l], PHIu[i, j, k, l], Av[i, j, k, l], and
        PHIv[i, j, k, l]. Defaults to False (i.e. no plot).

        Any number of dimensions are allowed as long as your computer
        resource can handle.

        SEMA: Semi-major axes, or the maximum speed.

        ECC:  Eccentricity, the ratio of semi-minor axis over the semi-major
              axis; its negative value indicates that the ellipse is traversed
              in clockwise direction.

        INC:  Inclination, the angles (in degrees) between the semi-major axes
              and u-axis.

        PHA:  Phase angles, the time (in angles and in degrees) when the tidal
              currents reach their maximum speeds,  (i.e.  PHA=omega*tmax).

              These four ep-parameters will have the same dimensionality (i.e.,
              vectors, or matrices) as the input ap-parameters.

        w:    A matrix whose rows allow for plotting ellipses and whose columns
              are for different ellipses corresponding columnwise to SEMA. For
              example, plot(np.real(w[0, :]), np.imag(w[0, :])) will let you
              see the first ellipse. You may need to use squeeze function when
              w is a more than two dimensional array. See example.py.

    Document:   tidal_ellipse.ps

    Revisions: May  2002, by Zhigang Xu,  --- adopting Foreman's northern semi
    major axis convention.

    For a given ellipse, its semi-major axis is undetermined by 180. If we
    borrow Foreman's terminology to call a semi major axis whose direction lies
    in a range of [0, 180) as the northern semi-major axis and otherwise as a
    southern semi major axis, one has freedom to pick up either northern or
    southern one as the semi major axis without affecting anything else.
    Foreman (1977) resolves the ambiguity by always taking the northern one as
    the semi-major axis. This revision is made to adopt Foreman's convention.
    Note the definition of the phase, PHA, is still defined as the angle
    between the initial current vector, but when converted into the maximum
    current time, it may not give the time when the maximum current first
    happens; it may give the second time that the current reaches the maximum
    (obviously, the 1st and 2nd maximum current times are half tidal period
    apart) depending on where the initial current vector happen to be and its
    rotating sense.

    Version 2, May 2002

    Converted to Python by Pierre Cazenave, October 2012.

    Authorship Copyright:

       The author retains the copyright of this program, while  you are welcome
    to use and distribute it as long as you credit the author properly and
    respect the program name itself. Particularly, you are expected to retain
    the original author's name in this original version or any of its modified
    version that you might make. You are also expected not to essentially
    change the name of the programs except for adding possible extension for
    your own version you might create, e.g. ap2ep_xx is acceptable.  Any
    suggestions are welcome and enjoy my program(s)!


    Author Info:
    _______________________________________________________________________
      Zhigang Xu, Ph.D.
      (pronounced as Tsi Gahng Hsu)
      Research Scientist
      Coastal Circulation
      Bedford Institute of Oceanography
      1 Challenge Dr.
      P.O. Box 1006                    Phone  (902) 426-2307 (o)
      Dartmouth, Nova Scotia           Fax    (902) 426-7827
      CANADA B2Y 4A2                   email xuz@dfo-mpo.gc.ca
    _______________________________________________________________________

    Release Date: Nov. 2000, Revised on May. 2002 to adopt Foreman's northern
    semi major axis convention.

    """

    # Assume the input phase lags are in degrees and convert them in radians.
    PHIu = PHIu / 180 * np.pi
    PHIv = PHIv / 180 * np.pi

    # Make complex amplitudes for u and v
    i = 1j
    u = Au * np.exp(-i * PHIu)
    v = Av * np.exp(-i * PHIv)

    # Calculate complex radius of anticlockwise and clockwise circles:
    wp = (u + i * v) / 2           # for anticlockwise circles
    wm = np.conj(u - i * v) / 2    # for clockwise circles
    # and their amplitudes and angles
    Wp = np.abs(wp)
    Wm = np.abs(wm)
    THETAp = np.angle(wp)
    THETAm = np.angle(wm)

    # calculate ep-parameters (ellipse parameters)
    SEMA = Wp + Wm                 # Semi Major Axis, or maximum speed
    SEMI = Wp - Wm                 # Semi Minor Axis, or minimum speed
    ECC = SEMI / SEMA              # Eccentricity

    PHA = (THETAm - THETAp) / 2    # Phase angle, the time (in angle) when
                                   # the velocity reaches the maximum
    INC = (THETAm + THETAp) / 2    # Inclination, the angle between the
                                   # semi major axis and x-axis (or u-axis).

    # convert to degrees for output
    PHA = PHA / np.pi*180
    INC = INC / np.pi*180
    THETAp = THETAp / np.pi*180
    THETAm = THETAm / np.pi*180

    # map the resultant angles to the range of [0, 360].
    PHA = np.mod(PHA + 360, 360)
    INC = np.mod(INC + 360, 360)

    # Mar. 2, 2002 Revision by Zhigang Xu    (REVISION_1)
    # Change the southern major axes to northern major axes to conform the tidal
    # analysis convention  (cf. Foreman, 1977, p. 13, Manual For Tidal Currents
    # Analysis Prediction, available in www.ios.bc.ca/ios/osap/people/foreman.htm)
    k = np.fix(INC / 180)
    INC = INC - k * 180
    PHA = PHA + k * 180
    PHA = np.mod(PHA, 360)


    ndot = np.prod(np.shape(SEMA))
    dot = 2 * np.pi / ndot
    ot = np.arange(0, 2 * np.pi, dot)
    w = wp.flatten() * np.exp(i * ot) + wm.flatten() * np.exp(-i * ot)
    w = np.reshape(w, np.shape(wp))

    return SEMA, ECC, INC, PHA, w

 

# input tidal constituents url
url = 'http://gamone.whoi.edu/thredds/dodsC/usgs/vault0/models/tides/ec2015/f54.ncml'


toexclude = ['UAmp', 'UPha','VAmp','VPha']
var_atts={'mesh':'adcirc_mesh', 'location':'node', 'coordinates':'lon lat', 'units':'m/s'}

with netCDF4.Dataset(url) as src, netCDF4.Dataset("m2_maj.nc", "w") as dst:
    # copy global attributes all at once via dictionary
    dst.setncatts(src.__dict__)
    # copy dimensions
    for name, dimension in src.dimensions.items():
        dst.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))
    # copy all file data except for the excluded
    for name, variable in src.variables.items():
        if name not in toexclude:
            x = dst.createVariable(name, variable.datatype, variable.dimensions)
            dst[name][:] = src[name][:]
            # copy variable attributes all at once via dictionary
            dst[name].setncatts(src[name].__dict__)
    if True:
        m2 = dst.createVariable('M2_UMaj', src['depth'].datatype, src['depth'].dimensions)
        var_atts['long_name']='M2 UMajor'
        m2.setncatts(var_atts)
       
        Au = src['UAmp'][0,:]
        PHIu = src['UPha'][0,:]
        Av = src['VAmp'][0,:]
        PHIv = src['VPha'][0,:]
 
        SEMA, ECC, INC, PHA, w = ap2ep(Au, PHIu, Av, PHIv) 

        m2[:] = SEMA

        dst['depth'][:] = src['depth'][:] * -1.0



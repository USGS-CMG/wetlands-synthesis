#!/usr/bin/env python

# Extract M2, S2, N2, O1, K1 tidal amplitudes from Amp array and write to new file

import netCDF4

toexclude = ['Amp', 'Pha','Re','Im']

url = 'http://gamone.whoi.edu/thredds/dodsC/usgs/vault0/models/tides/ec2015/f53.ncml'

var_atts={'mesh':'adcirc_mesh', 'location':'node', 'coordinates':'lon lat', 'units':'m'}

# this part is generic to copy all but specified variables to new file

with netCDF4.Dataset(url) as src, netCDF4.Dataset("out.nc", "w") as dst:
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

	if True:     #this part is custom to add/modify the generic part
        m2 = dst.createVariable('M2_Amp', src['depth'].datatype, src['depth'].dimensions)
        var_atts['long_name']='M2 Amplitude'
        m2.setncatts(var_atts)
        m2[:] = src['Amp'][0,:]
        
        n2 = dst.createVariable('N2_Amp', src['depth'].datatype, src['depth'].dimensions)
        var_atts['long_name']='N2 Amplitude'
        n2.setncatts(var_atts)
        n2[:] = src['Amp'][1,:]
        
        s2 = dst.createVariable('S2_Amp', src['depth'].datatype, src['depth'].dimensions)
        var_atts['long_name']='S2 Amplitude'
        s2.setncatts(var_atts)
        s2[:] = src['Amp'][2,:]

        o1 = dst.createVariable('O1_Amp', src['depth'].datatype, src['depth'].dimensions)
        var_atts['long_name']='O1 Amplitude'
        o1.setncatts(var_atts)
        o1[:] = src['Amp'][3,:]
        
        k1 = dst.createVariable('K1_Amp', src['depth'].datatype, src['depth'].dimensions)
        var_atts['long_name']='K1 Amplitude'
        k1.setncatts(var_atts)
        k1[:] = src['Amp'][4,:]

        dst['depth'][:] = src['depth'][:] * -1.0
 

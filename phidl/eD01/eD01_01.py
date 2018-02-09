#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 16:30:04 2017

@author: smb2
"""

from phidl import Device, Layer, quickplot2
import phidl.geometry as pg
import phidl.routing as pr
import basic_photonics as bp
import gdspy
import numpy as np

#==============================================================================
#Define the layers used in the process
#==============================================================================

layers = {
        'si_partial' : Layer(gds_layer = 0, gds_datatype = 0, description = 'partial Si etch', color = 'gray'),
        'p' : Layer(gds_layer = 1, gds_datatype = 0, description = 'p dopants', color = 'black'),
        'n' : Layer(gds_layer = 2, gds_datatype = 0, description = 'n dopants', color = 'green'),
        'pad'  : Layer(gds_layer = 3, gds_datatype = 0, description = 'pads', color = 'gold', alpha = 0.2),
        'nitride_open' : Layer(gds_layer = 4, gds_datatype = 0, description = 'nitride opening', color = 'red'),
        'void' : Layer(gds_layer = 98, gds_datatype = 0, description = 'for blank areas that I want to align to', color = 'gray', alpha = 0.05),
        'chip' : Layer(gds_layer = 99, gds_datatype = 0, description = 'chip area', color = 'gray', alpha = 0.01)
         }

#==============================================================================
# Special functions for this layout
#==============================================================================

def resistance_test(length = 20.,
width = 4.,
nitride_offset = 0.2,
dopant = 'p',
pad_offset = 100.,
pad_overlap = 4.,
small_pad_length = 20.,
etch_offset = 25.,
layers = layers,
si_etch = 1,
pad_device = Device()):
    
    D = Device()
    SmallPad = pg.compass(size = (small_pad_length, width), layer = layers['pad'])
    NitrideEt = pg.compass(size = (pad_overlap - nitride_offset, width-nitride_offset), layer = layers['nitride_open'])
    
    pad1 = D.add_ref(pad_device)
    smallpad1 = D.add_ref(SmallPad)
    smallpad2 = D.add_ref(SmallPad)
    nitride_etch1 = D.add_ref(NitrideEt)
    nitride_etch2 = D.add_ref(NitrideEt)
    
    pad2 = D.add_ref(pad_device)
    
    if dopant == 'p':
        layer_dopant = layers['p']
    elif dopant == 'n':
        layer_dopant = layers['n']
       
    Res_test = pg.compass(size = (length, width), layer = layer_dopant)
    res_test = D.add_ref(Res_test)
    pad1.xmin = res_test.xmax + pad_offset
    pad2.xmax = res_test.xmin - pad_offset
    smallpad1.xmin = res_test.xmax - pad_overlap
    smallpad2.xmax = res_test.xmin + pad_overlap
    
    nitride_etch1.xmax = smallpad2.xmax-nitride_offset
    nitride_etch2.xmin = smallpad1.xmin+nitride_offset
    
    D.add_ref(pr.route_basic(port1 = pad1.ports['W'], port2 =smallpad1.ports['E'], layer = layers['pad'], path_type = 'straight'))
    D.add_ref(pr.route_basic(port1 = pad2.ports['E'], port2 =smallpad2.ports['W'], layer = layers['pad'], path_type = 'straight')) 
    
    if si_etch == 1:
        SiEtch = pg.compass(size = D.size + (etch_offset, etch_offset), layer = layers['si_partial'])
        sietch = D.add_ref(SiEtch)
    
    return D

def resistance_test_control(length = 20.,
width = 4.,
nitride_offset = 0.2,
dopant = 'p',
pad_offset = 100.,
pad_overlap = 4.,
small_pad_length = 20.,
etch_offset = 25.,
layers = layers,
si_etch = 1,
pad_device = Device()):
    
    D = Device()
    SmallPad = pg.compass(size = (small_pad_length, width), layer = layers['pad'])
    NitrideEt = pg.compass(size = (pad_overlap - nitride_offset, width-nitride_offset), layer = layers['nitride_open'])
        
    pad1 = D.add_ref(pad_device)
    smallpad1 = D.add_ref(SmallPad)
    smallpad2 = D.add_ref(SmallPad)
    nitride_etch1 = D.add_ref(NitrideEt)
    nitride_etch2 = D.add_ref(NitrideEt)
    
    pad2 = D.add_ref(pad_device)
    
    if dopant == 'p':
        layer_dopant = layers['p']
    elif dopant == 'n':
        layer_dopant = layers['n']
        
    Doping_pad = pg.compass(size = (small_pad_length, width), layer = layer_dopant)
    doping_pad1 = D.add_ref(Doping_pad)
    doping_pad2 = D.add_ref(Doping_pad)
    
    Res_test = pg.compass(size = (length, width), layer = layers['void'])
    res_test = D.add_ref(Res_test)
    pad1.xmin = res_test.xmax + pad_offset
    pad2.xmax = res_test.xmin - pad_offset
    smallpad1.xmin = res_test.xmax - pad_overlap
    smallpad2.xmax = res_test.xmin + pad_overlap
    doping_pad1.center = smallpad1.center
    doping_pad2.center = smallpad2.center
    
    
    nitride_etch1.xmax = smallpad2.xmax-nitride_offset
    nitride_etch2.xmin = smallpad1.xmin+nitride_offset
    
    
    
    D.add_ref(pr.route_basic(port1 = pad1.ports['W'], port2 =smallpad1.ports['E'], layer = layers['pad'], path_type = 'straight'))
    D.add_ref(pr.route_basic(port1 = pad2.ports['E'], port2 =smallpad2.ports['W'], layer = layers['pad'], path_type = 'straight')) 
    
    if si_etch == 1:
        SiEtch = pg.compass(size = D.size + (etch_offset, etch_offset), layer = layers['si_partial'])
        sietch = D.add_ref(SiEtch)
    
    return D

def resistance_test_nodoping(length = 20.,
width = 4.,
nitride_offset = 0.2,
pad_offset = 100.,
pad_overlap = 4.,
small_pad_length = 20.,
etch_offset = 25.,
layers = layers,
si_etch = 1,
dopant = 'p',
pad_device = Device()):
    
    D = Device()
    SmallPad = pg.compass(size = (small_pad_length, width), layer = layers['pad'])
    NitrideEt = pg.compass(size = (pad_overlap - nitride_offset, width-nitride_offset), layer = layers['nitride_open'])
    
    pad1 = D.add_ref(pad_device)
    smallpad1 = D.add_ref(SmallPad)
    smallpad2 = D.add_ref(SmallPad)
    nitride_etch1 = D.add_ref(NitrideEt)
    nitride_etch2 = D.add_ref(NitrideEt)
    
    pad2 = D.add_ref(pad_device)
    
       
    Res_test = pg.compass(size = (length, width), layer = layers['void'])
    res_test = D.add_ref(Res_test)
    pad1.xmin = res_test.xmax + pad_offset
    pad2.xmax = res_test.xmin - pad_offset
    smallpad1.xmin = res_test.xmax - pad_overlap
    smallpad2.xmax = res_test.xmin + pad_overlap
    
    nitride_etch1.xmax = smallpad2.xmax-nitride_offset
    nitride_etch2.xmin = smallpad1.xmin+nitride_offset
    
    D.add_ref(pr.route_basic(port1 = pad1.ports['W'], port2 =smallpad1.ports['E'], layer = layers['pad'], path_type = 'straight'))
    D.add_ref(pr.route_basic(port1 = pad2.ports['E'], port2 =smallpad2.ports['W'], layer = layers['pad'], path_type = 'straight')) 
    
    if si_etch == 1:
        SiEtch = pg.compass(size = D.size + (etch_offset, etch_offset), layer = layers['si_partial'])
        sietch = D.add_ref(SiEtch)
    
    return D



def chip_corners(length = 1000.,
width = 25.,
chipsize = (1000, 1000),
layers = layers):
   
    D = Device()
 
    length = 1000
    width = 100
    
    D = Device()
    Chpmrk = Device()
    
    corner2 = np.array([5000,-5000])
    corner3 = np.array([-5000,-5000])
    corner4 = np.array([-5000,5000])
    corner1 = np.array([5000,5000])
    
    mark = pg.compass(size = (length, width), layer = layers['pad'])
    mark1 = Chpmrk.add_ref(mark)
    mark2 = Chpmrk.add_ref(mark)
    mark2.rotate(angle = 90)
    mark1.ymax = mark2.ymax
    mark1.xmax = mark2.xmax
    
    chpmrk1 = D.add_ref(Chpmrk)
    chpmrk1.ymax = corner1[1]
    chpmrk1.xmax = corner1[0]
    
    chpmrk2 = D.add_ref(Chpmrk)
    chpmrk2.rotate(-90)
    chpmrk2.ymin = corner2[1]
    chpmrk2.xmax = corner2[0]
    
    chpmrk3 = D.add_ref(Chpmrk)
    chpmrk3.rotate(180)
    chpmrk3.ymin = corner3[1]
    chpmrk3.xmin = corner3[0]
    
    chpmrk4 = D.add_ref(Chpmrk)
    chpmrk4.rotate(90)
    chpmrk4.ymax = corner4[1]
    chpmrk4.xmin = corner4[0]
    
    return D

#==============================================================================
# Define some parameters
#==============================================================================

# Taper constants
resistance_test_parameters = {
    'length' : 20.,
    'width' : 4.,
    'nitride_offset' : 0.2,
    'dopant' : 'p',
    'pad_offset' : 100,
    'pad_overlap' : 10,
    'small_pad_length' : 20,
    'etch_offset' : 25,
    'layers' : layers,
    'si_etch' : 1,
            }

# Taper constants
pad_parameters = {
            'size' : (150, 300),
            'layer' : layers['pad'],
            }

#==============================================================================
# Make the basic structures
#==============================================================================

#The main device
D = Device()

# make a pad device for the devices. We can pass this structure to as many 
# devices as we like
Pad = Device(pg.compass, config = pad_parameters)
#==============================================================================
# Make the more complex structures 
#==============================================================================

#==============================================================================
# here is where we plot and make device. Could also make arrays here 
#==============================================================================

# define some more shared parameters

gapx = 800
gapy = 500
width_list = np.array([1, 10, 25, 50, 75, 100])
width_nom = 50
length_list = np.array([100, 200, 300, 500, 750, 1000])
xvec = (np.append(0, length_list) + np.append(length_list, 0))/2+gapx
yvec = ([500, 500, 500, 500, 500, 500])
#These are the devices that have doping and an etch

UD = Device()

for ii in range(0, np.size(width_list)):
    for jj in range(0, np.size(length_list)):
        RT = Device(resistance_test, config = resistance_test_parameters, pad_device = Pad, width = width_list[ii], length = length_list[jj])
        RT.x = xvec[jj]+np.sum(xvec[0:jj])
        RT.y = yvec[ii]+np.sum(yvec[0:ii])
        UD.add_ref(RT)
        
UD.center = (0,2050)

# devices that have doping under the contacts, but not in between for controls, etch

Control = Device()

for jj in range(0, np.size(length_list)):
     RT = Device(resistance_test_control, config = resistance_test_parameters, pad_device = Pad, width = width_nom, length = length_list[jj])
     RT.center = (xvec[jj]+np.sum(xvec[0:jj]), yvec[ii]+np.sum(yvec[0:ii]))
     Control.add_ref(RT)
Control.x = 0
Control.ymax = UD.bbox[1,1]+gapy

# devices that have no doping at all, etch

Nodoping = Device()

for jj in range(0, np.size(length_list)):
     RT = Device(resistance_test_nodoping, config = resistance_test_parameters, pad_device = Pad, width = width_nom, length = length_list[jj])
     RT.center = (xvec[jj]+np.sum(xvec[0:jj]), yvec[ii]+np.sum(yvec[0:ii]))
     Nodoping.add_ref(RT)
Nodoping.x = 0
Nodoping.ymax = Control.bbox[1,1]+gapy

# These are the devices that have doping and no etch

LD = Device()

for ii in range(0, np.size(width_list)):
    for jj in range(0, np.size(length_list)):
        RT = Device(resistance_test, config = resistance_test_parameters, pad_device = Pad, width = width_list[ii], length = length_list[jj], si_etch = 0)
        RT.center = (xvec[jj]+np.sum(xvec[0:jj]), yvec[ii]+np.sum(yvec[0:ii]))
        LD.add_ref(RT)
LD.center = (0,-2050)

label = D.add_ref(pg.text(text = 'SMB 20171018-eD01', size = 100, position=(0, 0), justify = 'center', layer = layers['pad']))

LControl = Device()

for jj in range(0, np.size(length_list)):
     RT = Device(resistance_test_control, config = resistance_test_parameters, pad_device = Pad, width = width_nom, length = length_list[jj], si_etch = 0)
     RT.center = (xvec[jj]+np.sum(xvec[0:jj]), yvec[ii]+np.sum(yvec[0:ii]))
     LControl.add_ref(RT)
LControl.x = 0
LControl.ymin = LD.bbox[0,1]-gapy

# devices that have no doping at all, etch

LNodoping = Device()

for jj in range(0, np.size(length_list)):
     RT = Device(resistance_test_nodoping, config = resistance_test_parameters, pad_device = Pad, width = width_nom, length = length_list[jj], si_etch = 0)
     RT.center = (xvec[jj]+np.sum(xvec[0:jj]), yvec[ii]+np.sum(yvec[0:ii]))
     LNodoping.add_ref(RT)
LNodoping.x = 0
LNodoping.ymin = LControl.bbox[0,1]-gapy

# Add references to all
D.add_ref(UD)
D.add_ref(LD)
D.add_ref(Control)
D.add_ref(Nodoping)
D.add_ref(LControl)
D.add_ref(LNodoping)

# Chip
chip = D.add_ref(pg.rectangle(size = (10000,10000), layer = 99))
chip.center = (0,0)
D.add_ref(chip_corners(length = 1000.,width = 25.,chipsize = (1000, 1000),
layers = layers))

#fill = pg.fill_rectangle(D, fill_size = (50.0,50.0), avoid_layers = [layers['si_partial'], layers['p'], layers['n'], layers['nitride_open'], layers['void']], include_layers = [layers['pad']],
#                    margin = 0, fill_layers = [layers['pad']], 
#                   fill_densities = [1], fill_inverted = [False], bbox = None)
#D.add_ref(fill)

# plot the structure
# quickplot2(D)

# make into a gds file
D.write_gds('eD_01.gds')


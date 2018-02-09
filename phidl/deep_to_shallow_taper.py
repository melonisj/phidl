#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 15:33:04 2017

@author: smb2
"""

from phidl import Device, Layer, quickplot2
import phidl.geometry as pg
import phidl.routing as pr
import basic_photonics as bp
import gdspy


#==============================================================================
# Special functions for this layout
#==============================================================================

def deep_shallow_deep_wg(Taper = Device()):
    D = Device()
    taper = D.add_ref(Taper)
    taper2 = D.add_ref(Taper)
    taper2.rotate(angle = 180)
    taper2.center = taper.center
    taper2.xmax = taper.xmin
    
    D.add_port(port = taper.ports[3], name = 1)
    D.add_port(port = taper2.ports[3], name = 2)
    
    return D

def multi_taper(ntapers = 3, length_deep_wg = 100, ring_gap = 0.2, Ring = Device(), TwoTaper = Device(), Grating = Device()):
    
    D = Device()
    length= float(length_deep_wg)/(ntapers+1)
   
    Multi = Device()
    ring = D.add_ref(Ring)
    grating = D.add_ref(Grating)             
    grating2 = D.add_ref(Grating)
    grating2.rotate(angle = 180)       
 
    if ntapers > 0:        
        for i in range(0, ntapers):
            taper = Multi.add_ref(TwoTaper)
            taper.x += i*(length+taper.size[0])
            p = Multi.add_port(name = i+1, midpoint = [taper.xmax+length, 0], width = taper.ports[1].width, orientation = 180)
            wg = Multi.add_ref(pr.route_basic(port1 = taper.ports[1], port2 = p))
                
        Multi.add_port(name = i+2, midpoint = [Multi.xmin, Multi.y], width = taper.ports[1].width, orientation = 180)
        multitaper = D.add_ref(Multi)   
        grating.xmax = multitaper.xmin-length
        end_wg = D.add_ref(pr.route_basic(port1 = Multi.ports[i+2], port2 = grating.ports[1]))
        ring.x = Multi.xmin - length/2
        ring.ymin = end_wg.ymax + ring_gap
    else:
        Multi.add_ref(pg.compass(size = (length, TwoTaper.ports[1].width), layer = 0))
        multitaper = D.add_ref(Multi)
        end_wg = multitaper
        grating.xmax = multitaper.xmin
        ring.x = Multi.xmin + length/2
        ring.ymin = end_wg.ymax + ring_gap
    
    grating2.xmin = multitaper.xmax

    return D


#==============================================================================
# Define the layers
#==============================================================================

layers = {
        'full_wg_small' : Layer(gds_layer = 0, gds_datatype = 0, description = 'fully etched Si waveguides small features', color = 'gray'),
        'full_wg_big' : Layer(gds_layer = 0, gds_datatype = 1, description = 'fully etched Si large features', color = 'gold'),
        'shallow_wg_small'  : Layer(gds_layer = 1, gds_datatype = 0, description = 'shallow etched Si waveguide small features', color = 'lightblue', alpha = 0.2),
        'shallow_wg_big' : Layer(gds_layer = 1, gds_datatype = 1, description = 'shallow etched Si waveguide large features', color = 'blue'),
        'chip' : Layer(gds_layer = 99, gds_datatype = 0, description = 'chip area', color = 'gray', alpha = 0.01)
         }

#==============================================================================
# Initialize device
#==============================================================================

D = Device()
D.add_ref(bp.chip_corners(chip_size = 1000, layer = layers['full_wg_big']))

#==============================================================================
# Define the constants
#==============================================================================

# Taper constants
taper_parameters = {
            'width_full_etch' : 3.0,
            'width_shallow_wg' : 0.7,
            'width_full_wg' : 0.4,
            'length_shallow_wg' : 10.0,
            'length_taper' : 2.0,
            'layers' : layers
            }

#Grating constants
grating_parameters = {'num_periods' : 20,
            'period' : 0.75,
            'fill_factor' : 0.5,
            'width_grating' : 5.0,
            'length_taper' : 10.0,
            'width' : 0.4,
            'partial_etch' : True
            }

#Ring constants
ring_parameters = {'radius' : 3.0,
                   'width' : 0.5,
                   'layer': 0
                   }

#Multi-taper constants
ntapers = [0, 1,2,3,4,5]
length_deep_wg = 100.0
length_shallow_wg = 100.0
ring_gap = 0.2

#==============================================================================
# Make the basic structures
#==============================================================================

Grating = Device(bp.grating, config = grating_parameters)
Ring = Device(pg.ring, config = ring_parameters)
Taper_nom = Device(bp.wg_fulltaper, config = taper_parameters) 
TwoTaper_nom = deep_shallow_deep_wg(Taper = Taper_nom)

#==============================================================================
# Make the more complex structures and draw layout
#==============================================================================



disty = (TwoTaper_nom.size[1]+Ring.size[1])*1.5

for i in ntapers:
    if i > 0:
        Taper = Device(bp.wg_fulltaper, config = taper_parameters, length_shallow_wg = length_shallow_wg/i) 
        TwoTaper = deep_shallow_deep_wg(Taper = Taper)
    else:
        Taper = Taper_nom
        TwoTaper = TwoTaper_nom
        
    multi = multi_taper(ntapers = i, length_deep_wg = 100, ring_gap = 0.2, Ring = Ring, TwoTaper = TwoTaper, Grating = Grating)
    multi.y += i*disty
    multi.xmin = 0 
    D.add_ref(multi)

    
quickplot2(D)

D.write_gds('aS01.gds')

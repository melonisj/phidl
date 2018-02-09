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

#==============================================================================
#Define the layers used in the process
#==============================================================================

layers = {
        'wg_gaas' : Layer(gds_layer = 0, gds_datatype = 0, description = 'fully etched Si waveguides small features', color = 'gray'),
        'wg_algaas' : Layer(gds_layer = 1, gds_datatype = 0, description = 'fully etched Si large features', color = 'black'),
        'nw' : Layer(gds_layer = 2, gds_datatype = 0, description = 'fully etched Si large features', color = 'green'),
        'pad'  : Layer(gds_layer = 3, gds_datatype = 0, description = 'shallow etched Si waveguide small features', color = 'gold', alpha = 0.2),
        'padopen_clad' : Layer(gds_layer = 4, gds_datatype = 0, description = 'shallow etched Si waveguide large features', color = 'red'),
        'padopen_spacer' : Layer(gds_layer = 5, gds_datatype = 0, description = 'shallow etched Si waveguide large features', color = 'blue'),
        'chip' : Layer(gds_layer = 99, gds_datatype = 0, description = 'chip area', color = 'gray', alpha = 0.01)
         }

#==============================================================================
# Special functions for this layout
#==============================================================================

def gaas_led(length_LED = 4,
length_taper = 4,
min_taper = 0.2,
width_wg = 0.8,
large_pad_overhang = 0.2,
pad_taper_distance = 200,
small_pad_overhang = 0.2, 
layers = layers,
Pad = Device()):

    LED = Device()
    
    
    pad = LED.add_ref(Pad)
    pad.rotate(90)
    landingpad1 = LED.add_ref(pg.compass(size =pad.size + large_pad_overhang, layer = layers['wg_gaas']))
    landingpad1.center = pad.center
    
    wg1 = LED.add_ref(pg.compass(size = [length_LED, width_wg], layer = layers['wg_gaas']))
    wg2 = LED.add_ref(pg.compass(size = [length_LED, width_wg], layer = layers['wg_algaas']))
    wg1.xmin = landingpad1.xmax + pad_taper_distance
    wg1.y = wg2.y
    wg2.center = wg1.center
    wgtaper1 = LED.add_ref(pr.route_basic(port1 = landingpad1.ports['E'], port2 = wg1.ports['W'], layer = layers['wg_gaas']))
    
    wg2_ramp = LED.add_ref(pg.ramp(length = length_taper, width1 = wg2.ports['E'].width, width2 = min_taper, layer = layers['wg_algaas']))
    wg2_ramp.xmin = wg2.xmax
    wg2_ramp.connect(port = 1, destination = wg2.ports['E'])
    
    wg1_straight = LED.add_ref(pg.compass(size = [length_taper, width_wg], layer = layers['wg_gaas']))
    wg1_straight.connect(port = 'W', destination = wg1.ports['E'])
    
    small_pad = LED.add_ref(pg.compass(size = [length_LED-small_pad_overhang, width_wg-small_pad_overhang], layer = layers['pad']))
    small_pad.center = wg1.center
    small_pad.xmin = wg1.xmin
    small_pad_opening = LED.add_ref(pg.compass(size = [length_LED-2*small_pad_overhang, width_wg-2*small_pad_overhang], layer = layers['padopen_spacer']))
    small_pad_opening.center = small_pad.center
    
    pad_route = LED.add_ref(pr.route_basic(port1 = small_pad.ports['W'], port2 = pad.ports[1], layer = layers['pad']))
    LED.add_port(wg1_straight.ports['E'])
    
    return LED

def gaas_led_pad(width = 100, height = 300, po_offset1 = 20, po_offset2 = 20, pad_layer = 2, po_layer1 = 3, po_layer2 = 4):
    D = Device()
    pad = D.add_ref(pg.compass(size = [width, height], layer = pad_layer))
    pad_opening1 = D.add_ref(pg.compass(size = [width-2*po_offset1, height-2*po_offset1], layer = po_layer1))
    pad_opening2 = D.add_ref(pg.compass(size = [width-2*po_offset2, height-2*po_offset2], layer = po_layer2))

    D.add_port(port=pad.ports['S'], name = 1)
    return D
#==============================================================================
# Define some parameters
#==============================================================================

# Taper constants
LED_parameters = {
            'length_LED' : 4.0,
            'length_taper' : 4.0,
            'min_taper' : 0.2,
            'width_wg' : 0.8,
            'large_pad_overhang' : 0.2,
            'pad_taper_distance' : 200,
            'small_pad_overhang':0.2,
            'layers' : layers,
            }

pad_snspd_parameters = {
            'width' : 100,
            'height' : 300,
            'po_offset' : 0.2,
            'pad_layer' : layers['pad'],
            'po_layer' : layers['padopen_clad']        
            }

pad_led_parameters = {
            'width' : pad_snspd_parameters['width'],
            'height' : pad_snspd_parameters['height'],
            'po_offset1' : pad_snspd_parameters['po_offset'],
            'po_offset2' : 0,
            'pad_layer' : layers['pad'],
            'po_layer1' : layers['padopen_clad'],
            'po_layer2' : layers['padopen_spacer'],   
            }

dbl_pad_snspd_parameters = {
            'gap' : 10,        
            }

SNSPD_parameters = {
        'meander_width' : 0.4,
        'meander_pitch' : 0.8, 
        'num_squares' : 1000, 
        'wg_nw_width' : 0.1, 
        'wg_nw_pitch' : 0.3, 
        'wg_nw_length' : 100,
        'pad_distance' : LED_parameters['pad_taper_distance'], 
        'landing_pad_offset' : 1.0, 
        'nw_layer' : layers['nw'], 
        'wg_layer' : layers['wg_gaas'], 
        'metal_layer' : layers['pad'], 
        }

led_snspd_distance = 10

#==============================================================================
# Make the basic structures
#==============================================================================

#The main device
D = Device()

# make a pad device for the nanowire. We can pass this structure to as many 
# nanowires as we like
Pad_snspd = Device(pg.pad, config = pad_snspd_parameters)

# make a double pad device for the nanowire. 
# We can pass this double pad to as many nanowires as we like)
DBLpad_snspd = Device(pg.dblpad, config = dbl_pad_snspd_parameters, pad_device = Pad_snspd)

# make a pad device for the nanowire. We can pass this structure to as many 
# nanowires as we like
Pad_snspd = Device(gaas_led_pad, config = pad_led_parameters)

#==============================================================================
# Make the more complex structures 
#==============================================================================

# Make the SNSPD. This is called from the basic_photonics script, the one in 
# the geometry script is currently broken.
SNSPD = Device(bp.wg_snspd, config = SNSPD_parameters, dblpad_device = DBLpad_snspd)

# Add the GaAs LED to the main device
led = Device(gaas_led, config = LED_parameters, Pad = Pad_snspd)
D.add_ref(led)

# Add the SNSPD to the main device
wgsnspd = D.add_ref(SNSPD)

# Adjust the position
wgsnspd.xmin = led.xmax + led_snspd_distance
led.y = wgsnspd.ports[1].midpoint[1]

# Make the waveguide that connects them
route = D.add_ref(pr.route_basic(port1 = led.ports['E'], port2 = wgsnspd.ports[1]))

#==============================================================================
# here is where we plot and make device. Could also make arrays here 
#==============================================================================

# plot the structure
quickplot2(D)

# make into a gds file
D.write_gds('GaAs_test.gds')


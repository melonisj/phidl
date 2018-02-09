#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 13:41:42 2017

@author: smb2
"""

from phidl import Device, Layer, quickplot2
import phidl.geometry as pg
import phidl.routing as pr
import gdspy
import numpy as np
import basic_photonics as bp

bs_interaction_length = 50

layers = {
        'full_wg_small' : Layer(gds_layer = 0, gds_datatype = 0, description = 'fully etched Si waveguides small features', color = 'gray'),
        'full_wg_big' : Layer(gds_layer = 0, gds_datatype = 1, description = 'fully etched Si large features', color = 'gold'),
        'shallow_wg_small'  : Layer(gds_layer = 1, gds_datatype = 0, description = 'shallow etched Si waveguide small features', color = 'lightblue', alpha = 0.2),
        'shallow_wg_big' : Layer(gds_layer = 1, gds_datatype = 1, description = 'shallow etched Si waveguide large features', color = 'blue'),
        'nw' : Layer(gds_layer = 2, gds_datatype = 0, description = 'fully etched Si large features', color = 'green'),
        'pad'  : Layer(gds_layer = 3, gds_datatype = 0, description = 'shallow etched Si waveguide small features', color = 'gold', alpha = 0.2),
        'padopen' : Layer(gds_layer = 4, gds_datatype = 0, description = 'shallow etched Si waveguide large features', color = 'red'),
        'ptype' : Layer(gds_layer = 5, gds_datatype = 0, description = 'shallow etched Si waveguide large features', color = 'orange'),
        'ntype' : Layer(gds_layer = 6, gds_datatype = 0, description = 'shallow etched Si waveguide large features', color = 'yellow'),
        'ec' : Layer(gds_layer = 7, gds_datatype = 0, description = 'shallow etched Si waveguide large features', color = 'green'),
        'chip' : Layer(gds_layer = 99, gds_datatype = 0, description = 'chip area', color = 'gray', alpha = 0.01)
         }

D = Device()
pad = pg.pad(width = 100, height = 200, po_offset = 2, pad_layer = layers['pad'], po_layer = layers['padopen'])
dblpad = pg.dblpad(gap = 10, pad_device = pad)


LED = bp.led(width=1, length_wg=10, width_dope_offset=0.2, width_dope=5, wE=1, width_taper = 0.4, length_taper = 10, 
        metal_inset = 0.2, pad_device_distance = [100,0], pad_wire_width = 0.5,
        wg_layer = layers['shallow_wg_small'], p_layer = layers['ptype'], n_layer = layers['ntype'], w_layer = layers['ec'], padtaper_layer = layers['pad'], dblpad_device = dblpad)

wgSNSPD = bp.wg_snspd(meander_width = 0.4, meander_pitch = 0.8, num_squares = 1000, 
            wg_nw_width = 0.1, wg_nw_pitch = 0.3, wg_nw_length = 100, 
            pad_distance = 100, landing_pad_offset = 10.0, nw_layer = layers['nw'], wg_layer = layers['full_wg_small'], metal_layer = layers['pad'], 
            dblpad_device = dblpad)
#led = D.add_ref(LED)


BS = bp.adiabatic_beamsplitter(interaction_length = bs_interaction_length, gap1 = 1, gap2 = 0.1, 
                          port_widths = (0.4, 0.5, 0.4, 0.4),
                          height_sines = LED.size[0]/2, min_radius = 100, 
                          port_devices = (LED,LED,None,None), wg_layer = layers['full_wg_small'])
bs_leds1 = D.add_ref(BS)
bs_leds2 = D.add_ref(BS)
bs_leds3 = D.add_ref(BS)
#bs_leds2 = D.add_ref(BS)
#bs_leds3 = D.add_ref(BS)

BS2 = bp.adiabatic_beamsplitter(interaction_length = bs_interaction_length*4, gap1 = 1, gap2 = 0.1, 
                          port_widths = (0.4, 0.5, 0.4, 0.4),
                          height_sines = LED.size[0]/2, min_radius = 100, 
                          port_devices = (None,None,None,None), wg_layer = layers['full_wg_small'])

ring = BS2.add_ref(pg.ring(radius = 30))
ring2 = BS2.add_ref(pg.ring(radius = 30))
ring3 = BS2.add_ref(pg.ring(radius = 30))
ring.movey(32)
ring.x = BS2.x
ring2.movey(32)
ring3.movey(32)
ring2.x = ring.x-70
ring3.x = ring.x+70

BS3 = bp.adiabatic_beamsplitter(interaction_length = bs_interaction_length, gap1 = 1, gap2 = 0.1, 
                          port_widths = (0.4, 0.5, 0.4, 0.4),
                          height_sines = LED.size[0]/2, min_radius = 100, 
                          port_devices = (None,None,wgSNSPD,wgSNSPD), wg_layer = layers['full_wg_small'])

bs_mid1 = D.add_ref(BS2)
bs_mid2 = D.add_ref(BS2)
bs_mid3 = D.add_ref(BS2)
bs_mid4 = D.add_ref(BS2)

bs_nw1 = D.add_ref(BS3)
bs_nw2 = D.add_ref(BS3)
bs_nw3 = D.add_ref(BS3)

bs_mid1.connect(destination = bs_leds1.ports[3], port = 3)
bs_mid2.connect(destination = bs_leds1.ports[4], port = 4)
bs_leds2.connect(destination = bs_mid2.ports[3], port = 3)
bs_mid3.connect(destination = bs_leds2.ports[4], port = 4)
bs_leds3.connect(destination = bs_mid3.ports[3], port = 3)
bs_mid4.connect(destination = bs_leds3.ports[4], port = 4)
bs_nw1.connect(destination = bs_mid1.ports[1], port = 1)
bs_nw2.connect(destination = bs_mid2.ports[1], port = 1)
bs_nw3.connect(destination = bs_mid3.ports[1], port = 1)

g = D.center
chip = D.add_ref(pg.compass(size = [D.size[0]+500, D.size[1]+100], layer = layers['chip']))
chip.center = g

quickplot2(D)
D.write_gds('concept_device.gds')
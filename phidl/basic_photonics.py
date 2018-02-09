#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 14:33:40 2017

@author: smb2
"""

from phidl import Device, Layer, quickplot
import phidl.geometry as pg
import phidl.routing as pr
import gdspy
import numpy as np

def grating(num_periods = 25, period = 0.6, fill_factor = 0.5, width_grating = 10, length_taper = 100, width = 0.4, minEbeamTaperWidth = 3, partial_etch = False, layerDeepEtchLargeFeatures = 0, layerDeepEtchSmallFeatures = 1, layerPartialEtchLargeFeatures = 2, layerPartialEtchSmallFeatures = 3):
    
    #returns a fiber grating
    G = Device('grating')

# make the deep etched grating
    if partial_etch is False:
        # make the grating teeth
        for i in range(num_periods):
            cgrating = G.add_ref(pg.compass(size=[period*fill_factor,width_grating], layer = layerDeepEtchSmallFeatures))
            cgrating.x+=i*period
            
# make a partially etched grating
    if partial_etch is True:
        # hard coded overlap
            partetch_overhang = 0.5
            # make the etched areas (opposite to teeth)
            for i in range(num_periods):
                cgrating = G.add_ref(pg.compass(size=[period*(1-fill_factor),width_grating+partetch_overhang*2], layer = layerPartialEtchSmallFeatures))
                cgrating.x+=i*period
            
                
        #draw the deep etched square around the grating
            deepbox = G.add_ref(pg.compass(size=[num_periods*period, width_grating], layer=layerDeepEtchLargeFeatures))
            deepbox.center = cgrating.center
            deepbox.xmax = cgrating.xmax + period*(1-fill_factor)
            
            
    # make the taper. First half of taper written on large feature layer, second half on small feature layer
    theta = np.arctan((width_grating-width)/2/length_taper)
    length2 = (minEbeamTaperWidth-width)/2/np.tan(theta)
    tgrating1 = G.add_ref(pg.taper(length = length_taper-length2, width1 = width_grating, width2 = minEbeamTaperWidth, port = None, layer = layerDeepEtchLargeFeatures))
    tgrating2 = G.add_ref(pg.taper(length = length2, width1 = minEbeamTaperWidth, width2 = width, port = None, layer = layerDeepEtchSmallFeatures))
    tgrating1.xmin = cgrating.xmax+(1-fill_factor)*period
    tgrating2.xmin = tgrating1.xmax
    # define the port of the grating
    p = G.add_port(port = tgrating2.ports[2], name = 1)
    
    return G

def wg_fulltaper(width_full_etch = 3, width_shallow_wg = 1, width_full_wg = 0.5, length_shallow_wg = 2, length_taper = 2, layers = {
        'full_wg_small' : Layer(gds_layer = 0, gds_datatype = 0, description = 'fully etched Si waveguides small features', color = 'gray'),
        'full_wg_big' : Layer(gds_layer = 0, gds_datatype = 1, description = 'fully etched Si large features', color = 'gold'),
        'shallow_wg_small'  : Layer(gds_layer = 1, gds_datatype = 0, description = 'shallow etched Si waveguide small features', color = 'lightblue', alpha = 0.2),
        'shallow_wg_big' : Layer(gds_layer = 1, gds_datatype = 1, description = 'shallow etched Si waveguide large features', color = 'blue'),
         }):
    
    # taper from a shallow etch ridge waveguide to a fully etched rib waveguide
    D = Device('wgfulltaper')
    DE = Device('dpetch')
    SE = Device('swetch')

    c_left_fulletch = D.add_ref(pg.connector(midpoint = [0,0], width = width_full_etch, orientation = 0))
    c_left_shallow_wg = D.add_ref(pg.connector(midpoint = [0,0], width = width_shallow_wg, orientation = 0))
    c_right_shallow_wg = D.add_ref(pg.connector(midpoint = [length_shallow_wg,0], width = width_shallow_wg, orientation = 0))
    c_right_fulletch = D.add_ref(pg.connector(midpoint = [length_shallow_wg,0], width = width_full_etch, orientation = 0))
    c_fulletch_wg = D.add_ref(pg.connector(midpoint = [length_shallow_wg+length_taper,0], width = width_full_wg, orientation = 0))
    
    DE.add_ref(pr.route_basic(port1 = c_left_fulletch.ports[1], port2 = c_right_fulletch.ports[2], layer = layers['full_wg_small']))
    SE.add_ref(pr.route_basic(port1 = c_left_shallow_wg.ports[1], port2 = c_right_shallow_wg.ports[2], layer = layers['shallow_wg_small']))
    SE.add_ref(pr.route_basic(port1 = c_right_shallow_wg.ports[2], port2 = c_fulletch_wg.ports[1], layer = layers['shallow_wg_small']))
    DE.add_ref(pr.route_basic(port1 = c_right_fulletch.ports[2], port2 = c_fulletch_wg.ports[1], layer = layers['full_wg_small']))
    
    D.add_ref(subtract(elementA = DE, elementB = SE, layer=layers['shallow_wg_small']))
    deep_etch_region = D.add_ref(DE)
    
    D.add_port(port = c_left_fulletch.ports[2], name = 1)
    D.add_port(port = c_left_shallow_wg.ports[2], name = 2)
    D.add_port(port = c_fulletch_wg.ports[1], name = 3)
        
    return D

def subtract(elementA, elementB, precision = 0.001, layer = 0):
    
    A = Device()
    B = Device()
    if isinstance(elementA, Device): A.add_ref(elementA)
    else: A.elements.append(elementA)
    if isinstance(elementB, Device): B.add_ref(elementB)
    else: B.elements.append(elementB)
        
    gds_layer, gds_datatype = pg._parse_layer(layer)
        
    operandA = A.get_polygons()
    operandB = B.get_polygons()
    p = gdspy.fast_boolean(operandA, operandB, operation = 'not', precision=precision,
                 max_points=199, layer=gds_layer, datatype=gds_datatype)
        
    D = Device()
    D.add_polygon(p, layer=layer)
    return D

def chip_corners(chip_size = 1000.0, width = 5.0, length = 30.0, layer = 0):
    D = Device()
    CR = Device()
    corner1 = pg.compass(size = (width, length), layer = layer)
    corner2 = pg.compass(size = (length-width, width), layer = layer)
    corner2.xmin = corner1.xmax
    corner2.ymax = corner1.ymax
    CR.add_ref(corner1)
    CR.add_ref(corner2)
    
    tl = D.add_ref(CR)
    tl.xmin = -chip_size/2
    tl.ymax = chip_size/2
    tr = D.add_ref(CR)
    tr.rotate(-90)
    tr.xmax = chip_size/2
    tr.ymax = chip_size/2
    bl = D.add_ref(CR)
    bl.rotate(90)
    bl.xmin = -chip_size/2
    bl.ymin = -chip_size/2
    br = D.add_ref(CR)
    br.rotate(180)
    br.xmax = chip_size/2
    br.ymin = -chip_size/2
    
    return D

def wg_snspd(meander_width = 0.4, meander_pitch = 0.8, num_squares = 1000, 
            wg_nw_width = 0.1, wg_nw_pitch = 0.3, wg_nw_length = 100, 
            pad_distance = 500, landing_pad_offset = 10.0, 
            nw_layer = 6, wg_layer = 1, metal_layer = 2, dblpad_device = Device()):
    
    # the length and width of the meander are chosen so that it is approximately 
    # square
    
    D = Device()
    meander_length = np.sqrt(num_squares)*meander_width*2
    num_squares_per_turn = (meander_length/2-(meander_width+meander_pitch))/meander_width + 1
    
    # the number of turns that we will actually make. Must be an odd number if the leads need to come out the same side
    nturns = np.floor(num_squares/num_squares_per_turn)-1
    
    meanderWidth = nturns*meander_pitch
    meanderOffset = meander_pitch/2
    SNSPD = pg.snspd(wire_width = meander_width, wire_pitch = meander_pitch, size = (meander_length,meanderWidth),
              terminals_same_side = True, layer = nw_layer)
    SNSPD.add_port(name = 3, midpoint = [SNSPD.xmin+meander_width/2, SNSPD.ymin], width = meander_width, orientation = -90)
    SNSPD.add_port(name = 4, midpoint = [SNSPD.xmin+meander_width/2, SNSPD.ymax], width = meander_width, orientation = 90)
    meander = D.add_ref(SNSPD)
    wgNw = D.add_ref(pg.optimal_hairpin(width = wg_nw_width, pitch = wg_nw_pitch, length = wg_nw_length, layer = nw_layer))
    wgNw.reflect(p1 = wgNw.ports[1].midpoint, p2 = wgNw.ports[2].midpoint)
    wgNw.xmax = meander.xmin - 2*meanderOffset
    wgNw.ymax = meander.ymin - meanderOffset
    
    meander_size = [meander.bbox[1][0]-meander.bbox[0][1], meander.bbox[0][1]-meander.bbox[1][1]]
    
    # connector between the hairpin and the meander
    D.add_ref(pr.route_basic(port1 = wgNw.ports[1],port2 = SNSPD.ports[2],path_type = 'straight', layer=nw_layer))


    # nw layer pads
    nwPad = D.add_ref(pg.compass(size = [meander_length, meander_width], layer = nw_layer))
    nwPad.ymax = wgNw.ymin - meanderOffset
    nwPad.xmin = wgNw.xmax + 2*meanderOffset 
    D.add_ref(pr.route_basic(port1 = wgNw.ports[2],port2 = nwPad.ports['W'],path_type = 'straight',layer=nw_layer))
     
    # vertical fill rectangles
  
           
    # horizontal fill rectangles
    
    # connectors between nw and pad and meander and pad
    
    
    # metal layer pads
    M1 = pg.compass(size = [SNSPD.size[0], 5], layer = metal_layer)
    metalWire1 = D.add_ref(M1)
    metalWire2= D.add_ref(M1)
    metalWire1.xmin = SNSPD.xmin
    metalWire1.ymin = SNSPD.ymax + meanderOffset
    metalWire2.ymax = nwPad.ymax
    metalWire2.xmin = nwPad.xmin
    NW1 = pg.compass(size = [SNSPD.size[0], 5], layer = nw_layer)
    
    nwBigPad1 = D.add_ref(NW1)
    nwBigPad2 = D.add_ref(NW1)
    nwBigPad1.center = metalWire1.center
    nwBigPad2.center = metalWire2.center
    meander2pad = D.add_ref(pg.compass(size = [meander_width, meanderOffset], layer = nw_layer))
    meander2pad.xmin = SNSPD.xmin
    meander2pad.ymin = SNSPD.ymax
    
    pads = D.add_ref(dblpad_device)
    pads.rotate(angle = -90)
    pads.xmin = meander.xmax + pad_distance
    pads.y = meander.y
    Route1 = pr.route_basic(port1 = pads.ports[1], port2 = metalWire1.ports['E'], path_type = 'straight', layer = metal_layer)
    Route2 = pr.route_basic(port1 = pads.ports[2], port2 = metalWire2.ports['E'], path_type = 'straight', layer = metal_layer)
    D.add_ref(Route1)
    D.add_ref(Route2)
    

    
    # wg layer wg
    
    wg = D.add_ref(pg.compass(size = [wg_nw_length + wg_nw_pitch, wg_nw_pitch*2], layer = wg_layer)) 
    wg.xmax = wgNw.xmax
    wg.y = wgNw.y
    D.add_port(name = 1, port = wg.ports['W'])
    
    # wg layer landing pad
    # SNSPD side
    landingPad1 = D.add_ref(pg.compass(size = [SNSPD.size[0] + 2*meanderOffset + landing_pad_offset, nwBigPad1.ymax - nwBigPad2.ymin +landing_pad_offset], layer = wg_layer))  
    landingPad1.xmin = wg.xmax
    landingPad1.ymin = nwBigPad2.ymin-(landing_pad_offset/2)

    
    #padside
    landingPad2 = D.add_ref(pg.compass(size = pads.size + [landing_pad_offset, landing_pad_offset], layer = wg_layer))
    landingPad2.center = pads.center   
    #routing
    D.add_ref(pr.route_basic(port1 = landingPad1.ports['E'], port2 = landingPad2.ports['W'], layer = wg_layer))
    
   # D.meta['num_squares'] = meander.meta['num_squares']
    return D

def adiabatic_beamsplitter(interaction_length = 10, gap1 = 1, gap2 = 0.1, 
                          port_widths = (0.4, 0.5, 0.4, 0.4),
                          height_sines = 3, min_radius = 10, 
                          port_devices = (None,None,None,None), wg_layer = 0):

    length_sines = np.sqrt((np.pi**2)*height_sines*min_radius/2)
    #
    D = Device("adiabatic beam splitter")
    #
    # start with the actual beamsplitter part
    xpts_upper = [0, 0, interaction_length, interaction_length]
    ypts_upper = [0, port_widths[0], port_widths[2], 0]
    
    wg_upper = D.add_ref(pg.polygon(xpts = xpts_upper, ypts = ypts_upper, layer = wg_layer))
    #
    xpts_lower = [0, 0, interaction_length, interaction_length]
    ypts_lower = [-gap1, -gap1-port_widths[1], -gap2-port_widths[3], -gap2]
    wg_lower = D.add_ref(pg.polygon(xpts = xpts_lower, ypts = ypts_lower, layer = wg_layer))
    
    #locate the straight sections after the sine bends
    P = Device('ports')
    P.add_port(name = 1, midpoint = [wg_upper.xmin-length_sines, wg_upper.center[1]+height_sines], width = port_widths[0], orientation = 0)
    P.add_port(name = 2, midpoint = [wg_lower.xmin-length_sines, wg_lower.center[1]-height_sines], width = port_widths[1], orientation = 0)
    P.add_port(name = 3, midpoint = [wg_upper.xmax+length_sines, wg_upper.center[1]+height_sines], width = port_widths[2], orientation = 180)
    P.add_port(name = 4, midpoint = [wg_lower.xmax+length_sines, wg_lower.center[1]-height_sines], width = port_widths[3], orientation = 180)
    route1 = D.add_ref(pr.route_basic(port1 = P.ports[1], port2 = wg_upper.ports['1'], path_type = 'sine', layer = wg_layer))
    route2 = D.add_ref(pr.route_basic(port1 = P.ports[2], port2 = wg_lower.ports['1'], path_type = 'sine', layer = wg_layer))
    route3 = D.add_ref(pr.route_basic(port1 = P.ports[3], port2 = wg_upper.ports['3'], path_type = 'sine', layer = wg_layer))
    route4 = D.add_ref(pr.route_basic(port1 = P.ports[4], port2 = wg_lower.ports['3'], path_type = 'sine', layer = wg_layer))

    # now we put either devices or ports on the 4 outputs
    dest_ports = [route1.ports[1], route2.ports[1], route3.ports[1], route4.ports[1]]

    for i, PD in enumerate(port_devices):
        if PD is None:
            D.add_port(port = dest_ports[i], name = i+1)
        else:
            pd = D.add_ref(PD)
            pd.connect(port = pd.ports[1], destination = dest_ports[i])
        
    return D

def led(width=1, length_wg=10, width_dope_offset=0.2, width_dope=5, wE=1, width_taper = 0.4, length_taper = 10, 
        metal_inset = 0.2, pad_device_distance = [50,0], pad_wire_width = 0.5,
        wg_layer = 0, p_layer = 1, n_layer = 2, w_layer = 3, padtaper_layer = 4, dblpad_device = Device()):

    D = Device("LED")
    
    wg = D.add_ref(pg.compass(size=[length_wg,width],layer=wg_layer))
    wRegion = D.add_ref(pg.compass(size=[length_wg,wE],layer=w_layer))
    pRegion = D.add_ref(pg.compass(size=[length_wg,width_dope],layer=p_layer))
    nRegion = D.add_ref(pg.compass(size=[length_wg,width_dope],layer=n_layer))
    mytaper = D.add_ref(pg.taper(length = length_taper, width1 = width, width2 = width_taper))
    PW = pg.compass(size = [length_wg, pad_wire_width], layer = padtaper_layer)
    padWire1 = D.add_ref(PW)
    padWire2 = D.add_ref(PW)
    pads = D.add_ref(dblpad_device)
    pads.rotate(angle = 90)
    
    mytaper.xmin = wg.xmax
    wg.connect(port = 'W', destination = mytaper.ports[1])
    wg.center = wRegion.center
    pRegion.ymin = wRegion.ymax + width_dope_offset
    pRegion.center[0] = wRegion.center[0]
    nRegion.ymax = wRegion.ymin - width_dope_offset
    nRegion.center[0] = wRegion.center[0]
    
    padWire1.ymax = pRegion.ymax - metal_inset
    padWire2.ymax = nRegion.ymin + metal_inset
    pads.center = wg.center
    pads.xmax = wg.xmin - pad_device_distance[0]
    pads.movey = pad_device_distance[1]
    D.add_ref(pr.route_basic(port1 = padWire1.ports['W'], port2 =pads.ports[2], layer = padtaper_layer,path_type = 'straight'))
    D.add_ref(pr.route_basic(port1 = padWire2.ports['W'], port2 =pads.ports[1], layer = padtaper_layer,path_type = 'straight'))
    
    D.add_port(port = mytaper.ports[2], name = 1)
    return D

#1x2 y-splitter device.
def y_junction_curved(width_wg=1.1,width_wg_out = 0.35, length_straight=3,y_offset=3,length_split=15,layer=0):
    D=Device()
    B=pg.taper(length=length_straight,width1=width_wg,width2=width_wg_out,layer=layer)
    bin=D.add_ref(pg.taper(length=length_straight,width1=width_wg_out,width2=width_wg,layer=layer))
    bup=D.add_ref(B)

    bup.xmin=bin.xmax+length_split
    bup.y=bin.y+y_offset/2

    bdown=D.add_ref(B)

    bdown.y=bup.y-y_offset
    bdown.x=bup.x

    R=pr.route_basic(port1=bin.ports[2],port2=bup.ports[1],layer=layer)
    r=D.add_ref(R)
    R=pr.route_basic(port1=bin.ports[2],port2=bdown.ports[1],layer=layer)
    r=D.add_ref(R)

    D.add_port(name=1,port=bin.ports[1])
    D.add_port(name=2,port=bup.ports[2])
    D.add_port(name=3,port=bdown.ports[2])

    return D

#paperclip device.  does not currently calculate total length...but that should be easy to implement
def s_cutback2(bend_radius=50,num_sections=2,width_wg=1.1,offset=50,length_straight=8500.0,layer=0):
    D=Device()

    length_reduction=3*bend_radius

    A_cw=pg.arc(radius=bend_radius,width=width_wg,theta=180.0,start_angle=-90,angle_resolution=0.25,layer=layer)
    A_ccw=pg.arc(radius=bend_radius,width=width_wg,theta=-180.0,start_angle=90,angle_resolution=0.25,layer=layer)
    B=pg.taper(length=length_straight,width1=width_wg,width2=width_wg,layer=layer)

    Bin=pg.taper(length=bend_radius*2,width1=width_wg,width2=width_wg,layer=layer)
    bin=D.add_ref(Bin)
    b1=D.add_ref(B)
    b1.connect(port=1,destination=bin.ports[2])
    prev_port=b1.ports[2]
    last_port=b1.ports[2]
    D.add_port(name=1,port=bin.ports[1])
    length_current=length_straight-length_reduction-50
    for x in range(num_sections):
        a = D.add_ref(A_cw)
        a.connect(port=1,destination=prev_port)

        B=pg.taper(length=length_current,width1=width_wg,width2=width_wg,layer=layer)
        b=D.add_ref(B)
        length_current-=length_reduction

        b.connect(port=1,destination=a.ports[2])
        b.y=b.y-2*bend_radius+offset
        b.xmax=a.xmin-length_reduction

        R=pr.route_basic(port1=a.ports[2],port2=b.ports[1],layer=layer)
        r=D.add_ref(R)

        a = D.add_ref(A_ccw)
        a.connect(port=1,destination=b.ports[2])

        if x==num_sections-1:
            B=pg.taper(length=length_current+length_reduction+bend_radius,width1=width_wg,width2=width_wg,layer=layer)
            b=D.add_ref(B)
            b.connect(port=1,destination=a.ports[2])
            prev_port=b.ports[2]
        else:
            B=pg.taper(length=length_current,width1=width_wg,width2=width_wg,layer=layer)
            b=D.add_ref(B)
            length_current-=length_reduction

            b.connect(port=1,destination=a.ports[2])
            b.y=b.y-2*bend_radius+offset
            b.xmin=a.xmax+length_reduction

            R=pr.route_basic(port1=a.ports[2],port2=b.ports[1],layer=layer)
            r=D.add_ref(R)

            prev_port=b.ports[2]

    Bout=pg.taper(length=10,width1=width_wg,width2=width_wg,layer=layer)
    bout=D.add_ref(Bout)
    bout.y=b1.y
    bout.xmin=b1.xmax+bend_radius*2+100
    R=pr.route_manhattan(port1=prev_port,port2=bout.ports[1],layer=layer,radius=bend_radius)
    r=D.add_ref(R)
    D.add_port(name=2,port=bout.ports[2])


    return D


def modeconverter(length = 100,
    w1 = 0.35,
    w2 = 0.8,
    max_shallow_width = 10,
    w3 = 25,
    layer_de = 0,
    layer_se = 1,
    ):
    
    F = Device()
    point1 = [-length/2, w1/2]
    point2 = [-w3/2, max_shallow_width/2]
    point3 = [w3/2, max_shallow_width/2]
    point4 = [length/2, w2/2]
    pointsA = [point1, point2, point3, point4]

    point1 = [-length/2, -w1/2]
    point2 = [-w3/2, -max_shallow_width/2]
    point3 = [w3/2, -max_shallow_width/2]
    point4 = [length/2, -w2/2]
    pointsB = [point1, point2, point3, point4]

    port1 = F.add_port(name = 1, midpoint = [-length/2,0],  width = w1, orientation = 180)
    port2 = F.add_port(name = 2, midpoint = [length/2,0],  width = w2, orientation = 0)

    wg = F.add_ref(pr.route_basic(port1 = port1, port2 = port2, layer = layer_de))
    se1 = F.add_polygon(pointsA, layer_se)
    se2 = F.add_polygon(pointsB, layer_se)
    F.ports[1].orientation = 180
    F.ports[2].orientation = 0
    
    return F
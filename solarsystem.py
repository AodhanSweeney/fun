#! /usr/bin/env python                                                                                                                            
"""This code is dedicated to creating an almost entirely accurate representation of the solar system, showing the sun, and its six
most inner planets. This solar system model has many different scaleings to make sure that the animation runs smothely. All planets 
are approximated as sphere and all orbits as circles. In addition, radii of planets have been adjusted so that they appear on the 
screen but have scaling facors introduced so that some planets are not over powered. Finally, the for the period is also done so that 
objects will rotate more quickly.



Aodhan Sweeney
PHZ 4151C
Feb 13th 2019
Solar System Modeling"""



from __future__ import print_function, division
import numpy as np
import vpython as vp


c1 = .5e-3 # This constant will be multiplied by the radii of our planetary bodies so that we can see what is happening better.
c2 = 50 # Constant multiplied by the angular velocity so that we will see more planetary rotation!

"""The multidimensional array below holds information about the 6 closest orbiting planets ordered from smallest to largest orbital radii. The first column represents 
the radius of the object in km. The next column is the distance from the center of the solar system (or the orbital radius) in millions of meters. The third column is 
the orbital period of the object in days around the sun. Finally the last column is a set of vectors used for the coloring of each planet"""

planetary_data = [[2440,57.9,88, [1,1,1]], 
                  [6052,108.2,224.7, [0,1,1]],
                  [6371,149.6,365.4, [0,0,1]],
                  [3386,227.9,687.0, [1,.2,.2]],
                  [69173,778.5,4331.6, [1,.5,.5]],
                  [57316,1433.4,10759.2, [1,1,1]]]

theta, z = 0.0, 0.0 # Point of initialization for all planets and sun in our system
x = np.cos(theta) # x position and y position must be parameterized with the theta for this polar system
y = np.sin(theta)

sun = vp.sphere(pos=vp.vector(0,0,0), color=vp.vector(1,1,0), radius=45,  make_trail=True, trail_type="points", retain=20) # Our central location for which all planets will orbit about. The sun doesnt need to move.

planet_bodies = []
"""We will be creating a sphere for each of the rows in the planetary data array. Because of this, we will be looping through each row, creating the object, and appending that object to 
a new array which will be the planet_bodies array. To do this, we have to get the x and y position from the orbital radius from the second column in planetary_data, create the sphere with a 
specificed radius, and then color that sphere. This is all done in the first line of the forloop below. """
for planet in planetary_data:
    body = vp.sphere(pos=vp.vector(planet[1]*x, planet[1]*y, 0), color=vp.vector(planet[3][0], planet[3][1], planet[3][2]), radius=(planet[0]*c1), make_trail=True, trail_type="points")
    planet_bodies.append(body)


"""Finally, we need to make sure that the objects move throughout the solar system. Because the objects positions are parameterized by theta, we can write out the coordinates of the objects in
a cartesian way by breaking the x and y components into cos() and sin(). The equations for the position will look something like x = cos(omega*theta) and y = sin(omega*theta). In all cases,
the z component of the position is 0 because we are rotating in the x-y plane. The omega is found in a make shift way and instead of writing in units of radians/sec, we instead do radians/day. Because
the scaling is the same for all bodies it will not make a difference in the relative velocities of the planets. This omega is then = 2*pi/(orbital period in days). We will also multiply in a fudge 
factor of c2 so that the orbits happen faster and the animation isnt boring."""
for theta in np.arange(0, 50*np.pi, 0.1):
    vp.rate(50) # Rate ensures computer will not change screen more than 50 times a second.
    # cartesian position of planet  = (x = radius of orbit*(cos(scaling factor * omega * theta))), (y = radius of orbit*(sin(scaling factor * omega * theta))), z = 0
    planet_bodies[0].pos = vp.vector((planetary_data[0][1]*np.cos(c2*theta*(2*np.pi)/planetary_data[0][2])), (planetary_data[0][1]*np.sin(c2*theta*(2*np.pi)/planetary_data[0][2])),0)
    planet_bodies[1].pos = vp.vector((planetary_data[1][1]*np.cos(c2*theta*(2*np.pi)/planetary_data[1][2])), (planetary_data[1][1]*np.sin(c2*theta*(2*np.pi)/planetary_data[1][2])),0)
    planet_bodies[2].pos = vp.vector((planetary_data[2][1]*np.cos(c2*theta*(2*np.pi)/planetary_data[2][2])), (planetary_data[2][1]*np.sin(c2*theta*(2*np.pi)/planetary_data[2][2])),0)
    planet_bodies[3].pos = vp.vector((planetary_data[3][1]*np.cos(c2*theta*(2*np.pi)/planetary_data[3][2])), (planetary_data[3][1]*np.sin(c2*theta*(2*np.pi)/planetary_data[3][2])),0)
    planet_bodies[4].pos = vp.vector((planetary_data[4][1]*np.cos(c2*theta*(2*np.pi)/planetary_data[4][2])), (planetary_data[4][1]*np.sin(c2*theta*(2*np.pi)/planetary_data[4][2])),0)
    planet_bodies[5].pos = vp.vector((planetary_data[5][1]*np.cos(c2*theta*(2*np.pi)/planetary_data[5][2])), (planetary_data[5][1]*np.sin(c2*theta*(2*np.pi)/planetary_data[5][2])),0)
exit()

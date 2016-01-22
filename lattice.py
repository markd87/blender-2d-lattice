import bpy
import math
import numpy as np
from math import acos
from mathutils import Vector

#define lattice constant
a0=3.18
c=1.51

#define lattice primitive vectors
a1=np.array([a0,0.,0.])
a2=np.array([a0/2.,math.sqrt(3.)/2.*a0,0.])


#nearest neighbors sublattice positions, needed to join with bonds
nb1t=np.array([a0/2.,a0/(2.*math.sqrt(3.)),c])
nb1b=np.array([a0/2.,a0/(2.*math.sqrt(3.)),-c])
nb2t=np.array([-a0/2.,a0/(2.*math.sqrt(3.)),c])
nb2b=np.array([-a0/2.,a0/(2.*math.sqrt(3.)),-c])


#dictionary of atoms positions in the unitcell
atoms={"Mo": {"color": [0.2, 0.2, 1], "radius": 0.6, "location": [0.0, 0.0, 0.0]}, \
       "S1": {"color": [1,1,0.1], "radius": 0.4, "location": [0,-a0/(2.*math.cos(math.pi/6.)),c]}, \
       "S2": {"color": [1,1,0.1], "radius": 0.4, "location": [0,-a0/(2.*math.cos(math.pi/6.)),-c]}}


#draw a single atom
def draw_atom(atom, radius, color, loc):
     x=loc[0]
     y=loc[1]
     z=loc[2]
     bpy.ops.mesh.primitive_uv_sphere_add(size = radius, location = (x,y,z))
     bpy.ops.object.shade_smooth()

     bpy.data.materials.new(name = atom)
     bpy.data.materials[atom].diffuse_color = color
     bpy.data.materials[atom].specular_intensity = 0.2

     #Applies correct material to the new mesh
     bpy.context.active_object.data.materials.append(bpy.data.materials[atom])


#draw a cylinder between neighboring atoms to represent a bond
def draw_bonds(first_loc,second_loc):
	diff = tuple([c2-c1 for c2, c1 in zip(first_loc, second_loc)])
	center = tuple([(c2+c1)/2 for c2, c1 in zip(first_loc, second_loc)])
	magnitude = pow(sum([(c2-c1)**2 for c1, c2 in zip(first_loc, second_loc)]), 0.5)
    
	Vaxis = Vector(diff).normalized()
	Vobj = Vector((0,0,1))
	Vrot = Vobj.cross(Vaxis)
	angle = acos(Vobj.dot(Vaxis))
        
	bpy.ops.mesh.primitive_cylinder_add(radius=0.05, depth=magnitude, location=center)
	bpy.ops.object.shade_smooth()
	bpy.ops.transform.rotate(value=(angle),axis=Vrot)
    
    
#draw lattice with width rows and height cols of unitcells.
def draw_lat(rows,cols):
    for m in np.arange(-rows,cols):
        for n in np.arange(-rows,cols):
            rr=m*a1+n*a2
            for key,atom in atoms.items():
                loc=atom['location']+rr
                draw_atom(key,atom['radius'],atom['color'],loc)
            
            #bonds inside unit cell
            draw_bonds(atoms["Mo"]["location"]+rr,atoms["S1"]["location"]+rr)
            draw_bonds(atoms["Mo"]["location"]+rr,atoms["S2"]["location"]+rr)
            #additional bonds to neighbors
            draw_bonds(atoms["Mo"]["location"]+rr,nb1t+rr)
            draw_bonds(atoms["Mo"]["location"]+rr,nb1b+rr)
            draw_bonds(atoms["Mo"]["location"]+rr,nb2t+rr)
            draw_bonds(atoms["Mo"]["location"]+rr,nb2b+rr)


draw_lat(3,4)


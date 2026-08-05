"""Headless Blender render of an MFC droplet isosurface PLY.
Usage: blender -b -P render_blender.py -- <in.ply> <out.png> [samples] [dir]
Nishita sky environment + water dielectric + auto-framed camera.
dir: view direction preset  '3q' (3/4), 'front', 'side', 'top'
"""

import math
import sys

import bpy
from mathutils import Vector

argv = sys.argv[sys.argv.index("--") + 1 :]
ply_in, png_out = argv[0], argv[1]
samples = int(argv[2]) if len(argv) > 2 else 160
view = argv[3] if len(argv) > 3 else "3q"
style = argv[4] if len(argv) > 4 else "glass"

bpy.ops.wm.read_factory_settings(use_empty=True)
scene = bpy.context.scene

# --- import mesh (keep the smooth gradient normals from marching-cubes) ---
bpy.ops.wm.ply_import(filepath=ply_in, import_colors="NONE")
obj = bpy.context.selected_objects[0]
bpy.ops.object.origin_set(type="ORIGIN_GEOMETRY", center="BOUNDS")
obj.location = (0, 0, 0)
bpy.context.view_layer.objects.active = obj
bpy.ops.object.shade_smooth()  # custom split normals from the PLY drive shading

import os

bb = [obj.matrix_world @ Vector(c) for c in obj.bound_box]
zmin = min(v.z for v in bb)
r = max(v.length for v in bb)  # bounding-sphere radius (centered at origin)
# MFC_FIXED_R: lock camera + ground scale across an animation so the drop grows/
# shrinks in frame instead of the camera zooming to fit each frame.
if os.environ.get("MFC_FIXED_R"):
    r = float(os.environ["MFC_FIXED_R"])
    zmin = -r
size = 2 * r

# --- material (style: glass | frosted | opaque) ---
mat = bpy.data.materials.new("Drop")
mat.use_nodes = True
nt = mat.node_tree
nt.nodes.clear()
out = nt.nodes.new("ShaderNodeOutputMaterial")
bsdf = nt.nodes.new("ShaderNodeBsdfPrincipled")
nt.links.new(bsdf.outputs["BSDF"], out.inputs["Surface"])
if style == "glass":
    bsdf.inputs["Base Color"].default_value = (0.92, 0.96, 1.0, 1.0)
    bsdf.inputs["Roughness"].default_value = 0.0
    bsdf.inputs["IOR"].default_value = 1.333
    bsdf.inputs["Transmission Weight"].default_value = 1.0
    vol = nt.nodes.new("ShaderNodeVolumeAbsorption")
    vol.inputs["Color"].default_value = (0.45, 0.70, 0.92, 1.0)
    vol.inputs["Density"].default_value = 0.12 / max(r, 1e-6)  # faint blue depth
    nt.links.new(vol.outputs["Volume"], out.inputs["Volume"])
elif style == "frosted":
    bsdf.inputs["Base Color"].default_value = (0.80, 0.90, 1.0, 1.0)
    bsdf.inputs["Roughness"].default_value = 0.14
    bsdf.inputs["IOR"].default_value = 1.333
    bsdf.inputs["Transmission Weight"].default_value = 0.9
elif style == "matte":  # diagnostic clay: pure diffuse, no specular
    bsdf.inputs["Base Color"].default_value = (0.55, 0.60, 0.68, 1.0)
    bsdf.inputs["Roughness"].default_value = 1.0
    bsdf.inputs["Specular IOR Level"].default_value = 0.0
else:  # opaque satin resin (water-blue)
    bsdf.inputs["Base Color"].default_value = (0.09, 0.38, 0.82, 1.0)
    bsdf.inputs["Roughness"].default_value = 0.30
obj.data.materials.append(mat)

# --- ground plane (soft studio floor: bounces light up onto the underside) ---
bpy.ops.mesh.primitive_plane_add(size=size * 40, location=(0, 0, zmin - 0.03 * size))
floor = bpy.context.active_object
fm = bpy.data.materials.new("Floor")
fm.use_nodes = True
fb = fm.node_tree.nodes["Principled BSDF"]
fb.inputs["Base Color"].default_value = (0.55, 0.57, 0.60, 1.0)
fb.inputs["Roughness"].default_value = 0.5
floor.data.materials.append(fm)

# --- uniform bright world: even dome light from all directions (clean studio) ---
world = bpy.data.worlds.new("W")
scene.world = world
world.use_nodes = True
wn = world.node_tree
wn.nodes.clear()
wout = wn.nodes.new("ShaderNodeOutputWorld")
bg = wn.nodes.new("ShaderNodeBackground")
bg.inputs["Color"].default_value = (0.62, 0.66, 0.72, 1.0)
bg.inputs["Strength"].default_value = 1.1
wn.links.new(bg.outputs["Background"], wout.inputs["Surface"])


# --- sun key + fill + rim (parallel light: even, no distance/power tuning) ---
def sun(name, rot_deg, energy, color=(1, 1, 1), angle=6):
    ld = bpy.data.lights.new(name, "SUN")
    ld.energy = energy
    ld.color = color
    ld.angle = math.radians(angle)  # soft shadows
    lo = bpy.data.objects.new(name, ld)
    scene.collection.objects.link(lo)
    lo.rotation_euler = tuple(math.radians(a) for a in rot_deg)
    return lo


sun("key", (52, 0, 25), 3.2, (1.0, 0.98, 0.94))  # upper, camera-left
sun("fill", (62, 0, 200), 1.1, (0.92, 0.95, 1.0))  # opposite, dimmer
sun("rim", (-25, 0, 150), 2.0, (1.0, 0.99, 0.96))  # low from behind for edge pop

# --- camera, auto-framed ---
dirs = {"3q": (1.1, -1.6, 0.75), "front": (0, -1, 0.12), "side": (1, 0, 0.12), "top": (0.001, -0.2, 1)}
lens, sensor, rx, ry = 50.0, 36.0, 1600, 1200
half = min(math.atan((sensor / 2) / lens), math.atan((sensor * ry / rx / 2) / lens))
d = 1.25 * r / math.sin(half)
cam_loc = Vector(dirs[view]).normalized() * d
bpy.ops.object.camera_add(location=cam_loc)
cam = bpy.context.active_object
scene.camera = cam
cam.data.lens = lens
bpy.ops.object.empty_add(location=(0, 0, 0))
tgt = bpy.context.active_object
tc = cam.constraints.new("TRACK_TO")
tc.target = tgt
tc.track_axis = "TRACK_NEGATIVE_Z"
tc.up_axis = "UP_Y"

# --- render settings ---
scene.render.engine = "CYCLES"
scene.cycles.device = "CPU"
scene.cycles.samples = samples
scene.cycles.use_denoising = True
scene.cycles.max_bounces = 32
scene.cycles.transmission_bounces = 32
scene.cycles.transparent_max_bounces = 32
scene.cycles.glossy_bounces = 16
scene.render.resolution_x = rx
scene.render.resolution_y = ry
scene.render.image_settings.file_format = "PNG"
scene.render.filepath = png_out
scene.view_settings.view_transform = "Standard"
bpy.ops.render.render(write_still=True)
print("RENDERED", png_out)

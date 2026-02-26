# -*- coding: mbcs -*-
# ============================================================
#  CSV -> Build Cuboid + Individual Fibers (surfaces per fiber)
#  -> BooleanMerge fibers (for cut tool only)
#  -> BooleanCut matrix -> Only_matrix
#  -> Resume original fiber instances (KEEP DEPENDENT for mesh)
#  -> Mesh Only_matrix + EACH fiber part (dependent instances inherit)
#  -> Create matrix vertex/edge/face node sets using CSV bounds
#  -> Apply KUBC BCs (6 load cases) and write 6 input files
# ============================================================

from abaqus import *
from abaqusConstants import *
import part, material, section, assembly
import math
import csv
import regionToolset
import mesh

# ============================================================
# --- SETTINGS ---
# ============================================================
csv_file = 'voronoi_fiber_summary_S1_AF.csv'   # path in your working dir
target_region_id = 59
fiber_radius = 4.0
grow_box_if_needed = True

# Mesh settings
mesh_size = 6.0
deviationFactor = 0.1
minSizeFactor = 0.1

# KUBC strain level
eps0 = 0.01

# ============================================================
# --- HELPERS ---
# ============================================================
def parse_vector_list(vector_string):
    if vector_string is None:
        return []
    s = str(vector_string).strip()
    s = s.replace('(', '').replace(')', '')
    s = s.replace('[', '').replace(']', '')
    chunks = [c.strip() for c in s.split(',') if c.strip()]
    out = []
    for c in chunks:
        parts = c.split()
        if len(parts) >= 3:
            out.append([float(parts[0]), float(parts[1]), float(parts[2])])
    return out

def parse_float_list(s):
    if s is None:
        return []
    return [float(x.strip()) for x in str(s).split(',') if x.strip()]

def bbox_from_fibers(P1_list, P2_list, radius):
    xmin =  1e99; ymin =  1e99; zmin =  1e99
    xmax = -1e99; ymax = -1e99; zmax = -1e99
    for (p1, p2) in zip(P1_list, P2_list):
        for p in (p1, p2):
            xmin = min(xmin, p[0]); ymin = min(ymin, p[1]); zmin = min(zmin, p[2])
            xmax = max(xmax, p[0]); ymax = max(ymax, p[1]); zmax = max(zmax, p[2])
    return [xmin - radius, ymin - radius, zmin - radius], [xmax + radius, ymax + radius, zmax + radius]

def cuboid_bounds(min_corner, size):
    return min_corner, [min_corner[0] + size[0], min_corner[1] + size[1], min_corner[2] + size[2]]

def compute_fiber_volume(lengths, radius):
    A = math.pi * radius * radius
    return A * sum(lengths)

def mesh_part(p, part_name):
    if len(p.cells) == 0:
        raise RuntimeError("Part '{}' has no cells. Boolean op may have produced surface-only geometry.".format(part_name))
    p.seedPart(size=mesh_size, deviationFactor=deviationFactor, minSizeFactor=minSizeFactor)
    p.setMeshControls(regions=p.cells, elemShape=TET, technique=FREE)
    elemType1 = mesh.ElemType(elemCode=C3D20R, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=C3D15,  elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D10,  elemLibrary=STANDARD)
    p.setElementType(regions=(p.cells,), elemTypes=(elemType1, elemType2, elemType3))
    p.generateMesh()
    print("Meshing complete:", part_name)

def set_safe_delete(a1, name):
    if name in a1.sets.keys():
        del a1.sets[name]

def deleteBC_if_exists(model, bcname):
    if bcname in model.boundaryConditions.keys():
        del model.boundaryConditions[bcname]

def apply_sym_bc(model, bcname, region, normal):
    deleteBC_if_exists(model, bcname)
    u1 = UNSET; u2 = UNSET; u3 = UNSET
    if normal.upper() == 'X': u1 = 0.0
    if normal.upper() == 'Y': u2 = 0.0
    if normal.upper() == 'Z': u3 = 0.0
    model.DisplacementBC(
        name=bcname, createStepName='Loading', region=region,
        u1=u1, u2=u2, u3=u3, ur1=UNSET, ur2=UNSET, ur3=UNSET,
        amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
        fieldName='', localCsys=None
    )

def apply_disp_bc(model, bcname, region, dof, value):
    deleteBC_if_exists(model, bcname)
    u1 = UNSET; u2 = UNSET; u3 = UNSET
    if dof == 1: u1 = value
    if dof == 2: u2 = value
    if dof == 3: u3 = value
    model.DisplacementBC(
        name=bcname, createStepName='Loading', region=region,
        u1=u1, u2=u2, u3=u3, ur1=UNSET, ur2=UNSET, ur3=UNSET,
        amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
        fieldName='', localCsys=None
    )

def apply_common_symmetry_except_loaded(model, caseName, x_back, x_front, y_back, y_front, z_back, z_front):
    for nm in ['BC_XBACK_XSYM','BC_XFRONT_XSYM','BC_YBACK_YSYM','BC_YFRONT_YSYM','BC_ZBACK_ZSYM','BC_ZFRONT_ZSYM']:
        deleteBC_if_exists(model, nm)

    if caseName == 'EXX':
        apply_sym_bc(model, 'BC_XBACK_XSYM',  x_back,  'X')
        apply_sym_bc(model, 'BC_YBACK_YSYM',  y_back,  'Y')
        apply_sym_bc(model, 'BC_YFRONT_YSYM', y_front, 'Y')
        apply_sym_bc(model, 'BC_ZBACK_ZSYM',  z_back,  'Z')
        apply_sym_bc(model, 'BC_ZFRONT_ZSYM', z_front, 'Z')

    elif caseName == 'EYY':
        apply_sym_bc(model, 'BC_YBACK_YSYM',  y_back,  'Y')
        apply_sym_bc(model, 'BC_XBACK_XSYM',  x_back,  'X')
        apply_sym_bc(model, 'BC_XFRONT_XSYM', x_front, 'X')
        apply_sym_bc(model, 'BC_ZBACK_ZSYM',  z_back,  'Z')
        apply_sym_bc(model, 'BC_ZFRONT_ZSYM', z_front, 'Z')

    elif caseName == 'EZZ':
        apply_sym_bc(model, 'BC_ZBACK_ZSYM',  z_back,  'Z')
        apply_sym_bc(model, 'BC_XBACK_XSYM',  x_back,  'X')
        apply_sym_bc(model, 'BC_XFRONT_XSYM', x_front, 'X')
        apply_sym_bc(model, 'BC_YBACK_YSYM',  y_back,  'Y')
        apply_sym_bc(model, 'BC_YFRONT_YSYM', y_front, 'Y')

    elif caseName == 'G12':
        apply_sym_bc(model, 'BC_XBACK_XSYM',  x_back,  'X')
        apply_sym_bc(model, 'BC_YBACK_YSYM',  y_back,  'Y')
        apply_sym_bc(model, 'BC_YFRONT_YSYM', y_front, 'Y')
        apply_sym_bc(model, 'BC_ZBACK_ZSYM',  z_back,  'Z')
        apply_sym_bc(model, 'BC_ZFRONT_ZSYM', z_front, 'Z')

    elif caseName == 'G13':
        apply_sym_bc(model, 'BC_XBACK_XSYM',  x_back,  'X')
        apply_sym_bc(model, 'BC_YBACK_YSYM',  y_back,  'Y')
        apply_sym_bc(model, 'BC_YFRONT_YSYM', y_front, 'Y')
        apply_sym_bc(model, 'BC_ZBACK_ZSYM',  z_back,  'Z')
        apply_sym_bc(model, 'BC_ZFRONT_ZSYM', z_front, 'Z')

    elif caseName == 'G23':
        apply_sym_bc(model, 'BC_ZBACK_ZSYM',  z_back,  'Z')
        apply_sym_bc(model, 'BC_ZFRONT_ZSYM', z_front, 'Z')
        apply_sym_bc(model, 'BC_YBACK_YSYM',  y_back,  'Y')
        apply_sym_bc(model, 'BC_XBACK_XSYM',  x_back,  'X')
        apply_sym_bc(model, 'BC_XFRONT_XSYM', x_front, 'X')

# ============================================================
# --- LOAD DATA FROM CSV ---
# ============================================================
region_data = None
with open(csv_file, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        if int(float(row['RegionID'])) == int(target_region_id):
            region_data = row
            break
if region_data is None:
    raise ValueError("Region ID {} not found in CSV!".format(target_region_id))

voronoi_volume = float(region_data.get('VoronoiVolume', '0.0'))
fiber_lengths = parse_float_list(region_data['FiberLengthsInside'])
fiber_P1 = parse_vector_list(region_data['FiberP1_inside'])
fiber_P2 = parse_vector_list(region_data['FiberP2_inside'])
cuboid_min_corner = parse_vector_list(region_data['MinCorner'])[0]
cuboid_size = parse_vector_list(region_data['CuboidSize'])[0]

if len(fiber_lengths) != len(fiber_P1) or len(fiber_lengths) != len(fiber_P2):
    raise ValueError("Mismatch: FiberLengthsInside ({}) vs FiberP1_inside ({}) vs FiberP2_inside ({})"
                     .format(len(fiber_lengths), len(fiber_P1), len(fiber_P2)))

# ============================================================
# --- EXPAND CUBOID TO FIT ALL FIBERS (if needed) ---
# ============================================================
old_min = list(cuboid_min_corner)
old_size = list(cuboid_size)
old_box_vol = old_size[0] * old_size[1] * old_size[2]

fiber_vol = compute_fiber_volume(fiber_lengths, fiber_radius)
old_vf = 100.0 * fiber_vol / old_box_vol if old_box_vol > 0 else 0.0

bbox_min, bbox_max = bbox_from_fibers(fiber_P1, fiber_P2, fiber_radius)
box_min, box_max = cuboid_bounds(old_min, old_size)

new_min = list(old_min)
new_max = list(box_max)

if grow_box_if_needed:
    new_min[0] = min(new_min[0], bbox_min[0])
    new_min[1] = min(new_min[1], bbox_min[1])
    new_min[2] = min(new_min[2], bbox_min[2])
    new_max[0] = max(new_max[0], bbox_max[0])
    new_max[1] = max(new_max[1], bbox_max[1])
    new_max[2] = max(new_max[2], bbox_max[2])

new_size = [new_max[0] - new_min[0], new_max[1] - new_min[1], new_max[2] - new_min[2]]
new_box_vol = new_size[0] * new_size[1] * new_size[2]
new_vf = 100.0 * fiber_vol / new_box_vol if new_box_vol > 0 else 0.0

print("============================================================")
print("RegionID:", target_region_id)
print("Fiber count:", len(fiber_lengths))
print("Fiber radius:", fiber_radius)
print("Fiber volume (sum pi r^2 L):", fiber_vol)
print("Old cuboid minCorner:", old_min)
print("Old cuboid size     :", old_size)
print("Old cuboid volume   :", old_box_vol)
print("Old fiber Vf (%)    :", old_vf)
print("New cuboid minCorner:", new_min)
print("New cuboid size     :", new_size)
print("New cuboid volume   :", new_box_vol)
print("New fiber Vf (%)    :", new_vf)
print("============================================================")

cuboid_min_corner = new_min
cuboid_size = new_size

# ============================================================
# --- CREATE MODEL ---
# ============================================================
model = mdb.models['Model-1']

# ---- Create cuboid part (local 0..L) then translate to minCorner
s1 = model.ConstrainedSketch(name='__cuboid__', sheetSize=1000.0)
s1.rectangle(point1=(0.0, 0.0), point2=(cuboid_size[0], cuboid_size[1]))
cuboid_part_name = 'CuboidRegion' + str(target_region_id)
cuboid_part = model.Part(name=cuboid_part_name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
cuboid_part.BaseSolidExtrude(sketch=s1, depth=cuboid_size[2])

# ============================================================
# --- MATERIALS + SECTIONS ---
# ============================================================
if 'Matrix' not in model.materials.keys():
    model.Material(name='Matrix')
    model.materials['Matrix'].Density(table=((1020.0, ), ))
    model.materials['Matrix'].Conductivity(table=((0.15, ), ))
    model.materials['Matrix'].SpecificHeat(table=((2000.0, ), ))
    model.materials['Matrix'].Expansion(table=((0.00012, ), ))
    model.materials['Matrix'].specificHeat.setValues(law=CONSTANTPRESSURE)
    model.materials['Matrix'].UserMaterial(mechanicalConstants=(
        2244.25, 2242.17, 2241.87, 0.298, 0.2988, 863.63, 863.395, 863.102))
    model.materials['Matrix'].Depvar(n=10)

if 'Fiber' not in model.materials.keys():
    model.Material(name='Fiber')
    model.materials['Fiber'].Elastic(table=((231000, 14600, 14600, 0.25, 0.25, 0.3, 55660, 55660, 5640), ),
                                     type=ENGINEERING_CONSTANTS)
    model.materials['Fiber'].Density(table=((1760.0, ), ))
    model.materials['Fiber'].Expansion(table=((1.5e-05, ), ))
    model.materials['Fiber'].Conductivity(table=((10.64, ), ))
    model.materials['Fiber'].SpecificHeat(table=((795.0, ), ))
    model.materials['Fiber'].specificHeat.setValues(law=CONSTANTPRESSURE)

if 'MatrixSection' not in model.sections.keys():
    model.HomogeneousSolidSection(material='Matrix', name='MatrixSection', thickness=None)
if 'FiberSection' not in model.sections.keys():
    model.HomogeneousSolidSection(material='Fiber', name='FiberSection', thickness=None)

# Step (for later KUBC)
if 'Loading' not in model.steps.keys():
    model.StaticStep(name='Loading', previous='Initial', initialInc=0.01, maxInc=0.1)

# Contact property (kept, even if you don’t use it immediately)
if 'IntProp-1' not in model.interactionProperties.keys():
    model.ContactProperty('IntProp-1')
    model.interactionProperties['IntProp-1'].TangentialBehavior(
        dependencies=0, directionality=ISOTROPIC, elasticSlipStiffness=None,
        formulation=PENALTY, fraction=0.005, maximumElasticSlip=FRACTION,
        pressureDependency=OFF, shearStressLimit=None, slipRateDependency=OFF,
        table=((1.0, ), ), temperatureDependency=OFF)
    model.interactionProperties['IntProp-1'].NormalBehavior(
        allowSeparation=ON, constraintEnforcementMethod=DEFAULT,
        pressureOverclosure=HARD)
    model.interactionProperties['IntProp-1'].HeatGeneration(
        conversionFraction=1.0, secondaryFraction=0.5)

# Assign matrix section
cuboid_part.SectionAssignment(region=(cuboid_part.cells,), sectionName='MatrixSection')

# ============================================================
# --- ASSEMBLY: create cuboid instance and translate ---
# ============================================================
a = model.rootAssembly
a.DatumCsysByDefault(CARTESIAN)
a.Instance(name='CuboidRegionInst', part=cuboid_part, dependent=ON)
a.translate(instanceList=('CuboidRegionInst',), vector=tuple(cuboid_min_corner))

# ============================================================
# --- CREATE + PLACE FIBERS (EACH AS OWN PART + 3 SURFACES) ---
# ============================================================
created_fiber_inst_names = []
created_fiber_part_names = []

for i in range(len(fiber_lengths)):
    length = fiber_lengths[i]
    P1 = fiber_P1[i]
    P2 = fiber_P2[i]

    direction_vector = [P2[j] - P1[j] for j in range(3)]
    magnitude = math.sqrt(direction_vector[0]**2 + direction_vector[1]**2 + direction_vector[2]**2)

    if magnitude < 1e-12 or length <= 0.0:
        print("Warning: Fiber {} has invalid length; skipping.".format(i+1))
        continue

    direction = [direction_vector[0]/magnitude, direction_vector[1]/magnitude, direction_vector[2]/magnitude]

    # Create fiber part aligned along +Z
    sketch = model.ConstrainedSketch(name='FiberSketch_' + str(i+1), sheetSize=1000.0)
    sketch.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(fiber_radius, 0.0))

    part_name = 'Fiber_' + str(target_region_id) + '_' + str(i+1)
    fiber_part = model.Part(name=part_name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    fiber_part.BaseSolidExtrude(sketch=sketch, depth=length)
    fiber_part.SectionAssignment(region=(fiber_part.cells,), sectionName='FiberSection')

    # Optional partitions (your approach)
    try:
        fiber_part.PartitionCellByPlaneThreePoints(
            cells=fiber_part.cells.getSequenceFromMask(('[#1 ]',), ),
            point1=fiber_part.InterestingPoint(fiber_part.edges[0], MIDDLE),
            point2=fiber_part.vertices[0],
            point3=fiber_part.vertices[1]
        )
        if len(fiber_part.edges) >= 7:
            fiber_part.PartitionCellByPlaneThreePoints(
                cells=fiber_part.cells.getSequenceFromMask(('[#3 ]',), ),
                point1=fiber_part.InterestingPoint(fiber_part.edges[6], MIDDLE),
                point2=fiber_part.InterestingPoint(fiber_part.edges[4], MIDDLE),
                point3=fiber_part.InterestingPoint(fiber_part.edges[5], MIDDLE)
            )
        print("  Added partitions to fiber {}".format(i+1))
    except:
        print("  Warning: Could not add partitions to fiber {}".format(i+1))

    # Orientation (same as yours)
    fiber_part.MaterialOrientation(
        region=regionToolset.Region(cells=fiber_part.cells),
        orientationType=SYSTEM,
        axis=AXIS_2,
        additionalRotationType=ROTATION_ANGLE,
        angle=-90.0,
        fieldName='',
        localCsys=None,
        stackDirection=STACK_3,
        additionalRotationField=''
    )

    # ------------------------------------------------------------
    # Create 3 SURFACES on each FIBER PART (local coordinates)
    #   Surf_Lateral : circumference
    #   Surf_End1    : z ~ 0
    #   Surf_End2    : z ~ L
    # ------------------------------------------------------------
    surf_tol = max(1e-6, 1e-6 * length)
    rpad = 1.25 * fiber_radius

    lat_faces = fiber_part.faces.getByBoundingBox(
        xMin=-rpad, yMin=-rpad, zMin=surf_tol,
        xMax= rpad, yMax= rpad, zMax=length - surf_tol
    )
    end1_faces = fiber_part.faces.getByBoundingBox(
        xMin=-rpad, yMin=-rpad, zMin=-surf_tol,
        xMax= rpad, yMax= rpad, zMax= surf_tol
    )
    end2_faces = fiber_part.faces.getByBoundingBox(
        xMin=-rpad, yMin=-rpad, zMin=length - surf_tol,
        xMax= rpad, yMax= rpad, zMax=length + surf_tol
    )

    if len(lat_faces) == 0:
        print("  Warning: lateral surface empty for fiber {}".format(i+1))
    if len(end1_faces) == 0:
        print("  Warning: end1 surface empty for fiber {}".format(i+1))
    if len(end2_faces) == 0:
        print("  Warning: end2 surface empty for fiber {}".format(i+1))

    # Safe overwrite if rerun
    if 'Surf_Lateral' in fiber_part.surfaces.keys(): del fiber_part.surfaces['Surf_Lateral']
    if 'Surf_End1' in fiber_part.surfaces.keys():    del fiber_part.surfaces['Surf_End1']
    if 'Surf_End2' in fiber_part.surfaces.keys():    del fiber_part.surfaces['Surf_End2']

    fiber_part.Surface(name='Surf_Lateral', side1Faces=lat_faces)
    fiber_part.Surface(name='Surf_End1',    side1Faces=end1_faces)
    fiber_part.Surface(name='Surf_End2',    side1Faces=end2_faces)

    # Instance + place (KEEP DEPENDENT=ON)
    inst_name = 'FiberInst_' + str(i+1)
    a.Instance(name=inst_name, part=fiber_part, dependent=ON)
    a.translate(instanceList=(inst_name,), vector=tuple(P1))

    # Rotate fiber to align Z-axis to fiber direction
    z_axis = [0.0, 0.0, 1.0]
    rot_axis = [
        z_axis[1]*direction[2] - z_axis[2]*direction[1],
        z_axis[2]*direction[0] - z_axis[0]*direction[2],
        z_axis[0]*direction[1] - z_axis[1]*direction[0]
    ]
    dot = z_axis[0]*direction[0] + z_axis[1]*direction[1] + z_axis[2]*direction[2]
    dot = min(max(dot, -1.0), 1.0)
    angle_deg = math.degrees(math.acos(dot))

    axis_mag = math.sqrt(rot_axis[0]**2 + rot_axis[1]**2 + rot_axis[2]**2)
    if angle_deg > 1e-6 and axis_mag > 1e-12:
        rot_axis = [rot_axis[0]/axis_mag, rot_axis[1]/axis_mag, rot_axis[2]/axis_mag]
        a.rotate(instanceList=(inst_name,), axisPoint=tuple(P1), axisDirection=tuple(rot_axis), angle=angle_deg)

    created_fiber_inst_names.append(inst_name)
    created_fiber_part_names.append(part_name)

print("Done: Cuboid + individual fibers created for RegionID {}".format(target_region_id))
a.regenerate()

# ============================================================
# --- BOOLEAN: MERGE fibers (tool body) -> CUT matrix -> Only_matrix
# ============================================================
if len(created_fiber_inst_names) == 0:
    raise RuntimeError("No fiber instances were created. Cannot proceed with BooleanCut.")

# 1) Merge all fiber instances into one solid (tool only)
fibers_merged = a.InstanceFromBooleanMerge(
    name='Fibers_tuple',
    instances=[a.instances[nm] for nm in created_fiber_inst_names],
    keepIntersections=ON,
    originalInstances=SUPPRESS,
    domain=GEOMETRY)

# Make merged fibers independent for robust cut
a.makeIndependent(instances=(a.instances['Fibers_tuple-1'],))

# 2) Cut matrix by merged fiber tool -> Only_matrix
matrix_cut = a.InstanceFromBooleanCut(
    name='Only_matrix',
    instanceToBeCut=a.instances['CuboidRegionInst'],
    cuttingInstances=(a.instances['Fibers_tuple-1'],),
    originalInstances=SUPPRESS
)

# Make matrix independent
#a.makeIndependent(instances=(a.instances['Only_matrix-1'],))
a.regenerate()

# 3) Resume original fiber instances (KEEP THEM DEPENDENT; DO NOT makeIndependent)
for nm in created_fiber_inst_names:
    feat = 'FiberInst_{}'.format(nm.split('_')[-1])  # maps FiberInst_12 -> feature FiberInst_12
    if feat in a.features.keys():
        a.features[feat].resume()

a.regenerate()
print("BooleanCut done: matrix -> Only_matrix, fibers resumed (dependent).")

# ============================================================
# --- MESH: Only_matrix + EACH fiber PART (dependent instances inherit)
# ============================================================
MATRIX_PART_NAME = 'Only_matrix'
if MATRIX_PART_NAME not in model.parts.keys():
    raise RuntimeError("Cannot find matrix part '{}' in model.parts".format(MATRIX_PART_NAME))

mesh_part(model.parts[MATRIX_PART_NAME], MATRIX_PART_NAME)

# Mesh each fiber part
for part_name in created_fiber_part_names:
    if part_name in model.parts.keys():
        mesh_part(model.parts[part_name], part_name)
    else:
        print("Warning: fiber part missing for meshing:", part_name)

a.regenerate()
print("Meshing complete: Only_matrix + all individual fibers.")

# ============================================================
# --- CREATE MATRIX SETS (vertices/edges/faces) using CSV bounds
#     IMPORTANT: matrix instance is Only_matrix-1
# ============================================================
a1 = model.rootAssembly
MATRIX_INST_NAME = 'Only_matrix-1'
if MATRIX_INST_NAME not in a1.instances.keys():
    raise RuntimeError("Cannot find matrix instance '{}' in assembly.".format(MATRIX_INST_NAME))
matrixInstance = a1.instances[MATRIX_INST_NAME]

# Bounds from known translation + size (works even without getBoundingBox)
x0, y0, z0 = cuboid_min_corner[0], cuboid_min_corner[1], cuboid_min_corner[2]
Lx, Ly, Lz = cuboid_size[0], cuboid_size[1], cuboid_size[2]

xmin, ymin, zmin = x0, y0, z0
xmax, ymax, zmax = x0 + Lx, y0 + Ly, z0 + Lz

wid, hei, thk = Lx, Ly, Lz
tolX = 1e-3 * wid if wid != 0 else 1e-6
tolY = 1e-3 * hei if hei != 0 else 1e-6
tolZ = 1e-3 * thk if thk != 0 else 1e-6

print("Matrix bounds (from CSV):")
print("xmin,xmax =", xmin, xmax)
print("ymin,ymax =", ymin, ymax)
print("zmin,zmax =", zmin, zmax)

# ---- Vertex sets (geom + nodes)
vertexloc = [
    (xmin, ymin, zmin),
    (xmax, ymin, zmin),
    (xmax, ymax, zmin),
    (xmin, ymax, zmin),
    (xmin, ymin, zmax),
    (xmax, ymin, zmax),
    (xmax, ymax, zmax),
    (xmin, ymax, zmax),
]
for i, pt in enumerate(vertexloc, start=1):
    setname = 'vertex_{}'.format(i)
    set_safe_delete(a1, setname)
    set_safe_delete(a1, setname + '_nodes')
    v = matrixInstance.vertices.findAt((pt,))
    a1.Set(vertices=v, name=setname)
    a1.Set(nodes=a1.sets[setname].nodes, name=setname + '_nodes')

# ---- Face sets (geom + nodes)
for nm in ['X_inf_face','X_sup_face','Y_inf_face','Y_sup_face','Z_inf_face','Z_sup_face',
           'X_inf_face_nodes','X_sup_face_nodes','Y_inf_face_nodes','Y_sup_face_nodes','Z_inf_face_nodes','Z_sup_face_nodes']:
    set_safe_delete(a1, nm)

faces = matrixInstance.faces
XinfF = faces.getByBoundingBox(xmin - tolX, ymin - 10*tolY, zmin - 10*tolZ,
                               xmin + tolX, ymax + 10*tolY, zmax + 10*tolZ)
XsupF = faces.getByBoundingBox(xmax - tolX, ymin - 10*tolY, zmin - 10*tolZ,
                               xmax + tolX, ymax + 10*tolY, zmax + 10*tolZ)
YinfF = faces.getByBoundingBox(xmin - 10*tolX, ymin - tolY, zmin - 10*tolZ,
                               xmax + 10*tolX, ymin + tolY, zmax + 10*tolZ)
YsupF = faces.getByBoundingBox(xmin - 10*tolX, ymax - tolY, zmin - 10*tolZ,
                               xmax + 10*tolX, ymax + tolY, zmax + 10*tolZ)
ZinfF = faces.getByBoundingBox(xmin - 10*tolX, ymin - 10*tolY, zmin - tolZ,
                               xmax + 10*tolX, ymax + 10*tolY, zmin + tolZ)
ZsupF = faces.getByBoundingBox(xmin - 10*tolX, ymin - 10*tolY, zmax - tolZ,
                               xmax + 10*tolX, ymax + 10*tolY, zmax + tolZ)

a1.Set(faces=XinfF, name='X_inf_face')
a1.Set(nodes=a1.sets['X_inf_face'].nodes, name='X_inf_face_nodes')
a1.Set(faces=XsupF, name='X_sup_face')
a1.Set(nodes=a1.sets['X_sup_face'].nodes, name='X_sup_face_nodes')

a1.Set(faces=YinfF, name='Y_inf_face')
a1.Set(nodes=a1.sets['Y_inf_face'].nodes, name='Y_inf_face_nodes')
a1.Set(faces=YsupF, name='Y_sup_face')
a1.Set(nodes=a1.sets['Y_sup_face'].nodes, name='Y_sup_face_nodes')

a1.Set(faces=ZinfF, name='Z_inf_face')
a1.Set(nodes=a1.sets['Z_inf_face'].nodes, name='Z_inf_face_nodes')
a1.Set(faces=ZsupF, name='Z_sup_face')
a1.Set(nodes=a1.sets['Z_sup_face'].nodes, name='Z_sup_face_nodes')

# ---- Edge sets (optional; can be messy after BooleanCut, but included)
for nm in ['edge_{}'.format(i) for i in range(1,13)] + ['edge_{}_nodes'.format(i) for i in range(1,13)]:
    set_safe_delete(a1, nm)

edges = matrixInstance.edges

def make_edge_set(i, xMin,xMax,yMin,yMax,zMin,zMax):
    ed = edges.getByBoundingBox(xMin,xMax,yMin,yMax,zMin,zMax)
    if len(ed) == 0:
        print("Warning: edge_{} empty".format(i))
        return
    a1.Set(edges=ed, name='edge_{}'.format(i))
    a1.Set(nodes=a1.sets['edge_{}'.format(i)].nodes, name='edge_{}_nodes'.format(i))

make_edge_set(1, xmin - tolX, xmax + tolX, ymin - tolY, ymin + tolY, zmin - tolZ, zmin + tolZ)
make_edge_set(3, xmin - tolX, xmax + tolX, ymax - tolY, ymax + tolY, zmin - tolZ, zmin + tolZ)
make_edge_set(2, xmax - tolX, xmax + tolX, ymin - tolY, ymax + tolY, zmin - tolZ, zmin + tolZ)
make_edge_set(4, xmin - tolX, xmin + tolX, ymin - tolY, ymax + tolY, zmin - tolZ, zmin + tolZ)

make_edge_set(5, xmin - tolX, xmax + tolX, ymin - tolY, ymin + tolY, zmax - tolZ, zmax + tolZ)
make_edge_set(7, xmin - tolX, xmax + tolX, ymax - tolY, ymax + tolY, zmax - tolZ, zmax + tolZ)
make_edge_set(6, xmax - tolX, xmax + tolX, ymin - tolY, ymax + tolY, zmax - tolZ, zmax + tolZ)
make_edge_set(8, xmin - tolX, xmin + tolX, ymin - tolY, ymax + tolY, zmax - tolZ, zmax + tolZ)

make_edge_set(9,  xmin - tolX, xmin + tolX, ymin - tolY, ymin + tolY, zmin - tolZ, zmax + tolZ)
make_edge_set(10, xmax - tolX, xmax + tolX, ymin - tolY, ymin + tolY, zmin - tolZ, zmax + tolZ)
make_edge_set(11, xmax - tolX, xmax + tolX, ymax - tolY, ymax + tolY, zmin - tolZ, zmax + tolZ)
make_edge_set(12, xmin - tolX, xmin + tolX, ymax - tolY, ymax + tolY, zmin - tolZ, zmax + tolZ)

# ============================================================
# --- KUBC BCs + 6 input files
# ============================================================
# Use face node sets on Only_matrix-1
x_back  = a1.sets['X_inf_face_nodes']
x_front = a1.sets['X_sup_face_nodes']
y_back  = a1.sets['Y_inf_face_nodes']
y_front = a1.sets['Y_sup_face_nodes']
z_back  = a1.sets['Z_inf_face_nodes']
z_front = a1.sets['Z_sup_face_nodes']

# Output requests
if 'F-Output-1' in model.fieldOutputRequests.keys():
    model.fieldOutputRequests['F-Output-1'].setValues(variables=('S','U','IVOL','E'))

L = wid  # characteristic length

cases = [
    ('EXX', eps0),
    ('EYY', eps0),
    ('EZZ', eps0),
    ('G12', eps0),
    ('G13', eps0),
    ('G23', eps0),
]

# for (caseName, eps) in cases:
    # apply_common_symmetry_except_loaded(model, caseName, x_back, x_front, y_back, y_front, z_back, z_front)

    # deleteBC_if_exists(model, 'BC_LOAD')
    # disp = eps * L

    # if caseName == 'EXX':
        # apply_disp_bc(model, 'BC_LOAD', x_front, 1, disp)
    # elif caseName == 'EYY':
        # apply_disp_bc(model, 'BC_LOAD', y_front, 2, disp)
    # elif caseName == 'EZZ':
        # apply_disp_bc(model, 'BC_LOAD', z_front, 3, disp)
    # elif caseName == 'G12':
        # apply_disp_bc(model, 'BC_LOAD', y_front, 1, disp)
    # elif caseName == 'G13':
        # apply_disp_bc(model, 'BC_LOAD', z_front, 1, disp)
    # elif caseName == 'G23':
        # apply_disp_bc(model, 'BC_LOAD', z_front, 2, disp)

    # jobname = 'Mesoscale_' + caseName
    # if jobname in mdb.jobs.keys():
        # del mdb.jobs[jobname]

    # mdb.Job(
        # name=jobname, model='Model-1', description='KUBC case: ' + caseName,
        # type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0,
        # queue=None, memory=50, memoryUnits=PERCENTAGE,
        # getMemoryFromAnalysis=False, explicitPrecision=SINGLE,
        # nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF,
        # contactPrint=OFF, historyPrint=OFF, userSubroutine='',
        # scratch='', multiprocessingMode=DEFAULT, numCpus=5, numDomains=5
    # )
    # mdb.jobs[jobname].writeInput(consistencyChecking=OFF)

# print("DONE: Individual fibers meshed (dependent), 3 surfaces per fiber part created, matrix cut + meshed, sets + 6 KUBC input files written.")
# print("Fiber surfaces usage example:")
# print("  a.instances['FiberInst_1'].surfaces['Surf_Lateral']")
# print("  a.instances['FiberInst_1'].surfaces['Surf_End1']")
# print("  a.instances['FiberInst_1'].surfaces['Surf_End2']")
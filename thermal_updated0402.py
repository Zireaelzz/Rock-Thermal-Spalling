# -*- coding: mbcs -*-
# ==============================================================================
# Script: Step 1 - Thermal Prediction (Updated mesh only; thermal physics unchanged)
# Logic:
# 1. Keep original thermal loading / cavity radiation / output unchanged.
# 2. Only replace fragile findAt-based geometry edge selection and mesh seeding.
# ==============================================================================

from abaqus import *
from abaqusConstants import *
import abaqusConstants
import mesh, part, material, section, assembly, step, interaction, load, regionToolset, job
import os, math

# ==============================================================================
# Helper functions: ONLY for robust edge detection / mesh seeding / surface creation
# ==============================================================================

def vertex_xy(vtx):
    p = vtx.pointOn[0]
    return p[0], p[1]

def edge_vertices_xy(part_or_inst, edge):
    vids = edge.getVertices()
    v1 = part_or_inst.vertices[vids[0]]
    v2 = part_or_inst.vertices[vids[1]]
    x1, y1 = vertex_xy(v1)
    x2, y2 = vertex_xy(v2)
    return x1, y1, x2, y2

def unique_edges(edge_list):
    seen = set()
    out = []
    for e in edge_list:
        if e.index not in seen:
            out.append(e)
            seen.add(e.index)
    return out

def as_edge_array(part_or_inst, edge_list):
    if len(edge_list) == 0:
        return part_or_inst.edges[0:0]
    arr = part_or_inst.edges[edge_list[0].index:edge_list[0].index+1]
    for e in edge_list[1:]:
        arr = arr + part_or_inst.edges[e.index:e.index+1]
    return arr

def classify_rock_edges(p_rock, split_x, rock_x_max, tol=1e-6):
    """
    Multi-segment compatible rock edge classifier.
    Returns edge LISTS:
      bot_left_edges, bot_right_edges, arc_left_edges, arc_right_edges,
      part_edges, right_edges
    """
    out = {
        'bot_left_edges':  [],
        'bot_right_edges': [],
        'arc_left_edges':  [],
        'arc_right_edges': [],
        'part_edges':      [],
        'right_edges':     [],
    }

    for e in p_rock.edges:
        x1, y1, x2, y2 = edge_vertices_xy(p_rock, e)
        xm, ym = 0.5*(x1+x2), 0.5*(y1+y2)

        # bottom edges
        if abs(y1) < tol and abs(y2) < tol:
            if xm < split_x:
                out['bot_left_edges'].append(e)
            else:
                out['bot_right_edges'].append(e)
            continue

        # partition line
        if abs(x1 - split_x) < tol and abs(x2 - split_x) < tol:
            out['part_edges'].append(e)
            continue

        # right boundary
        if abs(x1 - rock_x_max) < tol and abs(x2 - rock_x_max) < tol:
            out['right_edges'].append(e)
            continue

        # remaining upper edges = arc candidates
        if ym > 0.0:
            if xm < split_x:
                out['arc_left_edges'].append(e)
            else:
                out['arc_right_edges'].append(e)

    for k in out.keys():
        out[k] = unique_edges(out[k])
        print('[classify_rock_edges] {}: {}'.format(k, len(out[k])))

    return out

def classify_tree_edges(p_tree, tree_w, tree_h=2.0, tol=1e-6):
    """
    Tree edge classifier for robust Surf-Tree creation.
    """
    out = {
        'bottom_edges': [],
        'left_edges':   [],
        'right_edges':  [],
        'top_edges':    [],
    }

    for e in p_tree.edges:
        x1, y1, x2, y2 = edge_vertices_xy(p_tree, e)

        if abs(y1) < tol and abs(y2) < tol:
            out['bottom_edges'].append(e)
            continue

        if abs(x1) < tol and abs(x2) < tol:
            out['left_edges'].append(e)
            continue

        if abs(x1 - tree_w) < tol and abs(x2 - tree_w) < tol:
            out['right_edges'].append(e)
            continue

        if abs(y1 - tree_h) < tol and abs(y2 - tree_h) < tol:
            out['top_edges'].append(e)
            continue

    for k in out.keys():
        out[k] = unique_edges(out[k])
        print('[classify_tree_edges] {}: {}'.format(k, len(out[k])))

    return out


# --- 1. Setup ---
work_dir = r'C:\ABAQUS2025\temp\researchproject\researchyue\project1'
csv_file = os.path.join(work_dir, 'amplitude.csv')

if not os.path.exists(csv_file):
    if not os.path.exists(work_dir):
        try:
            os.makedirs(work_dir)
        except:
            pass
    try:
        with open(csv_file, 'w') as f:
            f.write('0.0, 20.0\n60.0, 800.0\n86400.0, 800.0\n')
    except:
        pass

model_name = 'Rock_Tree_Cavity_FinalRun'
job_name = 'Job_Cavity_FinalRun_24h'

if model_name in mdb.models:
    del mdb.models[model_name]
if job_name in mdb.jobs:
    del mdb.jobs[job_name]

myModel = mdb.Model(name=model_name)
myModel.setValues(absoluteZero=-273.15, stefanBoltzmann=5.67e-8)

# ==============================================================================
# 2. Geometry & Mesh & Surface (Part Level)
# ==============================================================================
tree_w, gap, rock_r = 0.4, 3.0, 4.0
tree_x_min, tree_x_max = 0.0, tree_w
rock_x_min, rock_x_max = tree_w + gap, tree_w + gap + rock_r
split_x = (rock_x_min + rock_x_max) / 2.0

# --- Part 1: Rock ---
s_rock = myModel.ConstrainedSketch(name='s_rock', sheetSize=20.0)
s_rock.Line(point1=(rock_x_min, 0.0), point2=(rock_x_max, 0.0))
s_rock.Line(point1=(rock_x_max, 0.0), point2=(rock_x_max, rock_r))
s_rock.ArcByCenterEnds(center=(rock_x_max, 0.0),
                       point1=(rock_x_max, rock_r),
                       point2=(rock_x_min, 0.0))
p_rock = myModel.Part(name='Rock', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
p_rock.BaseShell(sketch=s_rock)

# Partition line: only remove fragile findAt, keep partition itself unchanged
t = p_rock.MakeSketchTransform(
    sketchPlane=p_rock.faces[0],
    sketchUpEdge=p_rock.edges[0],
    sketchPlaneSide=SIDE1,
    origin=(0, 0, 0)
)
s_part = myModel.ConstrainedSketch(name='partition', sheetSize=20, transform=t)
s_part.Line(point1=(split_x, -1.0), point2=(split_x, rock_r + 1.0))
p_rock.PartitionFaceBySketch(faces=p_rock.faces[:], sketch=s_part)

# Create Rock Surface & Set
# Keep original meaning: Surf-Granite = the whole exposed rock arc
rock_edges = classify_rock_edges(
    p_rock=p_rock,
    split_x=split_x,
    rock_x_max=rock_x_max,
    tol=1e-6
)
arc_left_edges  = rock_edges['arc_left_edges']
arc_right_edges = rock_edges['arc_right_edges']
surf_granite_arr = as_edge_array(p_rock, arc_left_edges + arc_right_edges)

p_rock.Surface(side1Edges=surf_granite_arr, name='Surf-Granite')
p_rock.Set(faces=p_rock.faces, name='All')

# --- Part 2: Tree ---
s_tree = myModel.ConstrainedSketch(name='s_tree', sheetSize=20.0)
s_tree.rectangle(point1=(tree_x_min, 0.0), point2=(tree_x_max, 2.0))
p_tree = myModel.Part(name='Tree', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
p_tree.BaseShell(sketch=s_tree)

# Create Tree Surface & Set
# Keep original meaning: Surf-Tree = the fire-facing right edge of tree
tree_edges = classify_tree_edges(p_tree, tree_w, tree_h=2.0, tol=1e-6)
surf_tree_arr = as_edge_array(p_tree, tree_edges['right_edges'])

p_tree.Surface(side1Edges=surf_tree_arr, name='Surf-Tree')
p_tree.Set(faces=p_tree.faces, name='All')  # 保持原样

# --- Materials & Section ---
m_rock = myModel.Material(name='Granite_User')
m_rock.Density(table=((2650.,),))
m_rock.Conductivity(table=((2.5,),))
m_rock.SpecificHeat(table=((800.,),))
m_rock.Elastic(table=((67e9, 0.25),))
myModel.HomogeneousSolidSection(name='Sec-Rock', material='Granite_User', thickness=1.0)
p_rock.SectionAssignment(region=p_rock.sets['All'], sectionName='Sec-Rock')

m_tree = myModel.Material(name='Tree_Mat')
m_tree.Density(table=((800.,),))
m_tree.Conductivity(table=((0.1,),))
m_tree.SpecificHeat(table=((1500.,),))
myModel.HomogeneousSolidSection(name='Sec-Tree', material='Tree_Mat', thickness=1.0)
p_tree.SectionAssignment(region=p_tree.sets['All'], sectionName='Sec-Tree')

# ==============================================================================
# 纯热模型：只改网格抓边方式，保留原 mesh 参数
# ==============================================================================
elem_type = mesh.ElemType(elemCode=DC2D4)
p_rock.setMeshControls(regions=p_rock.faces[:], elemShape=QUAD, technique=STRUCTURED)
p_rock.setElementType(regions=p_rock.sets['All'], elemTypes=(elem_type,))

p_tree.setMeshControls(regions=p_tree.faces[:], elemShape=QUAD, technique=STRUCTURED)
p_tree.setElementType(regions=p_tree.sets['All'], elemTypes=(elem_type,))

print(">>> Applying precise local edge seeds for Thermal Model...")

bot_left_edges  = rock_edges['bot_left_edges']
bot_right_edges = rock_edges['bot_right_edges']
part_edges      = rock_edges['part_edges']
right_edges     = rock_edges['right_edges']

all_ok = (
    len(bot_left_edges)  > 0 and
    len(bot_right_edges) > 0 and
    len(arc_left_edges)  > 0 and
    len(arc_right_edges) > 0 and
    len(part_edges)      > 0 and
    len(right_edges)     > 0
)

if all_ok:
    bot_left_arr  = as_edge_array(p_rock, bot_left_edges)
    bot_right_arr = as_edge_array(p_rock, bot_right_edges)
    arc_left_arr  = as_edge_array(p_rock, arc_left_edges)
    arc_right_arr = as_edge_array(p_rock, arc_right_edges)
    part_arr      = as_edge_array(p_rock, part_edges)
    right_arr     = as_edge_array(p_rock, right_edges)

    # 保持你原来的局部布种参数不变
    p_rock.seedEdgeByBias(end1Edges=bot_left_arr, minSize=0.01, maxSize=0.03)
    p_rock.seedEdgeBySize(edges=arc_left_arr, size=0.03)
    p_rock.seedEdgeBySize(edges=part_arr, size=0.05)
    p_rock.seedEdgeByBias(end1Edges=arc_right_arr, minSize=0.03, maxSize=0.1)
    p_rock.seedEdgeBySize(edges=right_arr, size=0.1)
    p_rock.seedEdgeByBias(end1Edges=bot_right_arr, minSize=0.03, maxSize=0.1)
else:
    print("WARNING: rock edge classification incomplete; fallback to seedPart(size=0.08)")
    p_rock.seedPart(size=0.08)

# 生成网格
p_rock.generateMesh()

# 树模型保持均匀网格 (原样)
p_tree.seedPart(size=0.1)
p_tree.generateMesh()

# ==============================================================================
# 3. Assembly & Step   <-- 保持原样
# ==============================================================================
myAssembly = myModel.rootAssembly
inst_rock = myAssembly.Instance(name='Rock-1', part=p_rock, dependent=ON)
inst_tree = myAssembly.Instance(name='Tree-1', part=p_tree, dependent=ON)

surf_rock_ref = inst_rock.surfaces['Surf-Granite']
surf_tree_ref = inst_tree.surfaces['Surf-Tree']

step_name = 'Step-Thermal-Pred'
myModel.HeatTransferStep(
    name=step_name, previous='Initial', response=TRANSIENT,
    timePeriod=86400.0, maxNumInc=20000,
    initialInc=0.1, minInc=1e-5, maxInc=60.0, deltmx=50.0
)
myModel.FieldOutputRequest(name='F-Output-1', createStepName=step_name, variables=('NT',))

# ==============================================================================
# 4. Amplitudes & Boundary Conditions   <-- 保持原样
# ==============================================================================
def read_amplitude_csv(filename):
    data = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split(',')
                if len(parts) >= 2:
                    data.append((float(parts[0]), float(parts[1])))
    except:
        pass
    if not data:
        return ((0.0, 20.0), (60.0, 800.0), (86400.0, 800.0))
    return tuple(data)

amp_data = read_amplitude_csv(csv_file)
myModel.TabularAmplitude(name='Amp_Fire_CSV', data=amp_data, smooth=SOLVER_DEFAULT)

# 4.1 初始环境温度：全局 20.0 度
myModel.Temperature(
    name='Init_Temp',
    createStepName='Initial',
    region=myAssembly.Set(faces=inst_rock.faces + inst_tree.faces, name='All_Nodes'),
    distributionType=UNIFORM,
    magnitudes=(20.0,)
)

# 4.2 树的热源边界条件 (原样)
print(">>> Applying Fire Temperature BC to the Tree...")
myModel.TemperatureBC(
    name='BC-Tree-Fire',
    createStepName=step_name,
    region=inst_tree.sets['All'],
    distributionType=UNIFORM,
    magnitude=1.0,
    amplitude='Amp_Fire_CSV',
    fixed=OFF
)

# ==============================================================================
# 5. Interaction   <-- 保持原样
# ==============================================================================
print(">>> Creating Cavity Radiation Properties...")
myModel.CavityRadiationProp(name='Prop-Tree', property=((0.98,),))
myModel.CavityRadiationProp(name='Prop-Granite', property=((0.90,),))

print(">>> Configuring Cavity Radiation Interaction...")
myModel.CavityRadiation(
    name='Int-Cavity',
    createStepName=step_name,
    surfaces=(surf_tree_ref, surf_rock_ref),
    surfaceEmissivities=('Prop-Tree', 'Prop-Granite'),
    ambientTemp=20.0,
    viewfactorAccurTol=0.05,
    blocking=abaqusConstants.BLOCKING_ALL
)

# ==============================================================================
# 6. Submit & View   <-- 保持原样
# ==============================================================================
print(">>> Creating Job: {}...".format(job_name))
myJob = mdb.Job(name=job_name, model=model_name, numCpus=4, numDomains=4)

print(">>> Submitting Job...")
myJob.submit()
myJob.waitForCompletion()

try:
    session.viewports[session.currentViewportName].setValues(displayedObject=myAssembly)
    session.viewports[session.currentViewportName].assemblyDisplay.setValues(interactions=ON)
except:
    pass

print(">>> DONE! ODB generated: {}.odb".format(job_name))
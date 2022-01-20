"""
    Compute 2D LV map using conformal flattening and 2 spatial constraints: apex fixed in the center + aorta point in 
    the border (angle=pi in the disk)

    Launch GUI and ask the user to select 3 seeds:
    1. LV apex
    2. point in the base contour (or close to) that should be placed in pi in the disk
    3. point in the base contour (or close to) next to seed 2 and in anticlockwise when the RV is
    located on the left side of the LV. This point is used to find the correct contour direction.

    After the quasi-conformal with constraints flattening, get more uniform spatial distribution of points using the
    radial displacement presented in: "Patient independent representation of the detailed cardiac ventricular anatomy."
    Bruno Paun, et al. Medical image analysis (2017)

    Input: LV mesh
    Optional: RV mesh and n for radial displacement
    Output: flat LV mesh
    Example usage: python flat_LV.py --lv_meshfile data/lv_mesh.vtk --rv_meshfile data/rv_mesh.vtk 
"""

from aux_functions import *
from seedselector import *
import math, argparse, os, sys
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument('--lv_meshfile', type=str, metavar='PATH', help='path to LV input mesh')
parser.add_argument('--rv_meshfile', type=str, default='', metavar='PATH', help='path to RV input mesh')
parser.add_argument('--n', type=float, default=1/3.0, help='n for radial displacement')
args = parser.parse_args()

if os.path.isfile(args.lv_meshfile):
    surface_lv = readvtk(args.lv_meshfile)
else:
    sys.exit('ERROR: LV input file does not exist. Please, specify a valid path.')

fileroot = os.path.dirname(args.lv_meshfile)
filename = os.path.basename(args.lv_meshfile)
filenameroot = os.path.splitext(filename)[0]

# if RV mesh available, append polydatas
if os.path.isfile(args.rv_meshfile):
    surface_rv = readvtk(args.rv_meshfile)
    surface = append(surface_lv, surface_rv)
else:
    surface = surface_lv

# define disk parameters (external radio, apex and aorta point position)
rdisk = 1.0
ap_x0 = np.array([0.0])
ap_y0 = np.array([0.0])

seedsfile = os.path.join(fileroot, filenameroot + '_seeds.vtk')
nseeds = 3
labels = [0, 1, 2]
if not os.path.exists(seedsfile):
    print('Select 3 seeds in this order: \n 1. LV apex \n 2. Aorta point \n 3. Auxiliary point')
    seeds = seed_interactor(surface)
    if not seeds.GetNumberOfIds() == nseeds:
        print('You should select exactly', nseeds, ' seeds. Try again!')
        seeds = seed_interactor(surface)
    newpoints = vtk.vtkPoints()
    newvertices = vtk.vtkCellArray()
    labels_array = vtk.vtkDoubleArray()
    labels_array.SetName('seed_label')
    for s in range(seeds.GetNumberOfIds()):
        label = labels[s]
        point = surface.GetPoint(seeds.GetId(s))
        pid = newpoints.InsertNextPoint(point)
        labels_array.InsertNextValue(label)
        # Create the topology of the point (vertex)
        newvertices.InsertNextCell(1)
        newvertices.InsertCellPoint(pid)
    pointspd = vtk.vtkPolyData()
    pointspd.SetPoints(newpoints)
    pointspd.SetVerts(newvertices)
    pointspd.GetPointData().AddArray(labels_array)
    writevtk(pointspd, seedsfile)
seeds_poly = readvtk(seedsfile)

# Find corresponding seeds in the surface mesh
locator = vtk.vtkPointLocator()
locator.SetDataSet(surface_lv)
locator.BuildLocator()
id_ap = int(locator.FindClosestPoint(seeds_poly.GetPoint(0)))

# Detect edges -> base contour
# cont = extractboundaryedge(surface_lv)
cont = extractlargestregion(extractboundaryedge(m_no_base))   # Ensure it's only 1
edge_cont_ids = get_ordered_cont_ids_based_on_distance(cont).astype(int)

# find corresponding ordered points in the COMPLETE mesh. Use same locator as before
cont_base_ids = np.zeros(edge_cont_ids.shape[0]) - 1
for i in range(cont_base_ids.shape[0]):
    p = cont.GetPoint(edge_cont_ids[i])
    cont_base_ids[i] = locator.FindClosestPoint(p)

# find closest point in base contour to 'Aorta' seed
locator_base = vtk.vtkPointLocator()
locator_base.SetDataSet(cont)
locator_base.BuildLocator()
id_0 = locator_base.FindClosestPoint(seeds_poly.GetPoint(1))
id_aux = locator_base.FindClosestPoint(seeds_poly.GetPoint(2))
# in surface
id_base0 = locator.FindClosestPoint(cont.GetPoint(id_0))
id_base_aux = locator.FindClosestPoint(cont.GetPoint(id_aux))

# base -> external disk
# order cont_base_ids to start in 'aorta'
ref_base = int(locator.FindClosestPoint(surface_lv.GetPoint(id_base0)))
reordered_base_cont = np.append(cont_base_ids[int(np.where(cont_base_ids == ref_base)[0]): cont_base_ids.shape[0]],
                              cont_base_ids[0: int(np.where(cont_base_ids == ref_base)[0])]).astype(int)
# print('Reference point in base is ', ref_base)

# check if the list of ordered points corresponding to the base contours has to be flipped
pos_auxpoint = int(np.where(reordered_base_cont == id_base_aux)[0])
if pos_auxpoint > 20:     # 20 points far... empirical, it will depend on the mesh elements :(
    # Flip
    print('I ll flip the base ids')
    aux = np.zeros(reordered_base_cont.size)
    for i in range(reordered_base_cont.size):
        aux[reordered_base_cont.size - 1 - i] = reordered_base_cont[i]
    reordered_base_cont = np.append(aux[aux.size - 1], aux[0:aux.size - 1]).astype(int)

complete_circumf_t = np.linspace(np.pi, np.pi + 2*np.pi, len(reordered_base_cont), endpoint=False)  # starting in pi
x0_ext = np.cos(complete_circumf_t) * rdisk
y0_ext = np.sin(complete_circumf_t) * rdisk

m_disk = flat_w_constraints(surface_lv, reordered_base_cont, np.array([id_ap]), x0_ext, y0_ext, ap_x0.astype(float),
                            ap_y0.astype(float))

transfer_all_scalar_arrays_by_point_id(surface_lv, m_disk)
# writevtk(m_disk, os.path.join(fileroot, filenameroot + '_flat.vtk'))

# Apply Bruno's radial displacement to enlarge central part and get more uniform mesh
# Paun, Bruno, et al. "Patient independent representation of the detailed cardiac ventricular anatomy."
# Medical image analysis (2017)
npoints = m_disk.GetNumberOfPoints()
m_out = vtk.vtkPolyData()
points = vtk.vtkPoints()

points.SetNumberOfPoints(npoints)
n = args.n
for i in range(npoints):
    q = np.array(m_disk.GetPoint(i))
    p = np.copy(q)
    x = math.pow(math.pow(p[0], 2) + math.pow(p[1], 2), np.divide(n, 2.0))*math.cos(math.atan2(p[1], p[0]))
    y = math.pow(math.pow(p[0], 2) + math.pow(p[1], 2), np.divide(n, 2.0))*math.sin(math.atan2(p[1], p[0]))
    points.SetPoint(i, x, y, 0)
m_out.SetPoints(points)
m_out.SetPolys(m_disk.GetPolys())

transfer_all_scalar_arrays_by_point_id(m_disk, m_out)
writevtk(m_out, os.path.join(fileroot, filenameroot + '_flat_radial_displacement' + str(1/n) + ' .vtk'))
# save with 1/n name to avoid having 0.333333 or similar in the filename

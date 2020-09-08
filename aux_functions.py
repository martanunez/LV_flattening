import vtk
import math
import numpy as np
from scipy import sparse
import scipy.sparse.linalg as linalg_sp
from scipy.sparse import vstack, hstack, coo_matrix, csc_matrix
import seedselector


def seed_interactor(surface):
    """Interactor for seed selection. Needs VMTK"""
    computer = seedselector.vmtkPickPointSeedSelector()
    computer.SetSurface(surface)
    computer.Execute()
    return computer.GetSourceSeedIds()


###     Input/Output    ###
def readvtk(filename):
    """Read VTK file"""
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()


def writevtk(surface, filename, type='ascii'):
    """Write binary or ascii VTK file"""
    writer = vtk.vtkPolyDataWriter()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        writer.SetInputData(surface)
    else:
        writer.SetInput(surface)
    writer.SetFileName(filename)
    if type == 'ascii':
        writer.SetFileTypeToASCII()
    elif type == 'binary':
        writer.SetFileTypeToBinary()
    writer.Write()


def append(polydata1, polydata2):
    """Define new polydata appending polydata1 and polydata2"""
    appender = vtk.vtkAppendPolyData()
    appender.AddInputData(polydata1)
    appender.AddInputData(polydata2)
    appender.Update()
    return appender.GetOutput()


def extractboundaryedge(polydata):
    edge = vtk.vtkFeatureEdges()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        edge.SetInputData(polydata)
    else:
        edge.SetInput(polydata)
    edge.FeatureEdgesOff()
    edge.NonManifoldEdgesOff()
    edge.Update()
    return edge.GetOutput()


def find_create_path(mesh, p1, p2):
    """Get shortest path (using Dijkstra algorithm) between p1 and p2 on the mesh. Returns a polydata"""
    dijkstra = vtk.vtkDijkstraGraphGeodesicPath()
    if vtk.vtkVersion().GetVTKMajorVersion() > 5:
        dijkstra.SetInputData(mesh)
    else:
        dijkstra.SetInput(mesh)
    dijkstra.SetStartVertex(p1)
    dijkstra.SetEndVertex(p2)
    dijkstra.Update()
    return dijkstra.GetOutput()


def compute_geodesic_distance(mesh, id_p1, id_p2):
    """Compute geodesic distance from point id_p1 to id_p2 on surface 'mesh'
    It first computes the path across the edges and then the corresponding distance adding up point to point distances)"""
    path = find_create_path(mesh, id_p1, id_p2)
    total_dist = 0
    n = path.GetNumberOfPoints()
    for i in range(n-1):   # Ids are ordered in the new polydata, from 0 to npoints_in_path
        p0 = path.GetPoint(i)
        p1 = path.GetPoint(i+1)
        dist = math.sqrt(math.pow(p0[0]-p1[0], 2) + math.pow(p0[1]-p1[1], 2) + math.pow(p0[2]-p1[2], 2) )
        total_dist = total_dist + dist
    return total_dist, path


def transfer_all_scalar_arrays_by_point_id(m1, m2):
    """ Transfer all scalar arrays from m1 to m2 by point id"""
    for i in range(m1.GetPointData().GetNumberOfArrays()):
        print('Transferring scalar array: {}'.format(m1.GetPointData().GetArray(i).GetName()))
        m2.GetPointData().AddArray(m1.GetPointData().GetArray(i))


def get_ordered_cont_ids_based_on_distance(mesh):
    """ Given a contour, get the ordered list of Ids (not ordered by default).
    Open the mesh duplicating the point with id = 0. Compute distance transform of point 0
    and get a ordered list of points (starting in 0) """
    m = vtk.vtkMath()
    m.RandomSeed(0)
    # copy the original mesh point by point
    points = vtk.vtkPoints()
    polys = vtk.vtkCellArray()
    cover = vtk.vtkPolyData()
    nver = mesh.GetNumberOfPoints()
    points.SetNumberOfPoints(nver+1)

    new_pid = nver  # id of the duplicated point
    added = False

    for j in range(mesh.GetNumberOfCells()):
        # get the 2 point ids
        ptids = mesh.GetCell(j).GetPointIds()
        cell = mesh.GetCell(j)
        if (ptids.GetNumberOfIds() != 2):
            # print "Non contour mesh (lines)"
            break

        # read the 2 involved points
        pid0 = ptids.GetId(0)
        pid1 = ptids.GetId(1)
        p0 = mesh.GetPoint(ptids.GetId(0))   # returns coordinates
        p1 = mesh.GetPoint(ptids.GetId(1))

        if pid0 == 0:
            if added == False:
                # Duplicate point 0. Add gaussian noise to the original point
                new_p = [p0[0] + m.Gaussian(0.0, 0.0005), p0[1] + m.Gaussian(0.0, 0.0005), p0[2] + m.Gaussian(0.0, 0.0005)]
                points.SetPoint(new_pid, new_p)
                points.SetPoint(pid1, p1)
                polys.InsertNextCell(2)
                polys.InsertCellPoint(pid1)
                polys.InsertCellPoint(new_pid)
                added = True
            else:  # act normal
                points.SetPoint(ptids.GetId(0), p0)
                points.SetPoint(ptids.GetId(1), p1)
                polys.InsertNextCell(2)
                polys.InsertCellPoint(cell.GetPointId(0))
                polys.InsertCellPoint(cell.GetPointId(1))
        elif pid1 == 0:
            if added == False:
                new_p = [p1[0] + m.Gaussian(0.0, 0.0005), p1[1] + m.Gaussian(0.0, 0.0005), p1[2] + m.Gaussian(0.0, 0.0005)]
                points.SetPoint(new_pid, new_p)
                points.SetPoint(pid0, p0)
                polys.InsertNextCell(2)
                polys.InsertCellPoint(pid0)
                polys.InsertCellPoint(new_pid)
                added = True
            else:  # act normal
                points.SetPoint(ptids.GetId(0), p0)
                points.SetPoint(ptids.GetId(1), p1)
                polys.InsertNextCell(2)
                polys.InsertCellPoint(cell.GetPointId(0))
                polys.InsertCellPoint(cell.GetPointId(1))

        else:
            points.SetPoint(ptids.GetId(0), p0)
            points.SetPoint(ptids.GetId(1), p1)
            polys.InsertNextCell(2)
            polys.InsertCellPoint(cell.GetPointId(0))
            polys.InsertCellPoint(cell.GetPointId(1))

    if added == False:
        print('Warning: I have not added any point, list of indexes may not be correct.')
    cover.SetPoints(points)
    cover.SetPolys(polys)
    if not vtk.vtkVersion.GetVTKMajorVersion() > 5:
        cover.Update()
    # compute distance from point with id 0 to all the rest
    npoints = cover.GetNumberOfPoints()
    dists = np.zeros(npoints)
    for i in range(npoints):
        [dists[i], poly] = compute_geodesic_distance(cover, int(0), i)
    list_ = np.argsort(dists).astype(int)
    return list_[0:len(list_)-1]    # skip last one, duplicated


def ExtractVTKPoints(mesh):
    """Extract points from vtk structures. Return the Nx3 numpy.array of the vertices."""
    n = mesh.GetNumberOfPoints()
    vertex = np.zeros((n, 3))
    for i in range(n):
        mesh.GetPoint(i, vertex[i, :])
    return vertex


def ExtractVTKTriFaces(mesh):
    """Extract triangular faces from vtkPolyData. Return the Nx3 numpy.array of the faces (make sure there are only triangles)."""
    m = mesh.GetNumberOfCells()
    faces = np.zeros((m, 3), dtype=int)
    for i in range(m):
        ptIDs = vtk.vtkIdList()
        mesh.GetCellPoints(i, ptIDs)
        if ptIDs.GetNumberOfIds() != 3:
            raise Exception("Nontriangular cell!")
        faces[i, 0] = ptIDs.GetId(0)
        faces[i, 1] = ptIDs.GetId(1)
        faces[i, 2] = ptIDs.GetId(2)
    return faces


def ComputeLaplacian(vertex, faces):
    """Calculates the laplacian of a mesh
    vertex 3xN numpy.array: vertices
    faces 3xM numpy.array: faces"""
    n = vertex.shape[1]
    m = faces.shape[1]

    # compute mesh weight matrix
    W = sparse.coo_matrix((n, n))
    for i in np.arange(1, 4, 1):
        i1 = np.mod(i - 1, 3)
        i2 = np.mod(i, 3)
        i3 = np.mod(i + 1, 3)
        pp = vertex[:, faces[i2, :]] - vertex[:, faces[i1, :]]
        qq = vertex[:, faces[i3, :]] - vertex[:, faces[i1, :]]
        # normalize the vectors
        pp = pp / np.sqrt(np.sum(pp ** 2, axis=0))
        qq = qq / np.sqrt(np.sum(qq ** 2, axis=0))

        # compute angles
        ang = np.arccos(np.sum(pp * qq, axis=0))
        W = W + sparse.coo_matrix((1 / np.tan(ang), (faces[i2, :], faces[i3, :])), shape=(n, n))
        W = W + sparse.coo_matrix((1 / np.tan(ang), (faces[i3, :], faces[i2, :])), shape=(n, n))

    # compute laplacian
    d = W.sum(axis=0)
    D = sparse.dia_matrix((d, 0), shape=(n, n))
    L = D - W
    return L


def flat_w_constraints(m, boundary_ids, constraints_ids, x0_b, y0_b, x0_c, y0_c):
    """ Conformal flattening fitting boundary points to (x0_b,y0_b) coordinate positions
    and additional contraint points to (x0_c,y0_c).
    Solve minimization problem using quadratic programming: https://en.wikipedia.org/wiki/Quadratic_programming"""

    penalization = 1000
    vertex = ExtractVTKPoints(m).T    # 3 x n_vertices
    faces = ExtractVTKTriFaces(m).T
    n = vertex.shape[1]
    L = ComputeLaplacian(vertex, faces)
    L = L.tolil()
    L[boundary_ids, :] = 0.0     # Not conformal there
    for i in range(boundary_ids.shape[0]):
         L[boundary_ids[i], boundary_ids[i]] = 1

    L = L*penalization

    Rx = np.zeros(n)
    Ry = np.zeros(n)
    Rx[boundary_ids] = x0_b * penalization
    Ry[boundary_ids] = y0_b * penalization

    L = L.tocsr()
    # result = np.zeros((Rx.size, 2))

    nconstraints = constraints_ids.shape[0]
    M = np.zeros([nconstraints, n])   # M, zero rows except 1 in constraint point
    for i in range(nconstraints):
        M[i, constraints_ids[i]] = 1
    dx = x0_c
    dy = y0_c

    block1 = hstack([L.T.dot(L), M.T])

    zeros_m = coo_matrix(np.zeros([len(dx),len(dx)]))
    block2 = hstack([M, zeros_m])

    C = vstack([block1, block2])

    prodx = coo_matrix([L.T.dot(Rx)])
    dxx = coo_matrix([dx])
    cx = hstack([prodx, dxx])

    prody = coo_matrix([L.T.dot(Ry)])
    dyy = coo_matrix([dy])
    cy = hstack([prody, dyy])

    solx = linalg_sp.spsolve(C, cx.T)
    soly = linalg_sp.spsolve(C, cy.T)

    # print('There are: ', len(np.argwhere(np.isnan(solx))), ' nans')
    # print('There are: ', len(np.argwhere(np.isnan(soly))), ' nans')
    if len(np.argwhere(np.isnan(solx))) > 0:
        print('WARNING!!! matrix is singular. It is probably due to the convergence of 2 different division lines in the same point.')
        print('Trying to assign different 2D possition to same 3D point. Try to create new division lines or increase resolution of mesh.')

    pd = vtk.vtkPolyData()
    pts = vtk.vtkPoints()

    pts.SetNumberOfPoints(n)
    for i in range(n):
        pts.SetPoint(i, solx[i], soly[i], 0)

    pd.SetPoints(pts)
    pd.SetPolys(m.GetPolys())
    pd.Modified()
    return pd
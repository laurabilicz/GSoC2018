import vtk
import numpy as np
import random

"""
Pipeline to convert skeleton data from a nifti file to vtkPolyData.
Works only small array yet because python doesn't allow deep recursion.

"""

def def_globals():
    reader = vtk.vtkNIFTIImageReader()
    reader.SetFileName("data.nii.gz")
    reader.Update()

    array = reader.GetOutput().GetPointData().GetScalars()

    global a
    a = vtk_to_numpy(array)
    a.shape = reader.GetOutput().GetDimensions()

    global look_up
    look_up = np.zeros(a.shape)
    look_up = look_up - 1

    global points
    points = vtk.vtkPoints()


def in_bounds(x, y, z):
    if x >= 0 and y >= 0 and z >= 0 and x < a.shape[0] and y < a.shape[1] and z < a.shape[2]:
        return True
    return False


def neighbor_finder(i, j, k, it):
    id = -1
    if a[i, j, k] > 0 and look_up[i, j, k] == -1:
        print("new point: ", i, j, k)
        id = points.InsertNextPoint([i, j, k])
        look_up[i,j,k] = id
        if a[i, j, k] != 2 or it == 0:
            it += 1
            for x in range(-1, 2):
                for y in range(-1, 2):
                     for z in range(-1, 2):
                        if not(x == 0 and y == 0 and z == 0) and in_bounds(i+x, j+y, k+z) and look_up[i+x, j+y, k+z] == -1 and a[i+x, j+y, k+z] > 0:
                            print(i+x, j+y, k+z)
                            id = neighbor_finder(i+x, j+y, k+z, it)

    return id




def_globals()

# The indexes array contains the end point of each line so I can reconstruct the lines.
indexes = []

for l in range(a.shape[0]):
    for m in range(a.shape[1]):
        for n in range(a.shape[2]):
            if a[l, m, n] == 2:
                id = neighbor_finder(l, m, n, 0)
                if id != -1:
                    indexes.append(id)
                    print(indexes)


colors = vtk.vtkNamedColors()

cells = vtk.vtkCellArray()

polylines = [vtk.vtkPolyLine() for i in range(len(indexes))]

for i in range(0, len(indexes)):
    polylines[i] = vtk.vtkPolyLine()
    if i == 0:
        polylines[i].GetPointIds().SetNumberOfIds(indexes[i] + 1)
        for j in range(0, indexes[i]+1):
            #print(j, j)
            polylines[i].GetPointIds().SetId(j, j)
    else:
        print(indexes[i] - indexes[i - 1])
        polylines[i].GetPointIds().SetNumberOfIds(indexes[i] - indexes[i - 1])
        for j in range(indexes[i-1]+1, indexes[i]+1):
            #print(j-indexes[i-1]-1, j)
            polylines[i].GetPointIds().SetId(j-indexes[i-1]-1, j)
    cells.InsertNextCell(polylines[i])
    polyline = None

# RENDERING STUFF

polydata = vtk.vtkPolyData()
polydata.SetPoints(points)
polydata.SetLines(cells)

mapper = vtk.vtkPolyDataMapper()
mapper.SetInputData(polydata)


# ACTOR
actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetColor(colors.GetColor3d("Tomato"))
actor.GetProperty().SetLineWidth(4)


#TUBE
# create radius
tubeRadius = vtk.vtkDoubleArray()
n = polydata.GetNumberOfPoints()
tubeRadius.SetNumberOfTuples(n)
tubeRadius.SetName("TubeRadius")

for i in range(n):
    tubeRadius.SetTuple1(i, random.random())

# Add the scalars to the polydata
tubePolyData = polydata
tubePolyData.GetPointData().AddArray(tubeRadius)
tubePolyData.GetPointData().SetActiveScalars("TubeRadius")

tubeFilter = vtk.vtkTubeFilter()
tubeFilter.SetInputData(tubePolyData)
tubeFilter.SetNumberOfSides(50)
tubeFilter.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
tubeFilter.Update()

tubeMapper = vtk.vtkPolyDataMapper()
tubeMapper.SetInputConnection(tubeFilter.GetOutputPort())

# Create a mapper
tubeMapper = vtk.vtkPolyDataMapper()
tubeMapper.SetInputConnection(tubeFilter.GetOutputPort())

tubeActor = vtk.vtkActor()
tubeActor.SetMapper(tubeMapper)

# RENDERING
ren = vtk.vtkRenderer()

renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)

# INTERACTOR
iren = vtk.vtkRenderWindowInteractor()
iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
iren.SetRenderWindow(renWin)

ren.AddActor(actor)
ren.AddActor(tubeActor)
ren.SetBackground(0.1, 0.2, 0.4)
renWin.SetSize(400, 400)

iren.Initialize()

# We'll zoom in a little by accessing the camera and invoking a "Zoom"
# method on it.
ren.ResetCamera()
ren.GetActiveCamera().Zoom(1.5)
renWin.Render()

# Start the event loop.
iren.Start()

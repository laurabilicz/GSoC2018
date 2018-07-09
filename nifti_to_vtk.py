import vtk
import numpy as np
from scipy.interpolate import interp1d
import nibabel as nib


# ================================================== UTIL ==============================================================


def read_vtk(file):
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(file)
    reader.Update()

    point_mapper = vtk.vtkPolyDataMapper()
    point_mapper.SetInputData(reader.GetOutput())

    point_actor = vtk.vtkActor()
    point_actor.SetMapper(point_mapper)
    point_actor.GetProperty().SetPointSize(5)

    return point_actor


def transform_polydata(actor):
    mapper = actor.GetMapper()
    poly = mapper.GetInput()
    return poly


def write_vtk(actor, filename):
    polydata = transform_polydata(actor)
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(filename + ".vtk")
    writer.SetInputData(polydata)
    writer.Write()

# ==========================  SKELETON  ====================================


def in_bounds(a, x, y, z):
    if x >= 0 and y >= 0 and z >= 0 and x < a.shape[0] and y < a.shape[1] and z < a.shape[2]:
        return True
    return False


def new_coord(a, look_up, i, j, k):
    result = []
    for x in range(-1, 2):
        for y in range(-1, 2):
            for z in range(-1, 2):
                if (not (x == 0 and y == 0 and z == 0)) and in_bounds(a, i + x, j + y, k + z):
                    if look_up[i + x, j + y, k + z] == -1 and a[i + x, j + y, k + z] > 0:
                        result.append(i+x)
                        result.append(j+y)
                        result.append(k+z)
                        return result

    return result


def neighbor_finder(a, look_up, points_array, i, j, k):
    cont = True
    points = []
    while cont:
        cont = False
        if a[i, j, k] > 0:  # and look_up[i, j, k] == -1:
            # first iteration
            if look_up[i, j, k] == -1:
                print("new point: ", i, j, k)
                points.append([i, j, k])
                result = new_coord(a, look_up, i, j, k)
                look_up[i, j, k] = 1
                if len(result) == 3:
                    i = result[0]
                    j = result[1]
                    k = result[2]
                    cont = True
            else:
                points.append([i, j, k])

    points_array.append(points)

    return [i, j, k]


def end_finder(point, points_array, skeleton_points):
    i = point[0]
    j = point[1]
    k = point[2]
    found = False
    end_point = []
    for x in range(-1, 2):
        for y in range(-1, 2):
            for z in range(-1, 2):
                if (not (x == 0 and y == 0 and z == 0)) and in_bounds(points_array, i + x, j + y, k + z) \
                 and points_array[i + x, j + y, k + z] > 0:
                    end_point = [i + x, j + y, k + z]
                    found = True
    if found:
        skeleton_points[-1].append(end_point)


def make_skeleton(points_array, diameter_array):
    look_up = np.zeros(points_array.shape)
    look_up = look_up - 1

    skeleton_points = []

    for i in range(points_array.shape[0]):
        for j in range(points_array.shape[1]):
            for k in range(points_array.shape[2]):
                """
                for i in range(100,150):
                    for j in range(100, 150):
                        for k in range(100, 150):
                """
                if points_array[i, j, k] == 2:
                    print("___________")
                    end_point = neighbor_finder(points_array, look_up, skeleton_points, i, j, k)
                    if points_array[i, j, k] == 1 or end_point == [i, j, k]:
                        end_finder(end_point, points_array, skeleton_points)
    print("======================")

    tube_actor, line_actor = make_tubes(skeleton_points, diameter_array)

    return tube_actor, line_actor

# ============================== RENDERING ===========================================


def make_polydata(points, r):
    # fit a spline to points
    x_spline = vtk.vtkKochanekSpline()
    y_spline = vtk.vtkKochanekSpline()
    z_spline = vtk.vtkKochanekSpline()

    spline = vtk.vtkParametricSpline()
    spline.SetXSpline(x_spline)
    spline.SetYSpline(y_spline)
    spline.SetZSpline(z_spline)
    spline.SetPoints(points)

    function_source = vtk.vtkParametricFunctionSource()
    function_source.SetParametricFunction(spline)
    function_source.Update()

    # create radius
    tube_radius = vtk.vtkDoubleArray()
    n = function_source.GetOutput().GetNumberOfPoints()
    tube_radius.SetNumberOfTuples(n)
    tube_radius.SetName("TubeRadius")

    radius = interp1d(np.linspace(0, 1, num=len(r), endpoint=True), r)

    for i in range(n):
        tube_radius.SetTuple1(i, radius(i / n))

    # Add the scalars to the polydata
    tube_poly_data = function_source.GetOutput()
    tube_poly_data.GetPointData().AddArray(tube_radius)
    tube_poly_data.GetPointData().SetActiveScalars("TubeRadius")

    return tube_poly_data


def make_actors(tube_poly_data):
    tuber = vtk.vtkTubeFilter()
    tuber.SetInputData(tube_poly_data)
    tuber.SetNumberOfSides(20)
    tuber.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
    tuber.Update()

    # Setup actors and mappers
    line_mapper = vtk.vtkPolyDataMapper()
    line_mapper.SetInputData(tube_poly_data)
    line_mapper.SetScalarRange(tube_poly_data.GetScalarRange())

    tube_mapper = vtk.vtkPolyDataMapper()
    tube_mapper.SetInputConnection(tuber.GetOutputPort())
    tube_mapper.SetScalarRange(tube_poly_data.GetScalarRange())

    line_actor = vtk.vtkActor()
    line_actor.SetMapper(line_mapper)

    tube_actor = vtk.vtkActor()
    tube_actor.SetMapper(tube_mapper)

    return tube_actor, line_actor


def assembly(polydata_array):
    append_filter = vtk.vtkAppendPolyData()
    for polydata in polydata_array:
        append_filter.AddInputData(polydata)
    append_filter.Update()

    clean_filter = vtk.vtkCleanPolyData()
    clean_filter.SetInputConnection(append_filter.GetOutputPort())
    clean_filter.Update()

    return clean_filter.GetOutput()


def make_tubes(points_array, diameter_array):
    polydata_array = []
    for i in range(len(points_array)):
        if len(points_array[i]) > 1:
            points = vtk.vtkPoints()
            diameters = []

            for j in range(len(points_array[i])):
                if diameter_array[points_array[i][j][0]][points_array[i][j][1]][points_array[i][j][2]] < 0.3:
                    diameter_array[points_array[i][j][0]][points_array[i][j][1]][points_array[i][j][2]] = 0.3

                diameter_array[points_array[i][j][0]][points_array[i][j][1]][points_array[i][j][2]] *= 0.5
                points.InsertNextPoint(points_array[i][j])
                diameters.append(diameter_array[points_array[i][j][0]][points_array[i][j][1]][points_array[i][j][2]])

            if points.GetNumberOfPoints() > 1:
                tmp_polydata = make_polydata(points, diameters)
                polydata_array.append(tmp_polydata)

    polydata = assembly(polydata_array)
    tube_actor, line_actor = make_actors(polydata)

    return tube_actor, line_actor


def render(actors):
    renderer = vtk.vtkRenderer()

    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
    interactor.SetRenderWindow(render_window)

    for actor in actors:
        actor.GetProperty().SetRepresentationToWireframe()
        renderer.AddActor(actor)

    renderer.SetBackground(0.1, 0.2, 0.4)
    render_window.SetSize(800, 800)

    interactor.Initialize()

    renderer.ResetCamera()
    renderer.GetActiveCamera().Zoom(1.5)
    render_window.Render()

    interactor.Start()

# ============================================================================================7

skeleton_filenames = []
diameter_filenames = []


read_folder = "path"
write_folder = folder = "path"

for i in range(2, 12):
    skeleton_filenames.append("filename")
    diameter_filenames.append("filename")

for i in range(len(skeleton_filenames)):
    img = nib.load(read_folder + skeleton_filenames[i] + ".nii.gz")
    a = np.array(img.dataobj)

    img = nib.load(read_folder + diameter_filenames[i] + ".nii.gz")
    r = np.array(img.dataobj)

    tube_actor, line_actor = make_skeleton(a, r)

    write_vtk(tube_actor, write_folder + diameter_filenames[i])
    write_vtk(line_actor, write_folder + skeleton_filenames[i])

    #render([tube_actor, line_actor])

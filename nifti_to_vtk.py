import vtk
import numpy as np
from scipy.interpolate import interp1d, UnivariateSpline
import nibabel as nib
import matplotlib.pyplot as plt

# ============================================================================


def write_vtk(polydata, filename):
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(filename + ".vtk")
    writer.SetInputData(polydata)
    writer.Write()


def make_actors(vessels, dti):
    tuber = vtk.vtkTubeFilter()
    tuber.SetInputData(vessels)
    tuber.SetNumberOfSides(6)
    tuber.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
    tuber.Update()

    # Setup actors and mappers
    skeleton_mapper = vtk.vtkPolyDataMapper()
    skeleton_mapper.SetInputData(vessels)
    skeleton_mapper.SetScalarRange(vessels.GetScalarRange())

    vessel_mapper = vtk.vtkPolyDataMapper()
    vessel_mapper.SetInputConnection(tuber.GetOutputPort())
    vessel_mapper.SetScalarRange(vessels.GetScalarRange())

    skeleton_actor = vtk.vtkActor()
    skeleton_actor.SetMapper(skeleton_mapper)

    vessel_actor = vtk.vtkActor()
    vessel_actor.SetMapper(vessel_mapper)

    dti_mapper = vtk.vtkPolyDataMapper()
    dti_mapper.SetInputData(dti)
    dti_mapper.Update()

    dti_actor = vtk.vtkActor()
    dti_actor.SetMapper(dti_mapper)

    return vessel_actor, skeleton_actor, dti_actor


def make_axis_actor():
    points = vtk.vtkPoints()
    points.InsertNextPoint([0, 0, 0])  # 0: origo
    points.InsertNextPoint([100, 0, 0])  # 1: x
    points.InsertNextPoint([0, 100, 0])  # 2: y
    points.InsertNextPoint([0, 0, 100])  # 3: z

    x = vtk.vtkLine()
    x.GetPointIds().SetId(0, 0)
    x.GetPointIds().SetId(1, 1)

    y = vtk.vtkLine()
    y.GetPointIds().SetId(0, 0)
    y.GetPointIds().SetId(1, 2)

    z = vtk.vtkLine()
    z.GetPointIds().SetId(0, 0)
    z.GetPointIds().SetId(1, 3)

    cells = vtk.vtkCellArray()
    cells.InsertNextCell(x)
    cells.InsertNextCell(y)
    cells.InsertNextCell(z)

    colors = vtk.vtkUnsignedCharArray()
    colors.SetNumberOfComponents(3)
    colors.InsertNextTuple3(255, 0, 0)
    colors.InsertNextTuple3(0, 255, 0)
    colors.InsertNextTuple3(0, 0, 255)

    lines = vtk.vtkPolyData()
    lines.SetPoints(points)
    lines.SetLines(cells)
    lines.GetCellData().SetScalars(colors)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(lines)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    return actor


def render(actors):
    renderer = vtk.vtkRenderer()

    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
    interactor.SetRenderWindow(render_window)

    actors.append(make_axis_actor())

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


# ============================= STATS ==========================================

def make_stats(skeleton_array, diameter_array, sum_angles):
    length = []
    avg_diameter = []
    for line in skeleton_array:
        sum_diameter = 0
        for point in line:
            sum_diameter += diameter_array[point[0]][point[1]][point[2]]
        length.append(len(line))
        avg_diameter.append(sum_diameter / len(line))

    plt.figure(1)
    plt.subplot(221)
    plt.xlabel('Segment lengths')
    plt.ylabel('Average diameter in segment')
    plt.plot(length, avg_diameter, 'bo')

    plt.subplot(222)
    plt.xlabel('Number of segment')
    plt.xlabel('Length of segment')
    plt.plot(np.sort(length))

    plt.subplot(223)
    plt.ylabel('Sum difference from x (red)(left to right), y (green)(forward), z (blue)(up) axis')
    plt.plot(0, sum_angles[0], 'r*', 1, sum_angles[1], 'g*', 2, sum_angles[2], 'b*')
    plt.show()


# =========================== GENERATE ARRAYS =================================


def make_skeleton(points_array):
    look_up = np.zeros(points_array.shape)
    look_up = look_up - 1

    skeleton_points = []

    for i in range(points_array.shape[0]):
        for j in range(points_array.shape[1]):
            for k in range(points_array.shape[2]):
                if points_array[i, j, k] == 2:
                    print("___________")
                    end_point = neighbor_finder(points_array, look_up, skeleton_points, i, j, k)

                    if points_array[i, j, k] == 1 or end_point == [i, j, k]:
                        end_finder(end_point, points_array, skeleton_points)

    print("======================")

    return skeleton_points


def neighbor_finder(array, look_up, skeleton_array, i, j, k):
    cont = True
    points = []
    while cont:
        cont = False
        if array[i, j, k] > 0:  # and look_up[i, j, k] == -1:
            if look_up[i, j, k] == -1:
                print("new point: ", i, j, k)
                points.append([i, j, k])
                result = new_coord(array, look_up, i, j, k)
                look_up[i, j, k] = 1
                if len(result) == 3:
                    i = result[0]
                    j = result[1]
                    k = result[2]
                    cont = True
            else:
                points.append([i, j, k])

    skeleton_array.append(points)
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


def in_bounds(array, x, y, z):
    if x >= 0 and y >= 0 and z >= 0 and x < array.shape[0] and y < array.shape[1] and z < array.shape[2]:
        return True
    return False


def new_coord(array, look_up, i, j, k):
    result = []
    for x in range(-1, 2):
        for y in range(-1, 2):
            for z in range(-1, 2):
                if (not (x == 0 and y == 0 and z == 0)) and in_bounds(array, i + x, j + y, k + z):
                    if look_up[i + x, j + y, k + z] == -1 and array[i + x, j + y, k + z] > 0:
                        result.append(i+x)
                        result.append(j+y)
                        result.append(k+z)
                        return result

    return result


# ============================================================


def make_tubes(points_array, diameter_array):
    clean_points_array = []
    thresholded_diameters_array = []
    colors_array = []

    x_angles = 0
    y_angles = 0
    z_angles = 0

    for i in range(len(points_array)):
        diameters = []
        colors = []
        points = []
        if len(points_array[i]) > 3:
            for j in range(len(points_array[i])):
                if diameter_array[points_array[i][j][0]][points_array[i][j][1]][points_array[i][j][2]] < 0.3:
                    diameter_array[points_array[i][j][0]][points_array[i][j][1]][points_array[i][j][2]] = 0.3

                diameter_array[points_array[i][j][0]][points_array[i][j][1]][points_array[i][j][2]] *= 0.5
                diameters.append(diameter_array[points_array[i][j][0]][points_array[i][j][1]][points_array[i][j][2]])

                points.append(points_array[i][j])

                if j > 0:
                    angles = calculate_angles(points_array[i][j - 1], points_array[i][j])
                else:
                    angles = calculate_angles(points_array[i][j], points_array[i][j + 1])

                rgb = angle_to_rgb(angles)
                x_angles += angles[0]
                y_angles += angles[1]
                z_angles += angles[2]
                colors.append(rgb)
        clean_points_array.append(points)
        thresholded_diameters_array.append(diameters)
        colors_array.append(colors)

    smoothed_lines = []
    interpolated_diameters = []
    interpolated_colors = []
    for line in clean_points_array:
        if len(line) > 3:
            smoothed_line = []
            transposed_points = np.transpose(line)
            for dim in transposed_points:
                smoothed_line.append(fit_univariated_spline_on_points(dim))
            smoothed_lines.append(np.transpose(smoothed_line))
    for diameter in thresholded_diameters_array:
        if len(diameter) > 3:
            interpolated_diameters.append(interpolate_points(diameter))
    for color in colors_array:
        if len(color) > 3:
            interpolated_color = []
            transposed_color = np.transpose(color)
            for component in transposed_color:
                interpolated_color.append(interpolate_points(component))
            interpolated_colors.append(np.transpose(interpolated_color))

    vessels_array = []
    dti_array = []

    for i in range(len(smoothed_lines)):
        points = vtk.vtkPoints()
        cells = vtk.vtkCellArray()
        cells.InsertNextCell(len(smoothed_lines[i]))

        diameters = vtk.vtkDoubleArray()
        diameters.SetNumberOfTuples(len(smoothed_lines[i]))
        diameters.SetName("VesselDiameter")

        colors = vtk.vtkUnsignedCharArray()
        colors.SetNumberOfComponents(3)
        colors.SetName("Colors")

        for j in range(len(smoothed_lines[i])):
            points.InsertNextPoint(smoothed_lines[i][j])
            cells.InsertCellPoint(j)
            diameters.SetTuple1(j, interpolated_diameters[i][j])
            colors.InsertNextTuple3(interpolated_colors[i][j][0], interpolated_colors[i][j][1], interpolated_colors[i][j][2])

        tmp_line = vtk.vtkPolyData()
        tmp_line.SetPoints(points)
        tmp_line.SetLines(cells)

        tmp_vessels = vtk.vtkPolyData()
        tmp_vessels.DeepCopy(tmp_line)
        tmp_vessels.GetPointData().AddArray(diameters)
        tmp_vessels.GetPointData().SetActiveScalars("VesselDiameter")

        tmp_dti = vtk.vtkPolyData()
        tmp_dti.DeepCopy(tmp_line)
        tmp_dti.GetPointData().SetScalars(colors)

        vessels_array.append(tmp_vessels)
        dti_array.append(tmp_dti)

    vessels = assembly(vessels_array)
    dti = assembly(dti_array)

    return vessels, dti, [x_angles, y_angles, z_angles]


def fit_univariated_spline_on_points(points):
    x = np.linspace(0, 1, num=len(points), endpoint=True)
    spline = UnivariateSpline(x, points)
    new_points = spline(np.linspace(0, 1, len(points)*10))
    return new_points


def interpolate_points(points):
    x = np.linspace(0, 1, num=len(points), endpoint=True)

    f = interp1d(x, points)
    return f((np.linspace(0, 1, num=len(points)*10, endpoint=True)))


# ===========================================================


def assembly(polydata_array):
    append_filter = vtk.vtkAppendPolyData()
    for polydata in polydata_array:
        append_filter.AddInputData(polydata)
    append_filter.Update()

    clean_filter = vtk.vtkCleanPolyData()
    clean_filter.SetInputConnection(append_filter.GetOutputPort())
    clean_filter.Update()

    return clean_filter.GetOutput()


# ===========================================================

def magnitude(p):
    summa = 0
    for x in p:
        summa += x * x
    return np.sqrt(summa)


def dot_product(u, v):
    summa = 0
    for i in range(len(u)):
        summa += u[i] * v[i]
    return summa


def calculate_angles(p1, p2):
    vector = []
    for i in range(len(p1)):
        vector.append(p1[i] - p2[i])
    theta_x = np.rad2deg((np.arccos(dot_product(vector, [1, 0, 0]) / magnitude(vector))) / 2)
    theta_y = np.rad2deg((np.arccos(dot_product(vector, [0, 1, 0]) / magnitude(vector))) / 2)
    theta_z = np.rad2deg((np.arccos(dot_product(vector, [0, 0, 1]) / magnitude(vector))) / 2)
    return [theta_x, theta_y, theta_z]


def angle_to_rgb(angles):
    new = []
    for angle in angles:
        new.append(int((angle / 90) * 255))
    return new


def scale(value, max_current, max_new=1):
    return (value / max_current) * max_new


# ============================================================


skeleton_filenames = []
diameter_filenames = []

read_folder = "read/folder/"
write_folder = "write/folder/"

for n in range(2, 12):
    skeleton_filenames.append("skeleton_files")
    diameter_filenames.append("diameter_files")

for n in range(len(skeleton_filenames)):
    img = nib.load(read_folder + skeleton_filenames[n] + ".nii.gz")
    skeleton_data = np.array(img.dataobj)

    img = nib.load(read_folder + diameter_filenames[n] + ".nii.gz")
    diameter_data = np.array(img.dataobj)

    skeleton = make_skeleton(skeleton_data)

    vessels_polydata, dti_polydata, sum_angles = make_tubes(skeleton, diameter_data)

    write_vtk(vessels_polydata, write_folder + diameter_filenames[n] + "_approx_spline")
    write_vtk(dti_polydata, write_folder + skeleton_filenames[n] + "_approx_spline_dti")

    vessels_actor, skeleton_actor, dti_actor = make_actors(vessels_polydata, dti_polydata)

    make_stats(skeleton, diameter_data, sum_angles)

    render([vessels_actor, skeleton_actor, dti_actor])


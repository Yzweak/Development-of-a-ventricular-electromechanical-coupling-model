import vtk

mesh_reader = vtk.vtkSTLReader()
mesh_reader.SetFileName("left_ventricle.stl")
mesh_reader.Update()
mesh = mesh_reader.GetOutput()

data = mesh
zoffset = 55.0
height =  data.GetBounds()[5] - zoffset
plane = vtk.vtkPlane()
plane.SetNormal(0, 0, -1)
plane.SetOrigin(0, 0, height)

clip = vtk.vtkClipPolyData()
clip.SetClipFunction(plane)
if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
    clip.SetInputData(data)
else:
    clip.SetInput(data)
clip.Update()
clipped_data = clip.GetOutput()

smoother = vtk.vtkSmoothPolyDataFilter()
smoother.SetInputData(clipped_data)
smoother.SetNumberOfIterations(20)  # 迭代次数，可以根据需要进行调整
smoother.SetRelaxationFactor(0.1)  # 平滑参数，可以根据需要进行调整

smoother.Update()
smoothed_polydata = smoother.GetOutput()

stl_writer = vtk.vtkSTLWriter()
stl_writer.SetFileName("cliped_data.stl")
if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
    stl_writer.SetInputData(smoothed_polydata)
else:
    stl_writer.SetInput(smoothed_polydata)
stl_writer.Update()
stl_writer.Write()

mapper = vtk.vtkPolyDataMapper()
mapper.SetInputData(smoothed_polydata)
actor = vtk.vtkActor()
actor.SetMapper(mapper)

renderer = vtk.vtkRenderer()
renderer.AddActor(actor)

render_window = vtk.vtkRenderWindow()
render_window.AddRenderer(renderer)
render_window_interactor = vtk.vtkRenderWindowInteractor()
render_window_interactor.SetRenderWindow(render_window)
render_window_interactor.Start()
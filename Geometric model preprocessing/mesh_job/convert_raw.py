import vtk

# Read 3D RAW image
reader = vtk.vtkImageReader()
reader.SetDataScalarType(vtk.VTK_UNSIGNED_CHAR)  # unsigned int8
reader.SetFileName('left_ventricle.raw')
reader.SetNumberOfScalarComponents(1)
reader.SetFileDimensionality(3)
reader.SetDataByteOrderToLittleEndian()
reader.SetDataExtent(0, 324, 0, 324, 0, 424)  # image size 488*488*332
reader.SetDataSpacing(1.0, 1.0, 1.0)  # Volume Pixel
reader.Update()

# Visualization
contour=vtk.vtkMarchingCubes()  # vtk.vtkContourFilter()
contour.SetInputConnection(reader.GetOutputPort())
contour.ComputeNormalsOn()
contour.SetValue(0,1)

mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(contour.GetOutputPort())
mapper.ScalarVisibilityOff()

actor = vtk.vtkActor()
actor.SetMapper(mapper)

renderer=vtk.vtkRenderer()
renderer.SetBackground([0.5, 0.5, 0.5])
renderer.AddActor(actor)

window = vtk.vtkRenderWindow()
window.SetSize(640, 640)
window.AddRenderer(renderer)

# Create interactor, add window & add observers
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(window)

# Start renderer & interactor
window.Render()
interactor.Initialize()
interactor.Start()


stl_writer = vtk.vtkSTLWriter()
stl_writer.SetFileName("left_ventricle.stl")
if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
    stl_writer.SetInputData(contour.GetOutput())
else:
    stl_writer.SetInput(contour.GetOutput())
stl_writer.Update()
stl_writer.Write()

# # Write in STL
# triangle = vtk.vtkTriangleFilter()
# triangle.SetInputConnection(contour.GetOutputPort())
# triangle.PassVertsOff()
# triangle.PassLinesOff()
#
# decimation=vtk.vtkQuadricDecimation()
# decimation.SetInputConnection(triangle.GetOutputPort())
#
# clean=vtk.vtkCleanPolyData()
# clean.SetInputConnection(triangle.GetOutputPort())
#
# triangle2 = vtk.vtkTriangleFilter()
# triangle2.SetInputConnection(clean.GetOutputPort())
# triangle2.PassVertsOff()
# triangle2.PassLinesOff()
#
# stlWriter = vtk.vtkSTLWriter()
# stlWriter.SetInputConnection(triangle2.GetOutputPort())
# stlWriter.SetFileName("left_ventricle.stl")
# stlWriter.Write()
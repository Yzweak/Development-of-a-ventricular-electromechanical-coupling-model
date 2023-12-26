# @Author: charlesmann
# @Date:   2021-09-20T19:22:52-04:00
# @Last modified by:   charlesmann
# @Last modified time: 2021-10-20T13:17:58-04:00

from __future__ import division
import sys
sys.path.append("/home/fenics/shared/Myosim/dependencies/")
sys.path.append("/home/fenics/shared/Myosim/source_code/")
import os as os
from dolfin import *
import numpy as np
from forms import Forms
from nsolver import NSolver as NSolver
import math
import Python_MyoSim.half_sarcomere.half_sarcomere as half_sarcomere
import Python_MyoSim.half_sarcomere.implement as implement
from cell_ion_module import cell_ion_driver
from edgetypebc import *
import pandas as pd
import copy
from methods import mesh_import
from methods.mesh_import import mesh_import as mesh_import
from methods.assign_initial_hsl import assign_initial_hsl as assign_hsl
from methods.assign_local_coordinate_system import assign_local_coordinate_system as lcs
from methods.assign_heterogeneous_params import assign_heterogeneous_params as assign_params
from methods.assign_heterogeneous_params import initialize_dolfin_functions as initialize_dolfin_functions
from methods.set_boundary_conditions import set_bcs as set_bcs
from methods.circulatory_module import circulatory_module as cm
from methods.update_boundary_conditions import update_boundary_conditions
from methods.save_solution import save_solution
from methods.load_solution import load_solution as load_sol
from methods.grow_mesh import grow_mesh
import recode_dictionary
import json
import timeit
import scipy.interpolate as interpol

# f = open("myprint.txt","w+")
# sys.stdout = f

def fenics(sim_params):
    # declare these as indices for things like as_tensor()
    i,j = indices(2)
    m,k = indices(2)
    ## Assign simulation parameters
    sim_geometry = sim_params["simulation_geometry"][0] #unit_cube, cylinder, ventricle
    if sim_geometry == "unit_cube":
        geo_options = {}
    else:
        geo_options = sim_params["geometry_options"]

    sim_protocol = sim_params["protocol"] # contains simulation dependent options
    sim_timestep = sim_protocol["simulation_timestep"][0]
    sim_duration = sim_protocol["simulation_duration"][0] # can be overwritten in ventricle protocol once number of cycles and heartrate specified
    no_of_time_steps = int(sim_duration/sim_timestep)
    t = np.linspace(0,sim_duration,no_of_time_steps)
    save_cell_output = sim_params["save_cell_output"][0] # myosim output
    save_visual_output = sim_params["save_visual_output"][0] # paraview files for visualization
    if "load_solution" in sim_params.keys():
	load_solution = sim_params["load_solution"][0]
	load_solution_dir = sim_params["load_solution"][1]
    else:
	load_solution = 0
    if "save_solution" in sim_params.keys():
	save_solution_flag = sim_params["save_solution"][0]
    else:
	save_solution_flag = 0
    output_path = sim_params["output_path"][0]
    print "output path: ", output_path

    # assign amount of random variation in f0 (cube and cylinder simulations, 0 means normal alignment)
    gaussian_width = sim_params["fiber_orientation"]["fiber_randomness"][0]

    if "save_fx_only" in sim_params.keys():
        save_fx_only = sim_params["save_fx_only"][0]
    else:
        save_fx_only = 0
    if 'growth_params' in globals():
        print "loading growth params"
        if "eccentric_growth" in growth_params.keys():
            ecc_growth_rate = growth_params["eccentric_growth"]["time_constant"][0]
            set_point = growth_params["eccentric_growth"]["passive_set_point"][0]
            k_myo_damp = Constant(growth_params["eccentric_growth"]["k_myo_damp"][0])
        if "fiber_reorientation" in growth_params.keys():
            print "ASSIGNING FIBER REMODELING LAW PARAMS"
            ordering_law = growth_params["fiber_reorientation"]["law"][0]
            kroon_time_constant = growth_params["fiber_reorientation"]["time_constant"][0]
            print "KROON TIME CONSTANT =",kroon_time_constant
            reorient_start_time = growth_params["fiber_reorientation"]["reorient_t_start"][0]
            stress_name = growth_params["fiber_reorientation"]["stress_type"][0]
    else:
        k_myo_damp = 0.0


#------------------------------------------------------------------------------
#           Initialize myosim info for active stress calculation
#------------------------------------------------------------------------------
    fiber_state = hs_params["myofilament_parameters"]["fiber_state"][0]
    no_of_states = hs_params["myofilament_parameters"]["num_states"][0]
    N0 = float(hs_params["myofilament_parameters"]["N0"][0])
    N1 = float(hs_params["myofilament_parameters"]["N1"][0])
    P0 = float(hs_params["myofilament_parameters"]["P0"][0])
    P1 = float(hs_params["myofilament_parameters"]["P1"][0])
    P2 = float(hs_params["myofilament_parameters"]["P2"][0])
    P3 = float(hs_params["myofilament_parameters"]["P3"][0])
    factor = float(hs_params["myofilament_parameters"]["force_factor"][0])
    LTRPNCa = float(hs_params["myofilament_parameters"]["LTRPNCa"][0])

    xfiber_fraction = hs_params["myofilament_parameters"]["xfiber_fraction"][0]
    n_array_length = no_of_states

# ------------------------------------------------------------------------------
#           Mesh Information
# ------------------------------------------------------------------------------
    mesh, lv_options = mesh_import.import_mesh(sim_geometry, geo_options)
    facetboundaries = MeshFunction('size_t', mesh, mesh.topology().dim() - 1)
    edgeboundaries = MeshFunction('size_t', mesh, mesh.topology().dim() - 2)
    subdomains = MeshFunction('int', mesh, 3)
    # define communicator, for running with multiple cores in parallel
    comm = mesh.mpi_comm()
    no_of_cells = len(subdomains.array())
    # from the mesh, define some things
    if sim_geometry == "cylinder" or sim_geometry == "unit_cube" or sim_geometry == "box_mesh" or sim_geometry == "gmesh_cylinder":
        no_of_int_points = 4 * np.shape(mesh.cells())[0]
        deg = 2
        ds = dolfin.ds(subdomain_data = facetboundaries)
        fx_rxn = np.zeros((no_of_time_steps))
    else:
        #ventricle modeling
        # Kurtis changing this 6/17/21 to test_1.5 accuracy/time
        #deg = 4
        #no_of_int_points = 14 * np.shape(mesh.cells())[0]
        deg = 2
        no_of_int_points = 4 * np.shape(mesh.cells())[0]
        #set surface id numbers
        topid = 4
        LVendoid = 2
        epiid = 1

    parameters["form_compiler"]["quadrature_degree"]=deg
    parameters["form_compiler"]["representation"] = "quadrature"
    #edgeboundaries = MeshFunction('size_t', mesh, mesh.topology().dim()-2)
    N = FacetNormal (mesh)
    X = SpatialCoordinate (mesh)
    dx = dolfin.dx(mesh,metadata = {"integration_order":2})
    isincomp = True
    # periodic boundary condition
    class PeriodicBoundary(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1E-14
            return on_boundary and abs(x[1]) < tol
        def map(self, x, y):
            # map coordinate from x on top to y on bottom
            y[0] = x[0]
            y[1] = x[1]-1.0
            y[2] = x[2]

    pbc = PeriodicBoundary()

#------------------------------------------------------------------------------
#           Storage data frames
#------------------------------------------------------------------------------
    ## Create files for saving information if needed
    output_file = XDMFFile(mpi_comm_world(), output_path + "output.xdmf")
    output_file.parameters.update({"functions_share_mesh": True, "rewrite_function_mesh": False, "flush_output": True})

    fiber_file = XDMFFile(mpi_comm_world(), output_path + "fiber.xdmf")
    fiber_file.parameters.update({"functions_share_mesh": True, "rewrite_function_mesh": False, "flush_output": True})
    ## Create files for saving information if needed
    if save_visual_output:
        f0_vs_time_array = np.zeros((no_of_int_points, 3, no_of_time_steps))
        # stress visualization
        pk2_passive_file = File(output_path + "pk2_passive.pvd")
        f0_deformed_file = File(output_path + "f0_deformed.pvd")
        shearfs_file = File(output_path + "shear_fs.pvd")
        shearfn_file = File(output_path + "shear_fn.pvd")
        shearsn_file = File(output_path + "shear_sn.pvd")
        stress_eigen_ds = pd.DataFrame(np.zeros((no_of_int_points, 3)), index=None)
        f_adjusted_ds = pd.DataFrame(np.zeros((no_of_int_points, 3)), index=None)

    if save_cell_output:
        # active stress
        active_stress_ds = pd.DataFrame(np.zeros(no_of_int_points),index=None)
        active_stress_ds = active_stress_ds.transpose()

        # saving f0 vectors
        f0_components = np.empty((no_of_int_points,3,no_of_time_steps))
        f0_components[:] = np.NaN

        calcium_ds = pd.DataFrame(np.zeros(no_of_int_points),index=None)
        calcium_ds = calcium_ds.transpose()

        # myosim populations
        dumped_populations_ds = pd.DataFrame(np.zeros((no_of_int_points,n_array_length)))
        dumped_populations = np.zeros((no_of_int_points,n_array_length))

        # passive stressed in myofiber
        p_f_array_ds = pd.DataFrame(np.zeros(no_of_int_points),index=None)
        p_f_array_ds = p_f_array_ds.transpose()

        # guccione bulk passive stress in fiber direction
        pgf_array_ds = pd.DataFrame(np.zeros(no_of_int_points),index=None)
        pgf_array_ds = pgf_array_ds.transpose()

        # guccione bulk passive stress in transverse directions
        pgt_array_ds = pd.DataFrame(np.zeros(no_of_int_points),index=None)
        pgt_array_ds = pgt_array_ds.transpose()

        # guccione bulk passive shear stress
        pgs_array_ds =pd.DataFrame(np.zeros(no_of_int_points),index=None)
        pgs_array_ds = pgs_array_ds.transpose()

        # thick and thin filament overlaps
        temp_overlap_ds = pd.DataFrame(np.zeros(no_of_int_points),index=None)
        temp_overlap_ds = temp_overlap_ds.transpose()

        # stretch
        alpha_array_ds = pd.DataFrame(np.zeros(no_of_int_points),index=None)
        alpha_array_ds = alpha_array_ds.transpose()

        # half-sarcomere length
        hsl_array_ds =pd.DataFrame(np.zeros(no_of_int_points),index=None)
        hsl_array_ds = hsl_array_ds.transpose()

        # change in half-sarcomere length
        delta_hsl_array_ds = pd.DataFrame(np.zeros(no_of_int_points),index=None)
        delta_hsl_array_ds = delta_hsl_array_ds.transpose()

    # Initialize some data holders that are necessary
    temp_overlap = np.zeros((no_of_int_points))
    calcium = np.zeros((no_of_time_steps, 1))
    rxn_force = np.zeros(no_of_time_steps + 21)
    avg_fdiff_x_array = np.zeros(no_of_time_steps)
    avg_fdiff_y_array = np.zeros(no_of_time_steps)
    avg_fdiff_z_array = np.zeros(no_of_time_steps)
    reorient_finish_array = np.zeros(no_of_time_steps)
    avg_fdiff_norm_array = np.zeros(no_of_time_steps)
    delta_hsl_array = np.zeros(no_of_int_points)
    traction_switch_flag = 0


    # Put any needed expressions here
    #--------------------------------
    # initialize LV cavity volume, only updated if windkessel is called
    if load_solution > 0:
        LVCavityvol = expressions_loaded["LVCavityvol"]
        u_D = expressions_loaded["u_D"]
        u_top = expressions_loaded["u_top"]
        u_front = expressions_loaded["u_front"]
        Press = expressions_loaded["Press"]
    else:
        LVCavityvol = Expression(("vol"), vol=0.0, degree=2)
        # displacement boundary expression for end of cell or fiber sims
        u_D = Expression(("u_D"), u_D = 0.0, degree = 0)
        # Forcing volume preserving for biaxial case
        u_top = Expression(("u_top"), u_top = 0.0, degree = 0)
        u_front = Expression(("u_front"), u_front = 0.0, degree = 0)
        # traction boundary condition for end of cell/fiber, could use this to apply a traction to epicardium or something
        Press = Expression(("P"), P=0.0, degree=0)

    sim_protocol["start_diastolic_pressure"] = Press.P

    end_systole = 0
    end_diastole = 0

    expressions = {
        "u_D":u_D,
        "u_top":u_top,
        "u_front":u_front,
        "Press":Press
    }

    if sim_protocol["simulation_type"][0] == "custom":
        custom_disp_array = np.load("./custom_displacement.npy")
    else:
        custom_disp_array = []

#-------------------------------------------------------------------------------
#           Initialize finite elements and function spaces
#-------------------------------------------------------------------------------

    # Vector element at gauss points (for fibers)
    VQuadelem = VectorElement("Quadrature", mesh.ufl_cell(), degree=deg, quad_scheme="default")
    VQuadelem._quad_scheme = 'default'

    # General quadrature element whose points we will evaluate myosim at
    Quadelem = FiniteElement("Quadrature", tetrahedron, degree=deg, quad_scheme="default")
    Quadelem._quad_scheme = 'default'

    # Vector element for displacement
    Velem = VectorElement("CG", mesh.ufl_cell(), 2, quad_scheme="default")
    Velem._quad_scheme = 'default'

    # Quadrature element for pressure
    Qelem = FiniteElement("CG", mesh.ufl_cell(), 1, quad_scheme="default")
    Qelem._quad_scheme = 'default'

    # Real element for rigid body motion boundary condition
    Relem = FiniteElement("Real", mesh.ufl_cell(), 0, quad_scheme="default")
    Relem._quad_scheme = 'default'

    Telem2 = TensorElement("Quadrature", mesh.ufl_cell(), degree=deg, shape=2*(3,), quad_scheme='default')
    Telem2._quad_scheme = 'default'
    for e in Telem2.sub_elements():
    	e._quad_scheme = 'default'

    # Mixed element for rigid body motion. One each for x, y displacement. One each for x, y, z rotation
    VRelem = MixedElement([Relem, Relem, Relem, Relem, Relem])


    # ------- Define function spaces on mesh using above elements --------------
    # Quadrature space for information needed at gauss points, such as hsl, cb_force, passive forces, etc.
    Quad = FunctionSpace(mesh, Quadelem)
    # go ahead and get coordinates of quadrature points
    gdim = mesh.geometry().dim()
    xq = Quad.tabulate_dof_coordinates().reshape((-1,gdim))
    np.save(output_path + 'quadrature_dof',xq)

    # Function space for myosim populations
    Quad_vectorized_Fspace = FunctionSpace(mesh, MixedElement(n_array_length*[Quadelem]))
    # Function space for local coordinate system (fiber, sheet, sheet-normal)
    fiberFS = FunctionSpace(mesh, VQuadelem)
    # Tensor function space
    TF = TensorFunctionSpace(mesh, 'DG', 1)
    TFQuad = FunctionSpace(mesh, Telem2)
    TF_kroon = TensorFunctionSpace(mesh,'DG',1)

    # Initializing functions to hold growth stimuli
    if load_solution > 0:
        concentric_growth_stimulus = functions_loaded["concentric_growth_stimulus"]
        eccentric_growth_stimulus = functions_loaded["eccentric_growth_stimulus"]
        W = functions_loaded["W"]
        x_dofs = W.sub(0).sub(0).dofmap().dofs() # will use this for x rxn forces later
    else:
        concentric_growth_stimulus = Function(FunctionSpace(mesh, "DG",1))
        eccentric_growth_stimulus = Function(FunctionSpace(mesh,"DG",1))
        if sim_geometry == "cylinder" or sim_geometry == "unit_cube" or sim_geometry == "box_mesh" or sim_geometry == "gmesh_cylinder":
            if sim_protocol["simulation_type"][0] == "ramp_and_hold_simple_shear":
                print "implementing periodic boundary condition"
                W = FunctionSpace(mesh, MixedElement([Velem,Qelem]),constrained_domain=PeriodicBoundary())
            else:
                #W = FunctionSpace(mesh, MixedElement([Velem,Qelem,Relem]))
                W = FunctionSpace(mesh, MixedElement([Velem,Qelem]))
            x_dofs = W.sub(0).sub(0).dofmap().dofs() # will use this for x rxn forces later
            #print "assigned W"
        else:
            W = FunctionSpace(mesh, MixedElement([Velem,Qelem,Relem,VRelem]))

    marker_space = FunctionSpace(mesh, 'CG', 1)

#-------------------------------------------------------------------------------
#           Initialize functions on the above spaces
#-------------------------------------------------------------------------------

    # fiber, sheet, and sheet-normal functions
    if load_solution > 0:
        f0 = functions_loaded["f0"]
        s0 = functions_loaded["s0"]
        n0 = functions_loaded["n0"]
    else:
        f0 = Function(fiberFS)
        s0 = Function(fiberFS)
        n0 = Function(fiberFS)

    x_dir = Function(VectorFunctionSpace(mesh,"CG",1))
    x_vec = Function(fiberFS)
    x2_shape = np.shape(x_dir.vector().get_local())

    for jj in np.arange(int(x2_shape[0]/3)):
        x_dir.vector()[jj*3] = 1.
        x_dir.vector()[jj*3+1] = 0.
        x_dir.vector()[jj*3+2] = 0.
    for jj in np.arange(no_of_int_points):
        x_vec.vector()[jj*3] = 1.
        x_vec.vector()[jj*3+1] = 0.
        x_vec.vector()[jj*3+1] = 0.

    # put these in a dictionary to pass to function for assignment
    coord_params = {
        "f0":f0,
        "s0":s0,
        "n0":n0,
        "fiberFS":fiberFS,
        "marker_space":marker_space,
        "sim_geometry":sim_geometry,
        "mesh":mesh,
        "Quad":Quad,
        "no_of_int_points":no_of_int_points,
        "geo_options":geo_options,
        "facetboundaries":facetboundaries,
        "edgeboundaries":edgeboundaries
    }

    # create lists of dictionaries that hold parameters for each gauss point
    # These will remain scalars. I anticipate this will be the case for
    # any parameters that will be passed into existing modules (myosim parameters,
    # possible calcium models, etc)
    hs_params_list = [{}]*no_of_int_points
    #passive_params_list = [{}]*no_of_int_points

    # Must make a deep copy so each item in the list is independent, and not linked
    # to the original paramter dictionary
    for jj in np.arange(np.shape(hs_params_list)[0]):
        hs_params_list[jj] = copy.deepcopy(hs_params) # because this is a copy, everything is initialized
        #passive_params_list[jj] = copy.deepcopy(passive_params)

    # create dictionary of parameters that will be initialized as dolfin functions
    # (anything used by dolfin in calculations will become a function, ex. passive SEF
    # or active stress calculation).
    # dolfin_functions is a nested dictionary
    dolfin_functions = {}
    dolfin_functions["passive_params"] = passive_params
    dolfin_functions["cb_number_density"] = hs_params["cb_number_density"]

    # Initialize all dolfin functions to take on their base value
    # Kurtis putting in conditional. For fiber sims, need coordinates of points to assign heterogeneity
    if (sim_geometry == "gmesh_cylinder"):
        # fiber simulation
        dolfin_functions = initialize_dolfin_functions.initialize_dolfin_functions(dolfin_functions,Quad)
    else:
        #print "element wise assignment of dolfin functions!!!!!!"
        # element wise assignment
        dolfin_functions = initialize_dolfin_functions.initialize_dolfin_functions(dolfin_functions,Quad)


    # functions for the weak form
    if load_solution > 0:
        w = functions_loaded["w"]
    else:
        w = Function(W)

    dw    = TrialFunction(W)
    wtest = TestFunction(W)

    if (sim_geometry == "ventricle") or (sim_geometry == "ellipsoid"):
        # Need the pressure function
        du,dp,dpendo,dc11 = TrialFunctions(W)
        (u,p,pendo,c11)   = split(w)
        (v,q,qendo,v11)   = TestFunctions(W)
        ventricle_params  = {
            "lv_volconst_variable": pendo,
            "lv_constrained_vol":LVCavityvol,
            "LVendoid": LVendoid,
            "LVendo_comp": 2,
        }

    else:
        du,dp = TrialFunctions(W)
        (u,p) = split(w)
        (v,q) = TestFunctions(W)

    # Initial and previous timestep half-sarcomere length functions
    if load_solution > 0:
        hsl0 = functions_loaded["hsl0"]
        y_vec = functions_loaded["y_vec"]
        #print "loaded y_vec, checking here"
        # print y_vec.vector().get_local()
    else:
        hsl0    = Function(Quad)
        y_vec = Function(Quad_vectorized_Fspace)
        hsl0 = assign_hsl.assign_initial_hsl(lv_options,hs_params,sim_geometry,hsl0)
        f0,s0,n0,geo_options = lcs.assign_local_coordinate_system(lv_options,coord_params,sim_params)

    # hdf5write = HDF5File(mesh.mpi_comm(),"fiber_direction.h5","w")
    # hdf5write.write(f0,"f0")
    # return


    hsl_old = Function(Quad)
    y_vec_array_new = y_vec.vector().get_local()[:]
    y_interp = np.zeros(np.shape(y_vec_array_new))

    pseudo_alpha = Function(Quad)
    pseudo_old = Function(Quad)
    pseudo_old.vector()[:] = 1.0
    hsl_diff_from_reference = Function(Quad)
    hsl_diff_from_reference.vector()[:] = 0.0


    y_vec_array = y_vec.vector().get_local()[:]

#-------------------------------------------------------------------------------
#           Assign function values
#-------------------------------------------------------------------------------
    hs_params_list,dolfin_functions = assign_params.assign_heterogeneous_params(sim_params,hs_params,hs_params_list,dolfin_functions,geo_options,no_of_int_points,no_of_cells)
    # Select fibers for visualization (exclude stiff regions)
    binary_mask = np.zeros((no_of_int_points),dtype=int)

    # temp_fcn_visualization = Function(Quad)
    # for mm in np.arange(no_of_int_points):
    #     temp_fcn_visualization.vector()[mm] = hs_params_list[mm]["myofilament_parameters"]["k_force"][0]
    # File(output_path + "k_force.pvd") << project(temp_fcn_visualization,FunctionSpace(mesh,"DG",0))
    # File(output_path + "c_param.pvd") << project(dolfin_functions["passive_params"]["c"][-1],FunctionSpace(mesh,"DG",0))
    # File(output_path + "cb_density.pvd") << project(dolfin_functions["cb_number_density"][-1],FunctionSpace(mesh,"DG",0))
#-------------------------------------------------------------------------------
#           Save initial values
#-------------------------------------------------------------------------------

    #save initial f0, s0, n0, hsl0
    if save_visual_output:
        hsl_temp = project(hsl0,FunctionSpace(mesh,'DG',0))
        hsl_temp.rename("hsl_temp","half-sarcomere length")
        output_file.write(hsl_temp, 0)
#-------------------------------------------------------------------------------
#           Initialize the solver and forms parameters, continuum tensors
#-------------------------------------------------------------------------------

    M1ij = project(as_tensor(f0[m]*f0[k], (m,k)),TF)
    M2ij = project(as_tensor(s0[m]*s0[k], (m,k)),TF)
    M3ij = project(as_tensor(n0[m]*n0[k], (m,k)),TF)

    Theta1 = Function(FunctionSpace(mesh,"DG",1))
    Theta1.vector()[:] = 1.0

    Theta2 = Function(FunctionSpace(mesh,"DG",1))
    Theta2.vector()[:] = 1.0

    Theta3 = Function(FunctionSpace(mesh,"DG",1))
    Theta3.vector()[:] = 1.0

    # Based on the material coordinates, we can define different Growth Tensor Construct
    Fg = Theta1*(M1ij) +  Theta2*M2ij + Theta3*M3ij #always created, only updated if growth


    # parameters for forms file
    params= {"mesh": mesh,
             "facetboundaries": facetboundaries,
             "facet_normal": N,
             "mixedfunctionspace": W,
             "mixedfunction": w,
             "displacement_variable": u,
             "pressure_variable": p,
             "fiber": f0,
             "sheet": s0,
             "sheet-normal": n0,
             "incompressible": isincomp,
             "hsl0": hsl0,
             "Kappa":Constant(1e5),
             "growth_tensor": Fg,
             "M1": M1ij,
             "M2": M2ij,
             "M3": M3ij,
             "TF": TF}


    params.update(dolfin_functions["passive_params"])

    # Initialize the forms module
    uflforms = Forms(params)
    # Get deformation gradient
    Fmat = uflforms.Fmat()
    Fe = uflforms.Fe()
    test_FS = FunctionSpace(mesh,Velem)
    Fmat2 = Function(TF)
    # Get right cauchy stretch tensor
    #Cmat = (Fmat.T*Fmat)
    Cmat = uflforms.Cmat()
    # Get Green strain tensor
    Emat = uflforms.Emat()
    # jacobian of deformation gradient
    J = uflforms.J()
    # facet normal in current config
    n = J*inv(Fmat.T)*N

#-------------------------------------------------------------------------------
#           Initialize boundary conditions
#-------------------------------------------------------------------------------
    # returns a dictionary of bcs and potentially a test_marker_fcn for work loops
    bc_output = set_bcs.set_bcs(sim_geometry,sim_protocol,geo_options,mesh,W,facetboundaries,expressions)
    bcs = bc_output["bcs"]
    bcright = bcs[-1]
    if (sim_protocol["simulation_type"][0] != "cycle"):
        test_marker_fcn = bc_output["test_marker_fcn"]

#-------------------------------------------------------------------------------
#           Active stress calculation
#-------------------------------------------------------------------------------
    # Start with active stress calculation here to validate mmoth_vent, then try
    # to move it to the forms file
    # Calculate a pseudo stretch, not based on deformation gradient
    hsl_old.vector()[:] = hsl0.vector()[:]
    hsl_diff_from_reference = (hsl_old - hsl0)/hsl0

    pseudo_alpha = pseudo_old*(1.-(k_myo_damp*(hsl_diff_from_reference)))
    alpha_f = sqrt(dot(f0, Cmat*f0)) # actual stretch based on deformation gradient
    hsl = pseudo_alpha*alpha_f*hsl0
    delta_hsl = hsl - hsl_old

    cb_force = Constant(0.0)

    y_vec_split = split(y_vec)
    Wp = uflforms.PassiveMatSEF(hsl)

    f01 = 0.15
    f12 = 0.5
    f23 = 0.35
    g01 = 0.1
    g12 = 0.2
    g23 = 0.3
    sumPath = g01 * g12 * g23 + f01 * g12 * g23 + f01 * f12 * g23 + f01 * f12 * f23
    Pmax1 = (f01 * g12 * g23) / sumPath
    Pmax2 = (f01 * f12 * g23) / sumPath
    Pmax3 = (f01 * f12 * f23) / sumPath
    SL = 2.15
    if SL < 2.2:
        alpha_SL = min(1.0, (SL - 2.0 * 1 - 0.1) / (1.5 - 0.1))
    else:
        alpha_SL = 1 - (SL - 2.2) / (1.5 - 0.1)
    Force = 0.1 * alpha_SL * (y_vec_split[3] + y_vec_split[1] + 2 * y_vec_split[4] + 3 * y_vec_split[5]) / (Pmax1 + 2 * Pmax2 + 3 * Pmax3)
    cb_force = Force * dolfin_functions["cb_number_density"][-1] * factor * 1e-8
    Pactive = cb_force * as_tensor(f0[m]*f0[k], (m,k))+ xfiber_fraction*cb_force * as_tensor(s0[m]*s0[k], (m,k))+ xfiber_fraction*cb_force * as_tensor(n0[m]*n0[k], (m,k))
    cb_f_array = project(cb_force, Quad).vector().get_local()[:]
#-------------------------------------------------------------------------------
#           Now hsl function is initiated, make sure all arrays are initialized
#-------------------------------------------------------------------------------
    #create hsl_array from projection
    hsl_array = project(hsl, Quad).vector().get_local()[:]
    #initialize y_vec_array to put all hads in SRX and all binding sites to off
    if load_solution < 1:
        for init_counter in range(0,n_array_length * no_of_int_points,n_array_length):
            y_vec_array[init_counter] = N0
            y_vec_array[init_counter + 1] = N1
            y_vec_array[init_counter + 2] = P0
            y_vec_array[init_counter + 3] = P1
            y_vec_array[init_counter + 4] = P2
            y_vec_array[init_counter + 5] = P3
            y_vec_array[init_counter + 6] = LTRPNCa

    #Get passive stress tensor, magnitude of myofiber passive stress
    PK2_passive,Sff = uflforms.stress(hsl)
    # need p_f_array for myosim
    temp_DG = project(Sff, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
    p_f = interpolate(temp_DG, Quad)
    p_f_array = p_f.vector().get_local()[:]


    F1 = derivative(Wp, w, wtest)*dx
    # active stress contribution (Pactive is PK2, transform to PK1)
    F2 = inner(Fmat*Pactive, grad(v))*dx
    F3 = inner(Press * N, v) * ds(2, domain=mesh)
    Ftotal = F1 + F2 - F3

    Jac1 = derivative(F1, w, dw)
    Jac2 = derivative(F2, w, dw)
    Jac3 = derivative(F3, w, dw)
    Jac = Jac1 + Jac2 - Jac3

    # Can use Dr. Lee's Nsolver if solver needs debugging
    solverparams = {"Jacobian": Jac,
                    "F": Ftotal,
                    "w": w,
                    "boundary_conditions": bcs,
                    "Type": 0,
                    "mesh": mesh,
                    "mode": 0
                    }

    solver= NSolver(solverparams)

#-------------------------------------------------------------------------------
#           Time Loop
#-------------------------------------------------------------------------------
    # Initialize half-sarcomere class. Methods used to calculate cross-bridges at gauss points
    hs = half_sarcomere.half_sarcomere(hs_params,1)
    CaList,VoltList = [],[]
    total_step = 1000
    if fiber_state == "normal":
        Ca_path = "/home/fenics/shared/Myosim/source_code/parameter_file/fiber/normal-calcium.dat"
        Volt_path = "/home/fenics/shared/Myosim/source_code/parameter_file/fiber/normal-volt.dat"
    elif fiber_state == "ischemia-1a":
        Ca_path = "/home/fenics/shared/Myosim/source_code/parameter_file/fiber/ischemia-1a-calcium.dat"
        Volt_path = "/home/fenics/shared/Myosim/source_code/parameter_file/fiber/ischemia-1a-volt.dat"
    elif fiber_state == "ischemia-1b":
        Ca_path = "/home/fenics/shared/Myosim/source_code/parameter_file/fiber/ischemia-1b-calcium.dat"
        Volt_path = "/home/fenics/shared/Myosim/source_code/parameter_file/fiber/ischemia-1b-volt.dat"
    elif fiber_state == "MI-short":
        Ca_path = "/home/fenics/shared/Myosim/source_code/parameter_file/fiber/MI-short-calcium.dat"
        Volt_path = "/home/fenics/shared/Myosim/source_code/parameter_file/fiber/MI-short-volt.dat"
    elif fiber_state == "MI-long":
        Ca_path = "/home/fenics/shared/Myosim/source_code/parameter_file/fiber/MI-long-calcium.dat"
        Volt_path = "/home/fenics/shared/Myosim/source_code/parameter_file/fiber/MI-long-volt.dat"
    elif fiber_state == "disarray_0.25":
        Ca_path = "/home/fenics/shared/Myosim/source_code/parameter_file/fiber/calcium_disarray_0.25.dat"
        Volt_path = "/home/fenics/shared/Myosim/source_code/parameter_file/fiber/volt_disarray_0.25.dat"
    elif fiber_state == "disarray_1.0":
        Ca_path = "/home/fenics/shared/Myosim/source_code/parameter_file/fiber/calcium_disarray_1.0.dat"
        Volt_path = "/home/fenics/shared/Myosim/source_code/parameter_file/fiber/volt_disarray_1.0.dat"
    elif fiber_state == "disarray_1.5":
        Ca_path = "/home/fenics/shared/Myosim/source_code/parameter_file/fiber/calcium_disarray_1.5.dat"
        Volt_path = "/home/fenics/shared/Myosim/source_code/parameter_file/fiber/volt_disarray_1.5.dat"
    elif fiber_state == "mix_single":
        Ca_path = "/home/fenics/shared/Myosim/source_code/parameter_file/fiber/SingleStim/calcium.dat"
        Volt_path = "/home/fenics/shared/Myosim/source_code/parameter_file/fiber/SingleStim/volt.dat"
    elif fiber_state == "mix_three":
        Ca_path = "/home/fenics/shared/Myosim/source_code/parameter_file/fiber/ThreeStim/calcium.dat"
        Volt_path = "/home/fenics/shared/Myosim/source_code/parameter_file/fiber/ThreeStim/volt.dat"
        total_step = 2000


    with open(Ca_path, "r") as f:
        for i in range(total_step):
            item = f.readline().split('\t')[1:-1]
            CaList.append(list(map(float, item)))

    with open(Volt_path, "r") as f:
        for i in range(total_step):
            item = f.readline().split('\t')[1:-1]
            VoltList.append(list(map(float, item)))

    for l in np.arange(no_of_time_steps):
        tic = timeit.default_timer()
        print "Time step number " + str(l)

        # At each gauss point, solve for cross-bridge distributions using myosim
        print "calling myosim"
        for mm in np.arange(no_of_int_points):
            temp_overlap[mm], y_interp[mm * n_array_length:(mm + 1) * n_array_length], y_vec_array_new[mm * n_array_length:(mm + 1) * n_array_length] = implement.update_simulation(hs, sim_timestep, delta_hsl_array[mm], hsl_array[mm], y_vec_array[mm*n_array_length:(mm+1)*n_array_length], p_f_array[mm], cb_f_array[mm], CaList[l][int(mm/4)], n_array_length, t,hs_params_list[mm])

        if save_cell_output:
            for  i in range(no_of_int_points):
                for j in range(n_array_length):
                    dumped_populations[i, j] = y_interp[i * n_array_length + j]

        # Update the populations
        y_vec_array = y_vec_array_new # for Myosim
        # Update the population function for fenics
        y_vec.vector()[:] = y_vec_array # for PDE
        # Update the array for myosim
        hsl_array_old = hsl_array
        # Update the hsl_old function for fenics
        hsl_old.vector()[:] = hsl_array_old[:]

        print "calling Newton Solver"
        # solve for displacement to satisfy balance of linear momentum
	try:
            solve(Ftotal == 0, w, bcs, J = Jac, form_compiler_parameters={"representation":"uflacs"},solver_parameters={"newton_solver":{"relative_tolerance":1e-8},"newton_solver":{"maximum_iterations":50},"newton_solver":{"absolute_tolerance":1e-8}})
        except:
            np.save(output_path + 'f0_vs_time.npy',f0_vs_time_array)
	    f0_dot_x_vec = inner(f0,x_vec)
	    f0_dot_x_vec_array = project(f0_dot_x_vec,Quad).vector().get_local()[:]
	    angles_array = np.arccos(f0_dot_x_vec_array)
	    np.save(output_path+"final_angles_array.npy",angles_array)

        print "PK2"
        PK2 = project(PK2_passive,TensorFunctionSpace(mesh,"DG",1),form_compiler_parameters={"representation":"uflacs"})
        # Update functions and arrays
        print "cb f array"
        cb_f_array[:] = project(cb_force, Quad).vector().get_local()[:]

        V0 = FunctionSpace(mesh, 'DG', 0)
        v0 = Function(V0)
        V_array = project(v0, V0).vector().get_local()[:]

        for i in range(len(V_array)):
            V_array[i] = VoltList[l][i]
        v0.vector().set_local(V_array)
        as_backend_type(v0.vector()).update_ghost_values()

        if save_visual_output:
            pk2temp = project(inner(f0,Pactive*f0),FunctionSpace(mesh,'DG',0),form_compiler_parameters={"representation":"uflacs"})
            pk2temp.rename("pk2_active","active_stress")
            v0.rename("volt", "Volt")
            output_file.write(v0, t[l])
            output_file.write(pk2temp,t[l])

        pseudo_old.vector()[:] = project(pseudo_alpha, Quad).vector().get_local()[:]
        hsl_array = project(hsl, Quad).vector().get_local()[:]           # for Myosim
        delta_hsl_array = project(sqrt(dot(f0, Cmat*f0))*hsl0, Quad).vector().get_local()[:] - hsl_array_old # for Myosim

        print "p_f_array"
        temp_DG = project(Sff, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        p_f = interpolate(temp_DG, Quad)
        p_f_array = p_f.vector().get_local()[:]

        for ii in range(np.shape(hsl_array)[0]):
            if p_f_array[ii] < 0.0:
                p_f_array[ii] = 0.0

        print "updating boundary conditions"
        bc_update_dict = update_boundary_conditions.update_bcs(bcs,sim_geometry,Ftotal,geo_options,sim_protocol,expressions,t[l],traction_switch_flag,x_dofs,test_marker_fcn,w,mesh,bcright,x_dir,l,W,facetboundaries,custom_disp_array)
        bcs = bc_update_dict["bcs"]
        if not (sim_geometry == "ventricle" or sim_geometry == "ellipsoid"):
            traction_switch_flag = bc_update_dict["traction_switch_flag"]
            rxn_force[l] = bc_update_dict["rxn_force"]
            u_D = bc_update_dict["expr"]["u_D"]
            Press = bc_update_dict["expr"]["Press"]
            print "current traction: ", Press.P


        # Save visualization info
        if save_visual_output:
            output_file.write(w.sub(0),t[l])
            hsl_temp = project(hsl, FunctionSpace(mesh, 'DG', 0))
            hsl_temp.rename("hsl_temp", "half-sarcomere length")
            output_file.write(hsl_temp, t[l])
            # Save fiber vectors associated with non-fibrotic regions separately
            f0_temp = project(f0, VectorFunctionSpace(mesh, "DG", 0))
            f0_temp.rename('f0', 'f0')
            fiber_file.write(f0_temp, t[l])
            s0_temp = project(s0, VectorFunctionSpace(mesh, "DG", 0))
            s0_temp.rename('s0', 's0')
            fiber_file.write(s0_temp, t[l])
            n0_temp = project(n0, VectorFunctionSpace(mesh, "DG", 0))
            n0_temp.rename('n0', 'n0')
            fiber_file.write(n0_temp, t[l])


        # Save cell info
        tic_save_cell = timeit.default_timer()
        if save_cell_output:
            active_stress_ds.iloc[0,:] = cb_f_array[:]
            active_stress_ds.to_csv(output_path + 'active_stress.csv',mode='a',header=False)

            #active_stress_ds = active_stress_ds.transpose()
            hsl_array_ds.iloc[0,:] = hsl_array[:]
            hsl_array_ds.to_csv(output_path + 'half_sarcomere_lengths.csv',mode='a',header=False)

            # calcium_ds.iloc[0,:] = calcium[l]
            # calcium_ds.to_csv(output_path + 'calcium.csv',mode='a',header=False)

            # for i in range(no_of_int_points):
            #     dumped_populations_ds.iloc[i,:] = dumped_populations[i,:]
            # dumped_populations_ds.to_csv(output_path + 'populations.csv',mode='a',header=False)

            #tarray_ds[l] = tarray[l]
            #tarray_ds.to_csv(output_path + 'time.csv',mode='a',header=False)
            # np.save(output_path+"time", t) # fix this

            p_f_array_ds.iloc[0,:] = p_f_array[:]
            p_f_array_ds.to_csv(output_path + 'myofiber_passive.csv',mode='a',header=False)

        if l == no_of_time_steps-1:
            for i in range(no_of_int_points):
                dumped_populations_ds.iloc[i,:] = dumped_populations[i,:]
            dumped_populations_ds.to_csv(output_path + 'populations.csv',mode='a',header=False)


        toc_save_cell = timeit.default_timer() - tic_save_cell
        print "time to save cell info = " + str(toc_save_cell)
        toc = timeit.default_timer() - tic
        print "time loop performance time = " + str(toc)

#-------------------------------------------------------------------------------
# for stand-alone testing
input_file = sys.argv[1]
#start_time = datetime.datetime.now()
start = timeit.default_timer()
# Load in JSON dictionary
with open(input_file, 'r') as json_input:
  input_parameters = json.load(json_input)

# Convert any unicode values to python strings so they work with some cpp libraries.
recode_dictionary.recode(input_parameters)

# Parse out the different types of parameters.
sim_params = input_parameters["simulation_parameters"]
passive_params = input_parameters["forms_parameters"]["passive_law_parameters"]
hs_params = input_parameters["myosim_parameters"]
cell_ion_params = input_parameters["electrophys_parameters"]["cell_ion_parameters"]
all_params = [sim_params,passive_params,hs_params,cell_ion_params]
#monodomain_params = input_parameters["electrophys_parameters"]["monodomain_parameters"]
if "windkessel_parameters" in input_parameters.keys():
    windkessel_params = input_parameters["windkessel_parameters"]
if "growth_and_remodeling" in input_parameters.keys():
    print "assigning growth_params"
    growth_params = input_parameters["growth_and_remodeling"]
    if 'growth_params' in locals():
        print "growth check worked"
    all_params.append(growth_params)
#optimization_params = input_parameters["optimization_parameters"]

fenics(sim_params)
sim_duration = timeit.default_timer() - start
#aplog.append_to_log(all_params,start_time,sim_duration,input_file)

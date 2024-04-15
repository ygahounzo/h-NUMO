OBJDIR = $(CURDIR)/objdir

include $(DEPEND_FILE)
include $(CONFIG_USER)

DIRS += 
.PHONY: $(DIRS)
.NOTPARALLEL:

SOURCEF += amain.F90
SOURCEF += apply_boundary_conditions.F90
SOURCEF += communicate_stencils.F90
SOURCEF += compute_element_size.F90
# SOURCEF += compute_norms.F90
SOURCEF += courant.F90
SOURCEF += create_boundary.F90
SOURCEF += create_face.F90
SOURCEF += create_global_rhs.F90
SOURCEF += create_grid_cube.F90
SOURCEF += create_grid_sphere_hex.F90
SOURCEF += create_grid_sphere_ico.F90
# SOURCEF += create_jump_penalty.F90
SOURCEF += create_local_global.F90
SOURCEF += create_mass.F90
SOURCEF += create_nbh_send_recv_local.F90
SOURCEF += create_normals.F90
SOURCEF += create_rhs_communicator.F90
SOURCEF += create_rhs_dynamics_flux.F90
SOURCEF += domain_decomp_geometric.F90
SOURCEF += domain_decomp_metis.F90
SOURCEF += filter.F90
SOURCEF += filter_init.F90
SOURCEF += gather_data.F90
SOURCEF += ico_quad.F90
SOURCEF += metrics.F90
SOURCEF += mod_adjncy.F90
SOURCEF += mod_array_util.F90
SOURCEF += mod_basis.F90
SOURCEF += mod_bc.F90
SOURCEF += mod_constants.F90
# SOURCEF += mod_convert_variables.F90
SOURCEF += mod_face.F90
SOURCEF += mod_filter.F90
SOURCEF += mod_global_grid.F90
SOURCEF += mod_grid.F90
SOURCEF += mod_initial.F90
SOURCEF += mod_input.F90
SOURCEF += mod_interface.F90
SOURCEF += mod_legendre.F90
# SOURCEF += mod_limiter.F90
SOURCEF += mod_metrics.F90
SOURCEF += mod_mpi_communicator.F90
# SOURCEF += mod_norms.F90
SOURCEF += mod_papi.F90
SOURCEF += mod_parallel.F90
SOURCEF += mod_p4est.F90
SOURCEF += mod_ref.F90
SOURCEF += mod_time_loop_mlswe.F90
SOURCEF += mod_viscosity.F90
SOURCEF += mod_types.F90
SOURCEF += mod_utilities.F90
SOURCEF += mod_mpi_utilities.F90
SOURCEF += mod_vtk_binary.F90
SOURCEF += mxm.F90
SOURCEF += nseam_out.F90
SOURCEF += outgraph.F90
SOURCEF += outmesh_2d_g.F90
SOURCEF += outmesh_level_g.F90
SOURCEF += outvtk_g_binary.F90
SOURCEF += outvtk_xml.F90
SOURCEF += print_diagnostics.F90
SOURCEF += print_header.F90
SOURCEF += qsort.F90
SOURCEF += read_gmsh.F90
SOURCEF += send_receive_bound.F90
SOURCEF += write_grid.F90
SOURCEF += write_output.F90

SOURCEC :=
SOURCEC += p8est.c
SOURCEC += p6est.c
SOURCEC += p4est.c

SOURCEH :=
#INCLUDE FILES FOR p8est.c
SOURCEH += p4est.c

#INCLUDE FILES FOR MLSWE
SOURCEF += mod_Tensorproduct.F90
SOURCEF += mod_initial_mlswe.F90
SOURCEF += mod_layer_terms.F90
SOURCEF += mod_barotropic_terms.F90
SOURCEF += metrics_quad.F90
SOURCEF += create_normals_quad.F90
SOURCEF += initial_conditions_mlswe.F90
SOURCEF += mod_create_rhs_mlswe.F90
SOURCEF += mod_splitting.F90
SOURCEF += ti_mlswe.F90
SOURCEF += diagnostics.F90
SOURCEF += mod_laplacian_quad.F90
SOURCEF += mod_rk_mlswe.F90
SOURCEF += ti_rk35_mlswe.F90


OBJS += $(addsuffix .o, $(basename $(SOURCEF) $(SOURCEC)))

ifdef OBJDIR
  OBJS := $(addprefix $(OBJDIR)/,$(OBJS))
endif

all: $(TARGET)

$(OBJDIR):
	@mkdir -p $(OBJDIR)

$(TARGET) : $(DIRS) $(OBJDIR) $(OBJS)
	$(MAKE) modmv
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $@

$(DIRS) :
	$(MAKE) -C $@ $(TARGET)

emptyrule:

clean :
	rm -rf $(OBJDIR)/*.o *.mod 
	for dir in $(DIRS); do $(MAKE) -C $$dir $@; done

deepclean :
	for dir in $(DIRS); do $(MAKE) -C $$dir $@; done
	rm -rf $(OBJDIR)/*.o *.mod 
ifdef OBJDIR
	rm -rf $(OBJDIR)
endif

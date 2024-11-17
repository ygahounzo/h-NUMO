/*
   p4est is a C library to manage a collection (a forest) of multiple
   connected adaptive quadtrees or octrees in parallel.

   Copyright (C) 2010 The University of Texas System
   Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

   p4est is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   p4est is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with p4est; if not, write to the Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

// #define VTK_OUTPUT

#include <sc.h>

#include <stdio.h>
#include <stdlib.h>

#include <p4est_bits.h>
#include <p4est_connectivity.h>
#include <p4est_extended.h>
#include <p4est_mesh.h>
#include <p6est.h>
#include <p6est_extended.h>
#include <p6est_ghost.h>
#include <p6est_lnodes.h>
#include <p6est_vtk.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))

// #include <p4est_connrefine.h>  // Yao Gahounzo (06/24/2024) comment out to use with new p4est

// define single/double precision
#ifdef SINGLE
#define real float
#else
#define real double
#endif

typedef struct
{
  int col_iref; // Marker for tracking refinement
  int col_oref; // original marker needed for flagging neighbors
                // (so corner connection do not propagate more than one element)
} user_data_t;

typedef struct
{
  sc_MPI_Comm mpicomm;
  int mpisize;
  int mpirank;
} mpi_context_t;

typedef struct
{
  int nop;  // nopx,nopy
  int nopz; // use different polynomial order in the vertical
  p4est_locidx_t npoin;
  p4est_locidx_t nelem;
  p4est_locidx_t npoin_dg;
  int num_nbh;
  int nz;    // number of points in a column
  int ncol;  // number of point-columns
  int nelxy; // number of elements in horizontal
  int nelz;  // number of elements in a column
  real *coord_dg;
  real *coord_cg;
  p4est_locidx_t *intma_table;
  p4est_locidx_t *D2C_mask;
  int num_send_recv_total;
  int *num_send_recv;
  int *nbh_proc;
  int *nbh_send_recv;
  p4est_locidx_t *nbh_send_recv_multi;
  p4est_locidx_t *nbh_send_recv_half;

  int *EToNC;
  int nNC;
  int *NC_face;
  int *NC_edge;
  int FACE_LEN;

  p4est_locidx_t nbsido;
  p4est_locidx_t nface;
  p4est_locidx_t nboun;
  p4est_locidx_t *bsido;
  p4est_locidx_t *face;
  p4est_locidx_t *face_type;
  p4est_locidx_t *node_column;
} p6esttonuma_t;

typedef struct
{
  int proc;   // sort key 1
  int ghostk; // sort key 2
  int nf;     // sort key 3
  int face;
} face_sort_t;

static int my_rank;
static const p4est_connect_type_t CONNECT_TYPE = P4EST_CONNECT_FULL;
static const p8est_connect_type_t p8est_connect_type = P8EST_CONNECT_FULL;
#ifdef VTK_OUTPUT
int p4est_output_count = 0;
#endif

/* This function was added here by MAK to create cubed-sphere topology
 */

static p4est_connectivity_t *p4est_connectivity_new_cubed_sphere(void)
{
  const p4est_topidx_t num_vertices = 8;
  const p4est_topidx_t num_trees = 6;
  const p4est_topidx_t num_ctt = 0;
  const double vertices[8 * 3] = {
      -1, -1, -1, 1, -1, -1, -1, 1, -1, 1, 1, -1,
      -1, -1, 1,  1, -1, 1,  -1, 1, 1,  1, 1, 1,
  };
  const p4est_topidx_t tree_to_vertex[6 * 4] = {
      0, 2, 1, 3, 2, 6, 3, 7, 0, 4, 2, 6, 4, 5, 6, 7, 0, 1, 4, 5, 1, 3, 5, 7,
  };
  const p4est_topidx_t tree_to_tree[6 * 4] = {
      4, 1, 2, 5, 0, 3, 2, 5, 0, 3, 4, 1, 2, 5, 4, 1, 2, 5, 0, 3, 4, 1, 0, 3,
  };
  const int8_t tree_to_face[6 * 4] = {
      2, 0, 0, 2, 1, 3, 3, 1, 2, 0, 0, 2, 1, 3, 3, 1, 2, 0, 0, 2, 1, 3, 3, 1,
  };

  return p4est_connectivity_new_copy(num_vertices, num_trees, 0, vertices,
                                     tree_to_vertex, tree_to_tree, tree_to_face,
                                     NULL, &num_ctt, NULL, NULL);
}

static void fill_coordinates_intma_p6est(int Nrp, int Nrpv, real *r, real *rz,
                                         p6est_t *p6est, p6est_lnodes_t *lnodes,
                                         p6esttonuma_t *p2n)
{

  p6est_connectivity_t *connectivity = p6est->connectivity;
  p4est_t *p4est = p6est->columns;
  sc_array_t *layers = p6est->layers;
  sc_array_t *trees = p4est->trees;
  const real intsize = 1.0 / P4EST_ROOT_LEN;
  double v[24];
  const p4est_topidx_t first_local_tree = p4est->first_local_tree;
  const p4est_topidx_t last_local_tree = p4est->last_local_tree;
  int xi, yi, j, k;
  int zi;
  int n, m, l, i, c;
  real h2, h2z, eta_x, eta_y, eta_z = 0.;
  real vxyz[P8EST_CHILDREN][3];
  size_t num_cols, zy, first, last;
  p4est_topidx_t jt;
  p4est_locidx_t quad_count, nelem_p6est;
  sc_array_t *columns;
  p4est_tree_t *tree;
  p4est_quadrant_t *col;
  p2est_quadrant_t *layer;

  int column_size = 0;

  real w[P8EST_CHILDREN];

  /* Allocate coordinates in DG storage */
  nelem_p6est = (int)lnodes->num_local_elements;
  p2n->coord_dg = (real *)malloc(sizeof(real) * 3 * p2n->npoin_dg);
  p2n->coord_cg = (real *)malloc(sizeof(real) * 3 * p2n->npoin);
  p2n->intma_table = malloc(sizeof(p4est_locidx_t) * p2n->npoin_dg);
  p2n->D2C_mask = malloc(sizeof(p4est_locidx_t) * p2n->npoin_dg);
  memset(p2n->D2C_mask, 0, sizeof(p4est_locidx_t) * p2n->npoin_dg);
  p4est_locidx_t *C2D_found =
      malloc(sizeof(p4est_locidx_t) * lnodes->owned_count);
  memset(C2D_found, 0, sizeof(p4est_locidx_t) * lnodes->owned_count);

  // loop over local trees
  for (jt = first_local_tree, quad_count = 0; jt <= last_local_tree; ++jt)
  {
    tree = p4est_tree_array_index(trees, jt);
    columns = &tree->quadrants;
    num_cols = columns->elem_count;
    p6est_tree_get_vertices(connectivity, jt, v);

    // loop over columns in each tree
    for (p4est_locidx_t q_v = 0; q_v < (p4est_locidx_t)num_cols; ++q_v)
    {
      col = p4est_quadrant_array_index(columns, q_v);
      P6EST_COLUMN_GET_RANGE(col, &first, &last);

      if (column_size == 0)
        column_size = (int)(last - first);
      else
        // Only allow AMR in the horizontal
        P4EST_ASSERT(column_size == (int)(last - first));
      P4EST_ASSERT(column_size > 0);

      // loop over elements (layers) in each column
      for (zy = first; zy < last; zy++, quad_count++)
      {
        layer = p2est_quadrant_array_index(layers, zy);
        h2 = .5 * intsize * P4EST_QUADRANT_LEN(col->level);
        h2z = .5 * intsize * P4EST_QUADRANT_LEN(layer->level);
        k = 0;

        // loop over element vertices
        for (zi = 0; zi < 2; ++zi)
        {
          for (yi = 0; yi < 2; ++yi)
          {
            for (xi = 0; xi < 2; ++xi)
            {
              P4EST_ASSERT(0 <= k && k < P8EST_CHILDREN);
              eta_x = intsize * col->x + h2 * (1. + (xi * 2 - 1));
              eta_y = intsize * col->y + h2 * (1. + (yi * 2 - 1));
              eta_z = intsize * layer->z + h2z * (1. + (zi * 2 - 1));

              // loop over x,y,z
              for (j = 0; j < 3; ++j)
              {
                vxyz[k][j] =
                    ((1. - eta_z) *
                         ((1. - eta_y) * ((1. - eta_x) * v[3 * 0 + j] +
                                          eta_x * v[3 * 1 + j]) +
                          eta_y * ((1. - eta_x) * v[3 * 2 + j] +
                                   eta_x * v[3 * 3 + j])) +
                     eta_z * ((1. - eta_y) * ((1. - eta_x) * v[3 * 4 + j] +
                                              eta_x * v[3 * 5 + j]) +
                              eta_y * ((1. - eta_x) * v[3 * 6 + j] +
                                       eta_x * v[3 * 7 + j])));
              } // loop over x,y,z
              ++k;
            } // xi
          }   // yi
        }     // loop over element vertices

        P4EST_ASSERT(k == P8EST_CHILDREN);

        /* Loop over element degrees of freedom to fill the coordinates */
        /* coord(3,1:npoin): x,y,z-coordinates of each point */

        for (m = 0; m < Nrp; ++m)
        {
          for (l = 0; l < Nrp; ++l)
          {
            for (n = 0; n < Nrpv; ++n)
            { // vertical coordinate varies fastest in p6est
              w[0] = (1 - r[l]) * (1 - r[m]) * (1 - rz[n]);
              w[1] = (1 + r[l]) * (1 - r[m]) * (1 - rz[n]);
              w[2] = (1 - r[l]) * (1 + r[m]) * (1 - rz[n]);
              w[3] = (1 + r[l]) * (1 + r[m]) * (1 - rz[n]);
              w[4] = (1 - r[l]) * (1 - r[m]) * (1 + rz[n]);
              w[5] = (1 + r[l]) * (1 - r[m]) * (1 + rz[n]);
              w[6] = (1 - r[l]) * (1 + r[m]) * (1 + rz[n]);
              w[7] = (1 + r[l]) * (1 + r[m]) * (1 + rz[n]);

              p4est_locidx_t offset =
                  m * Nrp * Nrpv + l * Nrpv + n; // different ordering in p6est
              p4est_locidx_t offset_intma = n * Nrp * Nrpv + m * Nrpv + l;
              const p4est_locidx_t dg_index =
                  Nrp * Nrp * Nrpv * quad_count + offset;
              const p4est_locidx_t cg_index = lnodes->element_nodes[dg_index];
              const p4est_locidx_t intma_index =
                  Nrp * Nrp * Nrpv * quad_count + offset_intma;
              real temp;

              for (i = 0; i < 3; ++i)
              {
                temp = 0;
                for (c = 0; c < P8EST_CHILDREN; c++)
                  temp += w[c] * vxyz[c][i];
                temp /= P8EST_CHILDREN;
                p2n->coord_dg[3 * intma_index + i] = temp;
                p2n->coord_cg[3 * cg_index + i] = temp;
              }

              /* intma_table(1:nglx,1:ngly,1:nglz,1:nelem): gives point number
                 (1...npoin) for each point in each element
                 (ngl*=nop*+1, nop: polynomial degree)*/
              p2n->intma_table[intma_index] = cg_index + 1;
              p2n->intma_table[intma_index] = cg_index + 1;

              /* If this is index is local, add it to the CG to DG map */
              if (cg_index < lnodes->owned_count && C2D_found[cg_index] == 0)
              {
                p2n->D2C_mask[intma_index] = 1;
                C2D_found[cg_index] = 1;
              }
              else
                p2n->D2C_mask[intma_index] = 0;
            }
          }
        } // loop over element dof

      } // loop over elements (layers) in each columns

    } // loop over columns in each tree

  } // loop over local trees

  p2n->nelz = column_size;
  // Works since we only allow horizontal AMR
  p2n->nelxy = column_size > 0 ? nelem_p6est / column_size : 0;
  P4EST_ASSERT(p2n->nelxy * column_size == nelem_p6est);

#ifdef P4EST_ENABLE_DEBUG
  /* If in debug mode check the mapping is complete */
  for (p4est_locidx_t n = 0; n < lnodes->owned_count; ++n)
    P4EST_ASSERT(C2D_found[n] > 0);
#endif
  free(C2D_found);
}

/* Function added by MAK from p6est/test/test_all.c */
void init_fn_p6est(p6est_t *p6est, p4est_topidx_t which_tree,
                   p4est_quadrant_t *col, p2est_quadrant_t *layer)
{
  user_data_t *data = (user_data_t *)layer->p.user_data;
  data->col_iref = 0;
  data->col_oref = 0;
}

void column_replace_fn(p6est_t *p6est, p4est_topidx_t which_tree,
                       int num_outcolumns, int num_outlayers,
                       p4est_quadrant_t *outcolumns[],
                       p2est_quadrant_t *outlayers[], int num_incolumns,
                       int num_inlayers, p4est_quadrant_t *incolumns[],
                       p2est_quadrant_t *inlayers[])
{
  if (num_outcolumns == 1)
  {
    SC_ASSERT(num_outcolumns == 1);
    SC_ASSERT(num_outlayers == 1);

    SC_ASSERT(num_incolumns == P4EST_CHILDREN);
    SC_ASSERT(num_inlayers == P4EST_CHILDREN);

    p2est_quadrant_t *out_layer = outlayers[0];
    user_data_t *out_data = (user_data_t *)out_layer->p.user_data;

    for (int i = 0; i < num_inlayers; ++i)
    {
      p2est_quadrant_t *in_layer = inlayers[i];
      user_data_t *in_data = (user_data_t *)in_layer->p.user_data;
      in_data->col_iref = out_data->col_iref;
    }
  }
  else
  {
    SC_ASSERT(num_outcolumns == P4EST_CHILDREN);
    SC_ASSERT(num_outlayers == 1);

    SC_ASSERT(num_incolumns == 1);
    SC_ASSERT(num_inlayers == 1);

    p2est_quadrant_t *out_layer = outlayers[0];
    user_data_t *out_data = (user_data_t *)out_layer->p.user_data;

    p2est_quadrant_t *in_layer = inlayers[0];
    user_data_t *in_data = (user_data_t *)in_layer->p.user_data;
    in_data->col_iref = out_data->col_iref;
#ifdef P4EST_ENABLE_DEBUG
    num_outlayers = 4; // do to bug in p4est noted above
    // In debug mode make sure we are uniform in horizontal
    for (int i = 0; i < num_outlayers; ++i)
    {
      p2est_quadrant_t *out_layer2 = outlayers[i];
      user_data_t *out_data2 = (user_data_t *)out_layer2->p.user_data;
      SC_ASSERT(out_data->col_iref = out_data2->col_iref);
    }
#endif
  }
}

/* If refine is column is marked for refinment */
static int coarsen_column_fn(p6est_t *p6est, p4est_topidx_t which_tree,
                             p4est_quadrant_t *columns[])
{

  for (p4est_locidx_t i = 0; i < P4EST_CHILDREN; i++)
  {
    p4est_quadrant_t *column = columns[i];
    size_t first, last;
    P6EST_COLUMN_GET_RANGE(column, &first, &last);
    p2est_quadrant_t *layer = p2est_quadrant_array_index(p6est->layers, first);
    user_data_t *data = (user_data_t *)layer->p.user_data;
#ifdef P4EST_ENABLE_DEBUG
    for (size_t q_v = first; q_v < last; ++q_v)
    {
      p2est_quadrant_t *layer2 = p2est_quadrant_array_index(p6est->layers, q_v);
      user_data_t *data2 = (user_data_t *)layer2->p.user_data;
      SC_ASSERT(data->col_iref == data2->col_iref);
      SC_ASSERT(data->col_oref == data2->col_oref);
    }
#endif
    if (data->col_iref >= 0)
      return 0;
  }
  return 1;
}

/* If refine is column is marked for refinment below current level */
static int refine_init_column_fn(p6est_t *p6est, p4est_topidx_t which_tree,
                                 p4est_quadrant_t *column)
{
  size_t first, last;
  P6EST_COLUMN_GET_RANGE(column, &first, &last);
  p2est_quadrant_t *layer = p2est_quadrant_array_index(p6est->layers, first);
  user_data_t *data = (user_data_t *)layer->p.user_data;
  int r = data->col_iref > column->level;
#ifdef P4EST_ENABLE_DEBUG
  // In debug mode make sure we are uniform in horizontal
  for (size_t q_v = first; q_v < last; ++q_v)
  {
    p2est_quadrant_t *layer2 = p2est_quadrant_array_index(p6est->layers, q_v);
    user_data_t *data2 = (user_data_t *)layer2->p.user_data;
    SC_ASSERT(data->col_iref == data2->col_iref);
  }
#endif
  return r;
}

static int refine_column_fn(p6est_t *p6est, p4est_topidx_t which_tree,
                            p4est_quadrant_t *column)
{
  size_t first, last;
  P6EST_COLUMN_GET_RANGE(column, &first, &last);
  p2est_quadrant_t *layer = p2est_quadrant_array_index(p6est->layers, first);
  user_data_t *data = (user_data_t *)layer->p.user_data;
  int r = data->col_iref > 0;
#ifdef P4EST_ENABLE_DEBUG
  // In debug mode make sure we are uniform in horizontal
  for (size_t q_v = first; q_v < last; ++q_v)
  {
    p2est_quadrant_t *layer2 = p2est_quadrant_array_index(p6est->layers, q_v);
    user_data_t *data2 = (user_data_t *)layer2->p.user_data;
    SC_ASSERT(data->col_iref == data2->col_iref);
  }
#endif
  return r;
}

/* Keep state of p4est */
static p6est_t *stored_p6est = NULL;
extern mpi_context_t mpi_context;

/* End of p6est section */

#ifdef ibmcompiler
void p6esttonuma_init
#else
void p6esttonuma_init_
#endif
    (int *nref_levs, int *nref_zlevs, int *is_non_conforming_flg,
     int *refine_columns, int *num_hroot, int *num_zroot, int *min_zlevel,
     int *is_cube, int *nnx, int *nny, int *nnz, int *xperiodic_flg,
     int *yperiodic_flg, int *iboundary, p6esttonuma_t **p2n)
{

  /* Here goes the p6est stuff */

  p6est_t *p6est;

  /* And here goes the p4est stuff */
  mpi_context_t *mpi = &mpi_context;

  p4est_locidx_t num_zroots = *num_zroot;

  /* Select p8est connectivity based on space_method*/
  p4est_connectivity_t *conn4 = NULL;

  /*
   * Here goes the p6est stuff
   */
  if (!stored_p6est)
  {
    double *topVertices = NULL;
    double height[3] = {0, 0, 1};
    int is_geometry_cube = *is_cube;
    p6est_connectivity_t *conn = NULL;
    if (!is_geometry_cube)
    {
      // create connectivity and forest structures
      p4est_connectivity_t *conn_init = p4est_connectivity_new_cubed_sphere();

      // conn4 = p4est_connectivity_new_cubed();
      conn4 = p4est_connectivity_refine(conn_init, *num_hroot);
      p4est_connectivity_destroy(conn_init);

      // create top vertices by scaling original vertices
      topVertices = malloc(sizeof(double) * conn4->num_vertices * 3);
      for (p4est_locidx_t i = 0; i < conn4->num_vertices * 3; ++i)
      {
        topVertices[i] = 2 * conn4->vertices[i];
      }
      // printf("%s","vertices created\n");
    }
    else
    {
      // here create flow-in-the-box mesh
      int nx = *nnx;
      int ny = *nny;
      int nz = *nnz;
      conn4 =
          p4est_connectivity_new_brick(nx, ny, *xperiodic_flg, *yperiodic_flg);
      num_zroots = nz; // redundant information in the input!!
    }

    conn = p6est_connectivity_new(conn4, topVertices, height);
    p4est_connectivity_destroy(conn4);
    p6est = p6est_new_ext(mpi->mpicomm, conn, 0, *nref_levs,
                          *min_zlevel + *nref_zlevs, num_zroots, 1,
                          sizeof(user_data_t), init_fn_p6est, NULL);
    if (topVertices)
    {
      free(topVertices);
      topVertices = NULL;
    }
  }
  else
    p6est = stored_p6est;

  /*
   * Here goes the refinement after initial mesh creation
   */
  if (*is_non_conforming_flg)
  {
    /* initial refinement data into the quadrant data */
    p4est_t *p4est = p6est->columns;
    p4est_ghost_t *ghost = p4est_ghost_new(p4est, CONNECT_TYPE);
    p4est_mesh_t *mesh = p4est_mesh_new_ext(p4est, ghost, 1, 0, CONNECT_TYPE);

    /* Mark the columns for refinement */
    for (p4est_locidx_t q_h = 0; q_h < p4est->local_num_quadrants; ++q_h)
    {
      /*  find which quadrant */
      p4est_quadrant_t *col = p4est_mesh_get_quadrant(p4est, mesh, q_h);
      size_t first, last;
      P6EST_COLUMN_GET_RANGE(col, &first, &last);
      for (size_t q_v = first; q_v < last; ++q_v)
      {
        p2est_quadrant_t *layer =
            p2est_quadrant_array_index(p6est->layers, q_v);
        user_data_t *data = (user_data_t *)layer->p.user_data;
        data->col_iref = refine_columns[q_h];
      }
    }

    p6est_refine_columns_ext(p6est, 1, -1, refine_init_column_fn, NULL,
                             column_replace_fn);

    p4est_mesh_destroy(mesh);
    mesh = NULL;
    p4est_ghost_destroy(ghost);
    ghost = NULL;
  }

  // balance after refinement
  p6est_balance(p6est, p8est_connect_type, NULL);

  // partition grid
  p6est_partition_ext(p6est, 1, NULL);

// Write vtk file
#ifdef VTK_OUTPUT
  char output[1024];
  sprintf(output, "p6est_mesh_%05d", p4est_output_count);
  p6est_vtk_write_file(p6est, output); /* Commented by FXG */
  ++p4est_output_count;
#endif

  stored_p6est = p6est;
}

#ifdef ibmcompiler
void p6esttonuma_fill_data
#else
void p6esttonuma_fill_data_
#endif
    (int *nnop, int *nnopz, real *xgl, real *xglz, int *is_cube, int *nis_dg,
     int *mFACE_LEN, int *iboundary, p6esttonuma_t **p2n)
{

  int is_geometry_cube = *is_cube;

  /* Here goes the p6est stuff */

  p6est_t *p6est = stored_p6est;
  SC_ASSERT(p6est);

  /* And here goes the p4est stuff */

  mpi_context_t *mpi = &mpi_context;

  int nop = *nnop;
  int nopz = *nnopz;
  int is_dg = *nis_dg;
  *p2n = malloc(sizeof(p6esttonuma_t));

  int FACE_LEN = (*p2n)->FACE_LEN = *mFACE_LEN;

  /* Create the ghost layer to learn about parallel neighbors. */
  p6est_ghost_t *ghost_p6est = p6est_ghost_new(p6est, CONNECT_TYPE);

  /* Create a node numbering for continuous linear finite elements. */
  p6est_lnodes_t *lnodes = p6est_lnodes_new(p6est, ghost_p6est, nop);

  p4est_t *p4est = p6est->columns;
  p4est_ghost_t *ghost_p4est = p4est_ghost_new(p4est, CONNECT_TYPE);

  p4est_lnodes_t *lnodes_p4est = p4est_lnodes_new(p4est, ghost_p4est, nop);
  p4est_lnodes_t *fnodes_p4est = p4est_lnodes_new(p4est, ghost_p4est, -1);
  p4est_lnodes_t *cnodes_p4est = is_dg ? fnodes_p4est : lnodes_p4est;

  p4est_mesh_t *mesh_p4est = p4est_mesh_new(p4est, ghost_p4est, CONNECT_TYPE);
  p4est_locidx_t nquads_h = mesh_p4est->local_num_quadrants;

  /* Get polynomial order */
  (*p2n)->nop = nop;
  (*p2n)->nopz = nopz;

  /* Get local (on rank) number of grid points */
  (*p2n)->npoin = lnodes->num_local_nodes;

  /* Get nelem: number of elements on the processor */
  int nelem_p6est = (*p2n)->nelem = (int)lnodes->num_local_elements;

  (*p2n)->npoin_dg = (*p2n)->nelem * (nop + 1) * (nop + 1) * (nopz + 1);

  fill_coordinates_intma_p6est(nop + 1, nopz + 1, xgl, xglz, p6est, lnodes,
                               *p2n);

  int *EToNC = (*p2n)->EToNC = malloc(sizeof(int) * nelem_p6est);

  int nNC = 0;
  for (int n = 0; n < nelem_p6est; ++n)
    if (lnodes->face_code[n])
    {
      EToNC[n] = nNC;
      ++nNC;
    }
    else
      EToNC[n] = -1;
  (*p2n)->nNC = nNC;

  /* For NC elements, determine how the faces/edges are hanging */
  int *NC_face = (*p2n)->NC_face = malloc(sizeof(int) * nNC * P8EST_FACES);
  int *NC_edge = (*p2n)->NC_edge = malloc(sizeof(int) * nNC * P8EST_EDGES);

  int faces_p6est[P8EST_FACES];
  int edges_p6est[P8EST_EDGES];

  int m = 0;
  for (int n = 0; n < nelem_p6est; ++n)
    if (lnodes->face_code[n])
    {
      p6est_lnodes_decode(lnodes->face_code[n], faces_p6est, edges_p6est);

      // Reorder faces and edges from p6est -> numa / p8est ordering
      NC_face[m * P8EST_FACES + 0] = faces_p6est[2];
      NC_face[m * P8EST_FACES + 1] = faces_p6est[3];
      NC_face[m * P8EST_FACES + 2] = faces_p6est[4];
      NC_face[m * P8EST_FACES + 3] = faces_p6est[5];
      NC_face[m * P8EST_FACES + 4] = faces_p6est[0];
      NC_face[m * P8EST_FACES + 5] = faces_p6est[1];

      NC_edge[m * P8EST_EDGES + 0] = edges_p6est[4];
      NC_edge[m * P8EST_EDGES + 1] = edges_p6est[6];
      NC_edge[m * P8EST_EDGES + 2] = edges_p6est[5];
      NC_edge[m * P8EST_EDGES + 3] = edges_p6est[7];

      NC_edge[m * P8EST_EDGES + 4] = edges_p6est[8];
      NC_edge[m * P8EST_EDGES + 5] = edges_p6est[10];
      NC_edge[m * P8EST_EDGES + 6] = edges_p6est[9];
      NC_edge[m * P8EST_EDGES + 7] = edges_p6est[11];

      NC_edge[m * P8EST_EDGES + 8] = edges_p6est[0];
      NC_edge[m * P8EST_EDGES + 9] = edges_p6est[1];
      NC_edge[m * P8EST_EDGES + 10] = edges_p6est[2];
      NC_edge[m * P8EST_EDGES + 11] = edges_p6est[3];

      // int h;
      // printf("\n");
      // for (h = 0; h < P8EST_FACES; ++h)
      //   printf("%d ", NC_face[m * P8EST_FACES + h]);
      // printf("\n");
      // for (h = 0; h < P8EST_EDGES; ++h)
      //   printf("%d ", NC_edge[m * P8EST_EDGES + h]);
      // printf("\n");

      ++m;
    }

  if (m != nNC)
    SC_ABORT("problem with counting face_code");

  // debug
  /*printf("%s","coord_cg\n");
    for(ll=0; ll<npoin_p6est; ++ll){
    printf("%d %f %f %f\n",ll,
    (*p2n)->coord[3*ll+0],(*p2n)->coord[3*ll+1],(*p2n)->coord[3*ll+2]);
    }*/

  /*
   *  Get column structure
   */

  (*p2n)->nz = (*p2n)->nelz * (*p2n)->nopz + 1;
  (*p2n)->ncol = (*p2n)->npoin / (*p2n)->nz;
  // Assumes AMR only in horizontal, e.g., all columns same size
  P4EST_ASSERT((*p2n)->ncol * (*p2n)->nz == (*p2n)->npoin);

  (*p2n)->node_column = malloc(sizeof(int) * (*p2n)->ncol * (*p2n)->nz);

  // loop over columns to fill node_column
  for (p4est_locidx_t zz = 0; zz < (*p2n)->npoin; ++zz)
    (*p2n)->node_column[zz] =
        zz + 1; // simplification due to inherent p6est node numbering scheme

  /*
   * Count boundary faces
   */
  p4est_locidx_t num_physical_boundary_faces = 0;
  p4est_locidx_t num_processor_boundary_faces = 0;

  /*
   * Each column has a top and bottom physical boundary
   */
  if (nopz > 0)
    num_physical_boundary_faces += 2 * (*p2n)->nelxy;

  /*
   * nbsido: number of (physical) domain boundary faces on the processor
   */

  // Count processor boundary faces
  for (p4est_locidx_t q_h = 0; q_h < nquads_h; ++q_h)
  {
    for (p4est_locidx_t f = 0; f < P4EST_FACES; ++f)
    {
      p4est_locidx_t nq_h = mesh_p4est->quad_to_quad[P4EST_FACES * q_h + f];
      p4est_locidx_t nv_h = mesh_p4est->quad_to_face[P4EST_FACES * q_h + f];
      p4est_locidx_t nf_h = nv_h % P4EST_FACES;

      if (nq_h == q_h && nf_h == f)
      {
        // Only the cube should have physical boundaries
        SC_ASSERT(is_geometry_cube);
        num_physical_boundary_faces += (*p2n)->nelz;
      }
      // If my neighbors face is negative I'm hanging
      else if (nv_h < 0)
      {
        p4est_locidx_t *nqs_h;
        int h;
        // Get neighbor half faces
        nqs_h = sc_array_index(mesh_p4est->quad_to_half, nq_h);
        // Check whether neighbor children are local or not
        for (h = 0; h < P4EST_HALF; ++h)
        {
          if (nqs_h[h] >= nquads_h)
            num_processor_boundary_faces += (*p2n)->nelz;
        }
      }
      // I'm conforming or a little guy, so see if my neighbor is local or not
      else if (nq_h >= nquads_h)
      {
        num_processor_boundary_faces += (*p2n)->nelz;
      }
    }
  }

  /*
   * nbsido: number of (physical) domain boundary faces on the processor
   */
  (*p2n)->nbsido = num_physical_boundary_faces;

  /*
   * nboun: processor boundary faces
   */
  (*p2n)->nboun = num_processor_boundary_faces;

  /*
   * nface: number of unique faces on the processor
   * (interior faces + processor boundary faces + domain boundary faces)
   *
   * Which since we only allow non-conforming in the horizontal is:
   */

  // horizontal faces
  (*p2n)->nface = (*p2n)->nelz * fnodes_p4est->num_local_nodes;

  // vertical faces
  if (nopz > 0)
    (*p2n)->nface += fnodes_p4est->num_local_elements * ((*p2n)->nelz + 1);

  my_rank = mpi->mpirank;

  /**********************
   * create faces
   ***********************/
  p4est_locidx_t ncolz = (*p2n)->nelz;

#define INTMA_H_FACE(q_h, q_v, f)                                              \
  ((fnodes_p4est->element_nodes[(q_h)*P4EST_FACES + (f)]) * ncolz + (q_v))

#define INTMA_V_FACE(q_h, q_v, f)                                              \
  (ncolz * fnodes_p4est->num_local_nodes + (ncolz + 1) * (q_h) + (q_v) +       \
   ((f)-P4EST_FACES))

#define INTMA_FACE(q_h, q_v, f)                                                \
  ((f < P4EST_FACES) ? INTMA_H_FACE(q_h, q_v, f) : INTMA_V_FACE(q_h, q_v, f))

  /*
   * bsido(1:6,1:nbsido): domain boundary data:
   * bsido(1:4,i): unused
   * bsido(5,i): element (1...nelem) to which face i belongs (only one because
   *             boundary face!)
   * bsido(6,i): boundary condition flag for face i (=4 for solid wall
   * boundary)
   */
#define BSIDO_LEN 6
#define BID_ELEM 4
#define BID_TYPE 5
  (*p2n)->bsido = malloc(sizeof(p4est_locidx_t) * BSIDO_LEN * (*p2n)->nbsido);
  memset((*p2n)->bsido, 0, sizeof(p4est_locidx_t) * BSIDO_LEN * (*p2n)->nbsido);

  /*
   * face(1:11, 1:nface): face data:
   * face(1:4,i): unused
   * face(5,i): local face relative to left element
   * face(6,i): local face relative to right element
   * face(7,i): element to left
   * face(8,i): element to the right or boundary condition flag (-bsido(6,sk))
   *
   * for non-conforming faces entries 8-11 store child numbers
   */
#define FID_LE_F 4
#define FID_RE_F 5
#define FID_LE_Q 6
#define FID_RE_Q 7
  (*p2n)->face = malloc(sizeof(p4est_locidx_t) * FACE_LEN * (*p2n)->nface);
  memset((*p2n)->face, 0, sizeof(p4est_locidx_t) * FACE_LEN * (*p2n)->nface);

  /*
   * face_type(i):
   * 1:  conforming internal face
   * 2:  parallel face
   * 4:  boundary face
   * 12: parllel NC face (parent is local)
   * 21: parllel NC face (parent is ghost)
   */
#define FT_LOCAL 1
#define FT_GHOST 2
#define FT_BOUND 4
#define FT_NC_LOCAL_PARENT 12
#define FT_NC_GHOST_PARENT 21
  (*p2n)->face_type = malloc(sizeof(p4est_locidx_t) * (*p2n)->nface);
  memset((*p2n)->face_type, 0, sizeof(p4est_locidx_t) * (*p2n)->nface);

  /*face orientation transformation array p4est->numa*/
  static const int numa_face[P8EST_FACES] = {4, 5, 2, 3, 0, 1};

#define HORZ_CONFORMING_LIMIT 8

  /*
   * parallel connectivity
   * nbh_proc: which ranks do I send to
   * num_send_recv: how many elements go to a rank
   */

  // Number of points in a column
  // for CG is number of element * polynomial order + 1
  // for DG is number of element
  p4est_locidx_t ncolpnts = is_dg ? ncolz : (ncolz * nopz + 1);

  p4est_locidx_t num_nbh;
  (*p2n)->num_nbh = num_nbh =
      MAX(0, (int)cnodes_p4est->sharers->elem_count - 1);
  (*p2n)->nbh_proc = malloc(sizeof(int) * num_nbh);
  (*p2n)->num_send_recv = malloc(sizeof(int) * num_nbh);
  (*p2n)->num_nbh = num_nbh;

  /*
   * Fill the sharing information for CG and DG
   */
  // Loop over neighbors
  p4est_locidx_t num_send_recv_total = 0;
  for (p4est_locidx_t iter = 0, num_nbh = 0;
       iter < (p4est_locidx_t)cnodes_p4est->sharers->elem_count; ++iter)
  {
    p4est_lnodes_rank_t *lrank =
        p4est_lnodes_rank_array_index(cnodes_p4est->sharers, iter);

    // If not myself: add the neighboring mpirank and count elements to comm
    if (lrank->rank != mpi->mpirank)
    {
      (*p2n)->nbh_proc[num_nbh] = lrank->rank + 1;
      (*p2n)->num_send_recv[num_nbh] =
          (int)lrank->shared_nodes.elem_count * ncolpnts;

      num_send_recv_total += (*p2n)->num_send_recv[num_nbh];

      num_nbh += 1;
    }
  }

  SC_ASSERT(num_nbh == (*p2n)->num_nbh);
  (*p2n)->num_send_recv_total = num_send_recv_total;

  /*
   * nbh_send_recv
   * for CG this is the the nodes to be sent / recv in global order
   * for DG this is the array of faces to send / recv in global order
   * first num_send_recv[0] elements go to rank nbh_proc[0]
   * next num_send_recv[1] elements go to rank nbh_proc[1]
   * etc.
   */
  (*p2n)->nbh_send_recv = malloc(sizeof(p4est_locidx_t) * num_send_recv_total);

  // Loop over the neighbors add elements / nodes to communication array
  p4est_locidx_t num_send_recv_offset = 0;
  for (p4est_locidx_t iter = 0;
       iter < (p4est_locidx_t)cnodes_p4est->sharers->elem_count; ++iter)
  {
    p4est_lnodes_rank_t *lrank =
        p4est_lnodes_rank_array_index(cnodes_p4est->sharers, iter);

    // if not myself, then add the nodes / elements to the list
    if (lrank->rank != mpi->mpirank)
    {

      // If DG loop over the horizontal points / faces
      p4est_locidx_t h_ter;
      for (h_ter = 0; h_ter < (p4est_locidx_t)lrank->shared_nodes.elem_count;
           ++h_ter)
      {
        // horizontal point ID
        p4est_locidx_t id_h =
            *(p4est_locidx_t *)sc_array_index(&lrank->shared_nodes, h_ter);
        // Loop up the elements and add the vertical elements
        // This works for CG since the columns are continuous storage
        for (p4est_locidx_t q_v = 0; q_v < ncolpnts; ++q_v)
        {
          // combine point ID
          p4est_locidx_t id_g = id_h * ncolpnts + q_v;
          (*p2n)->nbh_send_recv[num_send_recv_offset] = id_g + 1;
          ++num_send_recv_offset;
        }
      }
    }
  }

  /*
   * nbh_send_recv_multi
   * For each face / point count its multiplicity (used for non-conforming DG
   * to handle hanging faces)
   */
  (*p2n)->nbh_send_recv_multi =
      malloc(sizeof(p4est_locidx_t) * num_send_recv_total);
  memset((*p2n)->nbh_send_recv_multi, 0,
         sizeof(p4est_locidx_t) * num_send_recv_total);

  /*
   * nbh_send_recv_half
   * For each face I receive mark whether it is a child face (1) or not (0)
   */
  (*p2n)->nbh_send_recv_half =
      malloc(sizeof(p4est_locidx_t) * num_send_recv_total * P4EST_HALF);
  memset((*p2n)->nbh_send_recv_half, 0,
         sizeof(p4est_locidx_t) * num_send_recv_total * P4EST_HALF);

  /*
   * Build boundary and face maps
   */
  for (p4est_locidx_t q_h = 0, bc_iter = 0; q_h < nquads_h; ++q_h)
  {
    for (p4est_locidx_t f = 0; f < P4EST_FACES; ++f)
    {
      p4est_locidx_t f_numa = numa_face[f];

      p4est_locidx_t nq_h = mesh_p4est->quad_to_quad[P4EST_FACES * q_h + f];
      p4est_locidx_t nv_h = mesh_p4est->quad_to_face[P4EST_FACES * q_h + f];
      p4est_locidx_t nf_h = nv_h % P4EST_FACES;

      /*
       * Set the horizontal face maps
       */
      // Physical boundary for p4est -> horizontal boundary
      if (nq_h == q_h && nf_h == f)
      {
        // Only the cube should have physical boundaries
        SC_ASSERT(is_geometry_cube);
        for (p4est_locidx_t q_v = 0; q_v < ncolz; ++q_v)
        {
          // Global element number
          p4est_locidx_t q_g = q_h * ncolz + q_v;

          p4est_locidx_t fid = INTMA_H_FACE(q_h, q_v, f);
          (*p2n)->face_type[fid] = FT_BOUND;

          (*p2n)->face[FACE_LEN * fid + FID_LE_F] = f_numa + 1;
          (*p2n)->face[FACE_LEN * fid + FID_RE_F] = 0;
          (*p2n)->face[FACE_LEN * fid + FID_LE_Q] = q_g + 1;
          (*p2n)->face[FACE_LEN * fid + FID_RE_Q] = -iboundary[f_numa];

          (*p2n)->bsido[BSIDO_LEN * bc_iter + BID_ELEM] = q_g + 1;
          (*p2n)->bsido[BSIDO_LEN * bc_iter + BID_TYPE] = iboundary[f_numa];
          ++bc_iter;
        }
      }
      // Conforming face
      else if ((nv_h >= 0) && (nv_h < HORZ_CONFORMING_LIMIT))
      {
        // Neighbor is ghost
        if (nq_h >= nquads_h)
        {
          for (p4est_locidx_t q_v = 0; q_v < ncolz; ++q_v)
          {
            p4est_locidx_t q_g = q_h * ncolz + q_v;
            p4est_locidx_t sk = INTMA_H_FACE(q_h, q_v, f);
            (*p2n)->face_type[sk] = FT_GHOST;

            (*p2n)->face[FACE_LEN * sk + FID_LE_F] = f_numa + 1;
            (*p2n)->face[FACE_LEN * sk + FID_RE_F] = 0;
            (*p2n)->face[FACE_LEN * sk + FID_LE_Q] = q_g + 1;
            (*p2n)->face[FACE_LEN * sk + FID_RE_Q] = 0;

            // If this is DG, we need to fill the face parallel comm data
            if (is_dg)
            {
              // Get neighboring quad mpirank
              int nbh = mesh_p4est->ghost_to_proc[nq_h - nquads_h];

              // Should not be self
              SC_ASSERT(nbh != mpi->mpirank);

              // To check that we have found the neighbor
              int neighbor_found = 0;

              // start of the mpirank we are considering
              p4est_locidx_t rank_start = 0;
              // go over all neighbor procs
              for (p4est_locidx_t i = 0; i < num_nbh && !neighbor_found; ++i)
              {
                if ((*p2n)->nbh_proc[i] == nbh + 1)
                {

                  p4est_locidx_t h_ter;
                  for (h_ter = 0; h_ter < (*p2n)->num_send_recv[i]; ++h_ter)
                  {
                    // Check if we have found the face
                    if ((*p2n)->nbh_send_recv[h_ter + rank_start] == sk + 1)
                    {
                      ++(*p2n)->nbh_send_recv_multi[h_ter + rank_start];
                      ++neighbor_found;
                      break;
                    }
                  }
                }

                // Add on offset to get to next neighbor
                rank_start += (*p2n)->num_send_recv[i];
              }

              // Make sure we actually found the element!
              SC_ASSERT(neighbor_found);
            }
          }
        }
        // I am the left side of the face
        else if (q_h < nq_h || (q_h == nq_h && f < nf_h))
        {
          for (p4est_locidx_t q_v = 0; q_v < ncolz; ++q_v)
          {
            p4est_locidx_t q_g = q_h * ncolz + q_v;
            p4est_locidx_t sk = INTMA_H_FACE(q_h, q_v, f);
            (*p2n)->face_type[sk] = FT_LOCAL;

            (*p2n)->face[FACE_LEN * sk + FID_LE_F] = f_numa + 1;
            (*p2n)->face[FACE_LEN * sk + FID_LE_Q] = q_g + 1;
          }
        }
        // I am the right side of the face
        else
        {
          for (p4est_locidx_t q_v = 0; q_v < ncolz; ++q_v)
          {
            p4est_locidx_t q_g = q_h * ncolz + q_v;
            p4est_locidx_t sk = INTMA_H_FACE(q_h, q_v, f);
            (*p2n)->face_type[sk] = FT_LOCAL;

            (*p2n)->face[FACE_LEN * sk + FID_RE_F] = f_numa + 1;
            (*p2n)->face[FACE_LEN * sk + FID_RE_Q] = q_g + 1;
          }
        }
      }
      // non-conforming: I am the parent
      else if (nv_h < 0)
      {
        // Get child quadrant ids
        p4est_locidx_t *nqs_h = sc_array_index(mesh_p4est->quad_to_half, nq_h);

        // determine if we are a ghost or local face
        p4est_locidx_t nc_face_type = FT_LOCAL;
        p4est_locidx_t h;
        for (h = 0; h < P4EST_HALF; ++h)
          if (nqs_h[h] >= nquads_h)
            nc_face_type = FT_NC_LOCAL_PARENT;

        for (p4est_locidx_t q_v = 0; q_v < ncolz; ++q_v)
        {
          p4est_locidx_t q_g = q_h * ncolz + q_v;

          p4est_locidx_t sk = INTMA_H_FACE(q_h, q_v, f);

          // Parent (me) goes on the left
          (*p2n)->face[FACE_LEN * sk + FID_LE_F] = f_numa + 1;
          (*p2n)->face[FACE_LEN * sk + FID_LE_Q] = q_g + 1;

          (*p2n)->face_type[sk] = nc_face_type;

          // If this is DG, we need to fill the face parallel comm data
          if (is_dg)
          {

            // Loop through children and see if any are remote
            for (h = 0; h < P4EST_HALF; ++h)
            {
              if (nqs_h[h] >= mesh_p4est->local_num_quadrants)
              {
                // Get neighboring quad mpirank
                int nbh = mesh_p4est->ghost_to_proc[nqs_h[h] - nquads_h];

                // Should not be self
                SC_ASSERT(nbh != mpi->mpirank);

                // To check that we have found the neighbor
                int neighbor_found = 0;

                // start of the mpirank we are considering
                p4est_locidx_t rank_start = 0;

                // go over all neighbor procs
                for (p4est_locidx_t i = 0; i < num_nbh && !neighbor_found; ++i)
                {
                  if ((*p2n)->nbh_proc[i] == nbh + 1)
                  {
                    p4est_locidx_t h_ter;
                    for (h_ter = 0; h_ter < (*p2n)->num_send_recv[i]; ++h_ter)
                    {
                      p4est_locidx_t j = h_ter + rank_start;
                      // Check if we have found the face
                      if ((*p2n)->nbh_send_recv[j] == sk + 1)
                      {
                        ++(*p2n)->nbh_send_recv_multi[j];
                        ++neighbor_found;
                        (*p2n)->nbh_send_recv_half[j * P4EST_HALF + h] = 1;
                        break;
                      }
                    }
                  }

                  // Add on offset to get to next neighbor
                  rank_start += (*p2n)->num_send_recv[i];
                }

                // Make sure we actually found the element!
                SC_ASSERT(neighbor_found);
              }
            }
          }
        }
      }
      // non-conforming: I am the child
      else
      {
        // find my child id
        p4est_locidx_t h = nv_h / HORZ_CONFORMING_LIMIT - 1;

        for (p4est_locidx_t q_v = 0; q_v < ncolz; ++q_v)
        {
          p4est_locidx_t q_g = q_h * ncolz + q_v;

          p4est_locidx_t sk = INTMA_H_FACE(q_h, q_v, f);

          // Children go in the right element with shift for child id
          (*p2n)->face[FACE_LEN * sk + FID_RE_F] = f_numa + 1;
          (*p2n)->face[FACE_LEN * sk + FID_RE_Q + h] = q_g + 1;

          // If our neighbor is ghost, parent is ghost, thus mark
          // else parent is local, and they can mark the face correctly
          if (nq_h >= nquads_h)
          {
            (*p2n)->face_type[sk] = FT_NC_GHOST_PARENT;
            // If this is DG, we need to fill the face parallel comm data
            if (is_dg)
            {
              // Get neighboring quad mpirank
              int nbh = mesh_p4est->ghost_to_proc[nq_h - nquads_h];

              // Should not be self
              SC_ASSERT(nbh != mpi->mpirank);

              // To check that we have found the neighbor
              int neighbor_found = 0;

              // start of the mpirank we are considering
              p4est_locidx_t rank_start = 0;
              // go over all neighbor procs
              for (p4est_locidx_t i = 0; i < num_nbh && !neighbor_found; ++i)
              {
                if ((*p2n)->nbh_proc[i] == nbh + 1)
                {

                  p4est_locidx_t h_ter;
                  for (h_ter = 0; h_ter < (*p2n)->num_send_recv[i]; ++h_ter)
                  {
                    // Check if we have found the face
                    if ((*p2n)->nbh_send_recv[h_ter + rank_start] == sk + 1)
                    {
                      ++(*p2n)->nbh_send_recv_multi[h_ter + rank_start];
                      ++neighbor_found;
                      break;
                    }
                  }
                }

                // Add on offset to get to next neighbor
                rank_start += (*p2n)->num_send_recv[i];
              }

              // Make sure we actually found the element!
              SC_ASSERT(neighbor_found);
            }
          }
        }
      }
    }

    // top-bottom of cubes
    if (nopz > 0)
    {
      p4est_locidx_t f_bot = P4EST_FACES;
      p4est_locidx_t f_top = P4EST_FACES + 1;

      // bottom column face is boundary
      {
        p4est_locidx_t q_v = 0;
        p4est_locidx_t q_g = q_h * ncolz + q_v;
        p4est_locidx_t sk = INTMA_V_FACE(q_h, q_v, f_bot);
        (*p2n)->face_type[sk] = FT_BOUND;

        (*p2n)->face[FACE_LEN * sk + FID_LE_F] = numa_face[f_bot] + 1;
        (*p2n)->face[FACE_LEN * sk + FID_LE_Q] = q_g + 1;
        (*p2n)->face[FACE_LEN * sk + FID_RE_F] = 0;
        (*p2n)->face[FACE_LEN * sk + FID_RE_Q] = -iboundary[numa_face[f_bot]];

        (*p2n)->bsido[BSIDO_LEN * bc_iter + BID_ELEM] = q_g + 1;
        (*p2n)->bsido[BSIDO_LEN * bc_iter + BID_TYPE] =
            iboundary[numa_face[f_bot]];
        ++bc_iter;
      }

      // interior column faces
      for (p4est_locidx_t q_v = 0; q_v < ncolz - 1; ++q_v)
      {
        p4est_locidx_t q_g = q_h * ncolz + q_v;
        p4est_locidx_t sk = INTMA_V_FACE(q_h, q_v, f_top);
        (*p2n)->face_type[sk] = FT_LOCAL;

        (*p2n)->face[FACE_LEN * sk + FID_LE_F] = numa_face[f_top] + 1;
        (*p2n)->face[FACE_LEN * sk + FID_LE_Q] = q_g + 1;
        (*p2n)->face[FACE_LEN * sk + FID_RE_F] = numa_face[4] + 1;
        (*p2n)->face[FACE_LEN * sk + FID_RE_Q] = (q_g + 1) + 1;
      }

      // top column face is boundary
      {
        p4est_locidx_t q_v = ncolz - 1;
        p4est_locidx_t q_g = q_h * ncolz + q_v;
        p4est_locidx_t sk = INTMA_V_FACE(q_h, q_v, f_top);
        (*p2n)->face_type[sk] = FT_BOUND;

        (*p2n)->face[FACE_LEN * sk + FID_LE_F] = numa_face[f_top] + 1;
        (*p2n)->face[FACE_LEN * sk + FID_LE_Q] = q_g + 1;
        (*p2n)->face[FACE_LEN * sk + FID_RE_F] = 0;
        (*p2n)->face[FACE_LEN * sk + FID_RE_Q] = -iboundary[numa_face[f_top]];

        (*p2n)->bsido[BSIDO_LEN * bc_iter + BID_ELEM] = q_g + 1;
        (*p2n)->bsido[BSIDO_LEN * bc_iter + BID_TYPE] =
            iboundary[numa_face[f_top]];
        ++bc_iter;
      }
    }
  }

  /*
     printf("%s\n","--------------");
     printf("%s %d\n","nop",nop);
     printf("%s %d\n","npoin",npoin_p6est);
     printf("%s %d\n","nelem",nelem_p6est);
     printf("%s %d %d\n","ll",ll,nelem_p6est*(nop+1)*(nop+1)*(nop+1));
     printf("%s %d\n","num_nbh",num_nbh);
     printf("%s %d\n","num_send_recv_total",num_send_recv_total);
     printf("%s
     %d\n","num_physical_boundary_faces",num_physical_boundary_faces);
     printf("%s
     %d\n","num_processor_boundary_faces",num_processor_boundary_faces);
     printf("%s %d\n","nface",(*p2n)->nface);
     printf("%s\n","--------------");
   */
  p4est_mesh_destroy(mesh_p4est);
  p4est_ghost_destroy(ghost_p4est);
  p6est_lnodes_destroy(lnodes);
  p4est_lnodes_destroy(lnodes_p4est);
  p4est_lnodes_destroy(fnodes_p4est);
  p6est_ghost_destroy(ghost_p6est);

  /* End of p6est section */
}

#ifdef ibmcompiler
void p6esttonuma_get_mesh_scalars
#else
void p6esttonuma_get_mesh_scalars_
#endif
    (p6esttonuma_t **p2n, int *npoin, int *nelem, int *num_nbh,
     int *num_send_recv_total, int *nbsido, int *nface, int *nboun, int *nNC,
     int *ncol, int *nz, int *nelz)
{
  *npoin = (*p2n)->npoin;
  *nelem = (*p2n)->nelem;
  *num_nbh = (*p2n)->num_nbh;
  // num_send_recv_total: unique faces
  *num_send_recv_total = (*p2n)->num_send_recv_total;
  *nbsido = (*p2n)->nbsido;
  *nface = (*p2n)->nface;
  // nboun: duplicated parent faces
  *nboun = (*p2n)->nboun;
  *nNC = (*p2n)->nNC;
  *ncol = (*p2n)->ncol;
  *nz = (*p2n)->nz;
  *nelz = (*p2n)->nelz;
}

#ifdef ibmcompiler
void p6esttonuma_get_mesh_arrays
#else
void p6esttonuma_get_mesh_arrays_
#endif
    (p6esttonuma_t **p2n, real *coord, int *intma_table, int *NC_face,
     int *NC_edge, int *EToNC, int *face, int *face_type, int *num_send_recv,
     int *nbh_proc, int *nbh_send_recv, int *nbh_send_recv_multi,
     int *nbh_send_recv_half, int *bsido, int *node_column, int *nis_cgc,
     int *D2C_mask)
{
  int i;
  int is_cgc = *nis_cgc;

  if (is_cgc)
  {
    for (i = 0; i < 3 * (*p2n)->npoin; ++i)
      coord[i] = (*p2n)->coord_cg[i];
  }
  else
  {
    for (i = 0; i < 3 * (*p2n)->npoin_dg; ++i)
      coord[i] = (*p2n)->coord_dg[i];
  }

  for (i = 0; i < (*p2n)->npoin_dg; ++i)
    intma_table[i] = (*p2n)->intma_table[i];

  for (i = 0; i < (*p2n)->npoin_dg; ++i)
    D2C_mask[i] = (*p2n)->D2C_mask[i];

  for (i = 0; i < (*p2n)->num_nbh; ++i)
    num_send_recv[i] = (*p2n)->num_send_recv[i];

  for (i = 0; i < (*p2n)->num_nbh; ++i)
    nbh_proc[i] = (*p2n)->nbh_proc[i];

  for (i = 0; i < (*p2n)->num_send_recv_total; ++i)
    nbh_send_recv[i] = (*p2n)->nbh_send_recv[i];

  for (i = 0; i < (*p2n)->num_send_recv_total; ++i)
    nbh_send_recv_multi[i] = (*p2n)->nbh_send_recv_multi[i];

  for (i = 0; i < (*p2n)->num_send_recv_total * P4EST_HALF; ++i)
    nbh_send_recv_half[i] = (*p2n)->nbh_send_recv_half[i];

  for (i = 0; i < 6 * (*p2n)->nbsido; ++i)
    bsido[i] = (*p2n)->bsido[i];

  for (i = 0; i < (*p2n)->FACE_LEN * (*p2n)->nface; ++i)
    face[i] = (*p2n)->face[i];

  for (i = 0; i < (*p2n)->nface; ++i)
    face_type[i] = (*p2n)->face_type[i];

  for (i = 0; i < (*p2n)->nelem; ++i)
    EToNC[i] = (*p2n)->EToNC[i] + 1;

  for (i = 0; i < (*p2n)->nNC * P8EST_FACES; ++i)
    NC_face[i] = (*p2n)->NC_face[i];

  for (i = 0; i < (*p2n)->nNC * P8EST_EDGES; ++i)
    NC_edge[i] = (*p2n)->NC_edge[i];

  for (i = 0; i < ((*p2n)->nz) * ((*p2n)->ncol); ++i)
    node_column[i] = (*p2n)->node_column[i];
}

#ifdef ibmcompiler
void p6esttonuma_free
#else
void p6esttonuma_free_
#endif
    (p6esttonuma_t **p2n)
{
  free((*p2n)->coord_dg);
  free((*p2n)->coord_cg);
  free((*p2n)->intma_table);
  free((*p2n)->D2C_mask);
  free((*p2n)->num_send_recv);
  free((*p2n)->nbh_proc);
  free((*p2n)->nbh_send_recv);
  free((*p2n)->nbh_send_recv_multi);
  free((*p2n)->nbh_send_recv_half);
  free((*p2n)->bsido);
  free((*p2n)->face);
  free((*p2n)->face_type);
  free((*p2n)->node_column);
  free((*p2n)->EToNC);
  free((*p2n)->NC_face);
  free((*p2n)->NC_edge);
  free(*p2n);
}

/* START OF BFAM CODE */
/*
 * The following code is extracted from/based on the bfam project:
 *
 *    https://github.com/bfam/bfam
 *
 * Authors:
 *   Lucas C Wilcox
 *   Jeremy E Kozdon
 */

#ifdef ibmcompiler
void p6esttonuma_get_element_horiztonal_lvl
#else
void p6esttonuma_get_element_horiztonal_lvl_
#endif
    (int8_t *lvl, int32_t *N)
{
  p6est_t *p6est = stored_p6est;
  p4est_t *p4est = p6est->columns;

  p4est_ghost_t *ghost = p4est_ghost_new(p4est, CONNECT_TYPE);
  p4est_mesh_t *mesh = p4est_mesh_new_ext(p4est, ghost, 1, 0, CONNECT_TYPE);

  p4est_locidx_t nelem_v =
      (p4est_mesh_get_quadrant(p4est, mesh, 0))->p.piggy3.which_tree;

  p4est_locidx_t nelem_h = p4est->local_num_quadrants;

  SC_ASSERT(*N == nelem_v * nelem_h);

  /* Get the column mesh levels */
  for (p4est_locidx_t q_h = 0; q_h < nelem_h; ++q_h)
  {
    /*  find which quadrant */
    p4est_quadrant_t *col = p4est_mesh_get_quadrant(p4est, mesh, q_h);
    size_t first, last;
    P6EST_COLUMN_GET_RANGE(col, &first, &last);
    SC_ASSERT(last - first == (size_t)nelem_v);
    for (p4est_locidx_t q_v = 0; q_v < nelem_v; ++q_v)
      lvl[q_v + q_h * nelem_v] = col->level;
  }
  p4est_mesh_destroy(mesh);
  mesh = NULL;
  p4est_ghost_destroy(ghost);
  ghost = NULL;
}

#ifdef ibmcompiler
void p6esttonuma_mark_elements
#else
void p6esttonuma_mark_elements_
#endif
    (int *hadapt)
{
  p6est_t *p6est = stored_p6est;
  p4est_t *p4est = p6est->columns;

  p4est_ghost_t *ghost = p4est_ghost_new(p4est, CONNECT_TYPE);
  p4est_mesh_t *mesh = p4est_mesh_new_ext(p4est, ghost, 1, 0, CONNECT_TYPE);

  p4est_locidx_t nelem_v =
      (p4est_mesh_get_quadrant(p4est, mesh, 0))->p.piggy3.which_tree;

  p4est_locidx_t nelem_h = p4est->local_num_quadrants;

  /* Get the column mesh levels */
  for (p4est_locidx_t q_h = 0; q_h < nelem_h; ++q_h)
  {
    int col_flag = -1;
    for (p4est_locidx_t q_v = 0; q_v < nelem_v; ++q_v)
      col_flag = MAX(col_flag, hadapt[q_v + nelem_v * q_h]);

    /*  find which quadrant */
    p4est_quadrant_t *col = p4est_mesh_get_quadrant(p4est, mesh, q_h);
    size_t first, last;
    P6EST_COLUMN_GET_RANGE(col, &first, &last);
    SC_ASSERT(last - first == (size_t)nelem_v);
    for (size_t q_v = first; q_v < last; ++q_v)
    {
      p2est_quadrant_t *layer = p2est_quadrant_array_index(p6est->layers, q_v);
      user_data_t *data = (user_data_t *)layer->p.user_data;
      data->col_iref = col_flag;
    }
  }
  p4est_mesh_destroy(mesh);
  mesh = NULL;
  p4est_ghost_destroy(ghost);
  ghost = NULL;
}

static void mark_neighbors_volume_fun(p4est_iter_volume_info_t *info,
                                      void *user_data)
{
  p6est_t *p6est = stored_p6est;
  p4est_quadrant_t *column = info->quad;
  size_t first, last;
  P6EST_COLUMN_GET_RANGE(column, &first, &last);
  p2est_quadrant_t *layer = p2est_quadrant_array_index(p6est->layers, first);
  user_data_t *data = (user_data_t *)layer->p.user_data;
  data->col_oref = data->col_iref;
#ifdef P4EST_ENABLE_DEBUG
  // In debug mode make sure we are uniform in column
  for (size_t q_v = first; q_v < last; ++q_v)
  {
    p2est_quadrant_t *layer2 = p2est_quadrant_array_index(p6est->layers, q_v);
    user_data_t *data2 = (user_data_t *)layer2->p.user_data;
    data2->col_oref = data2->col_iref;
    SC_ASSERT(data->col_iref == data2->col_iref);
    SC_ASSERT(data->col_oref == data2->col_oref);
  }
#endif
}

static void mark_neighbors_corner_fun(p4est_iter_corner_info_t *info,
                                      void *user_data)
{
  p6est_t *p6est = stored_p6est;
  user_data_t *ghost_data = (user_data_t *)user_data;
  sc_array_t *sides = &(info->sides);
  int col_iref = 0;
  int col_oref = 0;
  int8_t lvl = 0;

  // loop over corner connections to determine whether to refine the corners
  for (p4est_locidx_t k = 0; k < (p4est_locidx_t)sides->elem_count; ++k)
  {
    p4est_iter_corner_side_t *side = p4est_iter_cside_array_index_int(sides, k);
    // If ghost in the forest used the ghost_data
    // otherwise read from the first quadrant layer
    if (side->is_ghost)
      col_oref = ghost_data[side->quadid].col_iref;
    else
    {
      p4est_quadrant_t *column = side->quad;
      size_t first, last;
      P6EST_COLUMN_GET_RANGE(column, &first, &last);
      p2est_quadrant_t *layer =
          p2est_quadrant_array_index(p6est->layers, first);
      user_data_t *data = (user_data_t *)layer->p.user_data;
      col_oref = data->col_oref;
    }

    if (col_oref > 0)
    {
      col_iref = 1;
      lvl = side->quad->level > lvl ? side->quad->level : lvl;
    }
  }

  // If one of the these elements should be refined, refine them all
  if (col_iref > 0)
    for (p4est_locidx_t k = 0; k < (p4est_locidx_t)sides->elem_count; ++k)
    {
      p4est_iter_corner_side_t *side =
          p4est_iter_cside_array_index_int(sides, k);
      if (side->is_ghost)
        continue;
      if (lvl >= side->quad->level)
      {
        p4est_quadrant_t *column = side->quad;
        size_t first, last;
        P6EST_COLUMN_GET_RANGE(column, &first, &last);
        p2est_quadrant_t *layer =
            p2est_quadrant_array_index(p6est->layers, first);
        user_data_t *data = (user_data_t *)layer->p.user_data;
        data->col_iref = 1;
#ifdef P4EST_ENABLE_DEBUG
        // In debug mode mark all elements in the column
        for (size_t q_v = first; q_v < last; ++q_v)
        {
          p2est_quadrant_t *layer2 =
              p2est_quadrant_array_index(p6est->layers, q_v);
          user_data_t *data2 = (user_data_t *)layer2->p.user_data;
          data2->col_iref = 1;
        }
#endif
      }
    }
}

#ifdef ibmcompiler
void p6esttonuma_mark_neighbors_p4est
#else
void p6esttonuma_mark_neighbors_p4est_
#endif
    (int32_t *num_iterations)
{
  p6est_t *p6est = stored_p6est;
  p4est_t *p4est = p6est->columns;

  p4est_ghost_t *ghost_p4est = p4est_ghost_new(p4est, CONNECT_TYPE);
  p4est_mesh_t *mesh_p4est =
      p4est_mesh_new_ext(p4est, ghost_p4est, 1, 0, CONNECT_TYPE);

  sc_array_t mirrors = ghost_p4est->mirrors;
  user_data_t *ghost_data = (user_data_t *)malloc(
      sizeof(user_data_t) * ghost_p4est->ghosts.elem_count);
  void **mirror_data =
      (void **)malloc(sizeof(user_data_t *) * mirrors.elem_count);

  // Fill mirror data
  for (p4est_locidx_t m = 0;
       m < (p4est_locidx_t)ghost_p4est->mirrors.elem_count; ++m)
  {
    // Get the mirror quadrant
    p4est_quadrant_t *mirror = p4est_quadrant_array_index(&mirrors, m);
    // Get the real quadrant from the p4est
    p4est_quadrant_t *column =
        p4est_mesh_get_quadrant(p4est, mesh_p4est, mirror->p.piggy3.local_num);

    size_t first, last;
    P6EST_COLUMN_GET_RANGE(column, &first, &last);
    p2est_quadrant_t *layer = p2est_quadrant_array_index(p6est->layers, first);
    mirror_data[m] = layer->p.user_data;
#ifdef P4EST_ENABLE_DEBUG
    // In debug mode make sure we are uniform in column
    user_data_t *data = (user_data_t *)layer->p.user_data;
    for (size_t q_v = first; q_v < last; ++q_v)
    {
      p2est_quadrant_t *layer2 = p2est_quadrant_array_index(p6est->layers, q_v);
      user_data_t *data2 = (user_data_t *)layer2->p.user_data;
      SC_ASSERT(data->col_iref == data2->col_iref);
      SC_ASSERT(data->col_oref == data2->col_oref);
    }
#endif
  }

  for (p4est_locidx_t k = 0; k < *num_iterations; ++k)
  {
    // exchange the ghost data
    p4est_ghost_exchange_custom(p4est, ghost_p4est, sizeof(user_data_t),
                                mirror_data, ghost_data);

    p4est_iterate(p4est, ghost_p4est, (void *)ghost_data,
                  mark_neighbors_volume_fun, NULL, mark_neighbors_corner_fun);
  }

  free(ghost_data);
  ghost_data = NULL;
  free(mirror_data);
  mirror_data = NULL;
  p4est_mesh_destroy(mesh_p4est);
  mesh_p4est = NULL;
  p4est_ghost_destroy(ghost_p4est);
  ghost_p4est = NULL;
}

#ifdef ibmcompiler
void p6esttonuma_coarsen_refine_p4est
#else
void p6esttonuma_coarsen_refine_p4est_
#endif
    (int32_t *num_dst)
{
  p6est_t *p6est = stored_p6est;

  p6est_coarsen_columns_ext(p6est, 0, 0, coarsen_column_fn, NULL,
                            column_replace_fn);
  p6est_refine_columns_ext(p6est, 0, -1, refine_column_fn, NULL,
                           column_replace_fn);

  p6est_balance_ext(p6est, p8est_connect_type, 0, 1, NULL, column_replace_fn);

#ifdef VTK_OUTPUT
  char output[1024];
  sprintf(output, "p6est_mesh_%05d", p4est_output_count);
  p6est_vtk_write_file(p6est, output); /* Commented by FXG */
  ++p4est_output_count;
#endif

  *num_dst = p6est->columns->local_num_quadrants;
}

#ifdef ibmcompiler
void p6esttonuma_repartition
#else
void p6esttonuma_repartition_
#endif
    (int64_t *qid_src, int64_t *qid_dst)
{
  p6est_t *p6est = stored_p6est;
  p4est_t *p4est = p6est->columns;

  for (int k = 0; k < p4est->mpisize + 1; ++k)
    qid_src[k] = p4est->global_first_quadrant[k];

  p6est_partition_ext(p6est, 1, NULL);

  for (int k = 0; k < p4est->mpisize + 1; ++k)
    qid_dst[k] = p4est->global_first_quadrant[k];

#ifdef VTK_OUTPUT
  char output[1024];
  sprintf(output, "p6est_mesh_%05d", p4est_output_count);
  p6est_vtk_write_file(p6est, output); /* Commented by FXG */
  ++p4est_output_count;
#endif
}
/* END OF BFAM CODE */

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

//#define VTK_OUTPUT

#include <sc.h>

#include <stdio.h>
#include <stdlib.h>

#include <p4est_bits.h>
#include <p4est_connectivity.h>
#include <p4est_mesh.h>
#include <p6est.h>
#include <p6est_extended.h>
#include <p6est_ghost.h>
#include <p6est_lnodes.h>
#include <p6est_vtk.h>

#include <p4est_connrefine.h>

// define single/double precision
#ifdef SINGLE
#define real float
#else
#define real double
#endif

typedef struct
{
  p4est_topidx_t a;
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
  p4est_locidx_t *intma;
  int num_send_recv_total;
  int *num_send_recv;
  int *nbh_proc;
  int *nbh_send_recv;
  int log_zroots;

  p4est_locidx_t nbsido;
  p4est_locidx_t nface;
  p4est_locidx_t nboun;
  p4est_locidx_t *bsido;
  p4est_locidx_t *face;
  p4est_locidx_t *face_type;
  p4est_locidx_t *node_column;
} p4esttonuma_t;

typedef struct
{
  int proc;   // sort key 1
  int ghostk; // sort key 2
  int nf;     // sort key 3
  int face;
} face_sort_t;

static int my_rank;

static int cmpfunc(const void *lhs, const void *rhs)
{
  int diff = (((face_sort_t *)lhs)->proc - ((face_sort_t *)rhs)->proc);
  if (diff == 0)
  {
    diff = ((face_sort_t *)lhs)->ghostk - ((face_sort_t *)rhs)->ghostk;
    if (diff == 0)
    {
      diff = ((face_sort_t *)lhs)->nf - ((face_sort_t *)rhs)->nf;
    }
  }
  return diff;
}

static int refine_level = 0;
static int refine_zlevel = 0;

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

/* This function fills COORD and INTMA arrays for NUMA */
static void fill_coordinates_p6est(int Nrp, int Nrpv, real *r, real *rz,
                                   p6est_t *p6est, p6est_lnodes_t *lnodes,
                                   p4esttonuma_t *p2n)
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
  size_t num_cols, zz, zy, first, last;
  p4est_topidx_t jt;
  p4est_locidx_t quad_count, Ntotal;
  sc_array_t *columns;
  p4est_tree_t *tree;
  p4est_quadrant_t *col;
  p2est_quadrant_t *layer;

  int column_size = 0;

  real w[P8EST_CHILDREN];

  /* Allocate coordinates in DG storage */
  Ntotal = (int)p6est->layers->elem_count;
  p2n->coord_dg = (real *)malloc(sizeof(real) * 3 * Ntotal * Nrp * Nrp * Nrpv);
  p2n->coord_cg = (real *)malloc(sizeof(real) * 3 * lnodes->num_local_nodes);
  p2n->intma = malloc(sizeof(p4est_locidx_t) * lnodes->vnodes * Ntotal);

  // loop over local trees
  for (jt = first_local_tree, quad_count = 0; jt <= last_local_tree; ++jt)
  {
    tree = p4est_tree_array_index(trees, jt);
    columns = &tree->quadrants;
    num_cols = columns->elem_count;
    p6est_tree_get_vertices(connectivity, jt, v);

    // loop over columns in each tree
    for (zz = 0; zz < num_cols; ++zz)
    {
      col = p4est_quadrant_array_index(columns, zz);
      P6EST_COLUMN_GET_RANGE(col, &first, &last);

      column_size = (int)(last - first);
      // printf("%s %d\n","-------> col_size",column_size);

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

              /* intma(1:nglx,1:ngly,1:nglz,1:nelem): gives point number
                 (1...npoin) for
                 each point in each element (ngl*=nop*+1, nop: polynomial
                 degree)*/
              p2n->intma[intma_index] = cg_index + 1;
            }
          }
        } // loop over element dof

      } // loop over elements (layers) in each columns

    } // loop over columns in each tree

  } // loop over local trees

  // FIXME: This assumes that all column_size are same size!
  // FIXME: Won't work for AMR
  P4EST_ASSERT(column_size > 0);
  p2n->nelz = column_size;
  p2n->nelxy = Ntotal / column_size;
}

/* Function added by MAK from p6est/test/test_all.c */

char test_data = 'x';
char *TEST_USER_POINTER = &test_data;

void init_fn_p6est(p6est_t *p6est, p4est_topidx_t which_tree,
                   p4est_quadrant_t *col, p2est_quadrant_t *layer)
{
  SC_CHECK_ABORT(p6est->user_pointer == TEST_USER_POINTER,
                 "user_pointer corruption\n");
}

/* To define a p6est_refine_column_t, all we have to do is take a p4est refine
 * function for uniform refinement... */
static int p4est_refine_fn(p4est_t *p4est, p4est_topidx_t which_tree,
                           p4est_quadrant_t *quadrant)
{

  /* printf("%s %d","---> refine_level",refine_level); */

  /* if (quadrant->level >= refine_level) { */
  /*   return 0; */
  /* } */

  return 1;
}

/* and wrap it.*/
static int refine_column_fn(p6est_t *p6est, p4est_topidx_t which_tree,
                            p4est_quadrant_t *column)
{
  return p4est_refine_fn(p6est->columns, which_tree, column);
}

static int refine_layer_fn(p6est_t *p6est, p4est_topidx_t which_tree,
                           p4est_quadrant_t *column, p2est_quadrant_t *layer)
{
  /* p4est_topidx_t      tohash[4]; */
  /* unsigned            hash; */

  /* tohash[0] = (p4est_topidx_t) column->x; */
  /* tohash[1] = (p4est_topidx_t) column->y; */
  /* tohash[2] = (p4est_topidx_t) layer->z; */
  /* tohash[3] = (((p4est_topidx_t) column->level) << 16) | */
  /*   ((p4est_topidx_t) layer->level); */

  /* hash = p4est_topidx_hash4 (tohash); */

  //  return (layer->level < refine_zlevel);// && !((int) hash % 3));
  return 1;
}

/* Keep state of p4est */
static p6est_t *stored_p6est = NULL;
extern mpi_context_t mpi_context;

static int FACE_LEN;
static int ORIENT;

/* End of p6est section */

#ifdef ibmcompiler
void p6esttonuma_init
#else
void p6esttonuma_init_
#endif
    (int *nref_levs, int *nref_zlevs, int *num_hroot, int *num_zroot, int *nnop,
     int *nnopz, real *xgl, real *xglz, int *is_cube, int *nnx, int *nny,
     int *nnz, int *nis_dg, int *xperiodic_flg, int *yperiodic_flg,
     int *mFACE_LEN, int *mORIENT, int* iboundary, p4esttonuma_t **p2n)
{

  /* Here goes the p6est stuff */

  p4est_connectivity_t *conn4, *conn_init;
  p6est_connectivity_t *conn;
  p6est_t *p6est;
  p6est_ghost_t *ghost_p6est;
  p6est_lnodes_t *lnodes_p6est;
  sc_array_t *columns;
  p4est_tree_t *tree;
  p4est_topidx_t jt;

  /* And here goes the p4est stuff */

  mpi_context_t *mpi = &mpi_context;
  p4est_t *p4est;
  p4est_ghost_t *ghost_p4est;
  p4est_mesh_t *mesh_p4est;

  int nref_levels, nref_zlevels;
  int nop, nopz, num_cols;
  int zz, i;
  int nx, ny, nz;
  int num_send_recv_total, num_send_recv_offset, num_nbh;
  p4est_locidx_t q, nq, sk, skb;
  int f, nf;
  int is_dg;
  p4est_locidx_t num_physical_boundary_faces;
  p4est_locidx_t num_processor_boundary_faces;

  nref_levels = *nref_levs;
  nref_zlevels = *nref_zlevs;
  int num_hroots = *num_hroot;
  int num_zroots = *num_zroot;
  nop = *nnop;
  nopz = *nnopz;
  nx = *nnx;
  ny = *nny;
  nz = *nnz;
  is_dg = *nis_dg;
  int xperiodic = *xperiodic_flg;
  int yperiodic = *yperiodic_flg;
  int is_geometry_cube = *is_cube;
  FACE_LEN = *mFACE_LEN;
  ORIENT = *mORIENT;

  *p2n = malloc(sizeof(p4esttonuma_t));

  /* Select p8est connectivity based on space_method*/
  p4est_connect_type_t connect_type;
  connect_type = P4EST_CONNECT_FULL;

  /*
   * Here goes the p6est stuff
   */
  {
    if (!stored_p6est)
    {
      if (!is_geometry_cube)
      {
        // create connectivity and forest structures
        conn_init = p4est_connectivity_new_cubed_sphere();

        // conn4 = p4est_connectivity_new_cubed();
        conn4 = p4est_connectivity_refine(conn_init, num_hroots);
        // printf("%s","conn4 created\n");

        // create top vertices by scaling original vertices
        double *topVertices = malloc(sizeof(double) * conn4->num_vertices * 3);
        for (i = 0; i < conn4->num_vertices * 3; ++i)
        {
          topVertices[i] = 2 * conn4->vertices[i];
        }
        // printf("%s","vertices created\n");

        // create p6est connectivity
        conn = p6est_connectivity_new(conn4, topVertices, NULL);
        // printf("%s","connectivity created\n");

        p4est_connectivity_destroy(conn4);
        p4est_connectivity_destroy(conn_init);
      }
      else
      {
        // here create flow-in-the-box mesh
        conn4 = p4est_connectivity_new_brick(nx, ny, xperiodic, yperiodic);
        double height[3] = {0., 0., 1.};
        conn = p6est_connectivity_new(conn4, NULL, height);
        p4est_connectivity_destroy(conn4);
        num_zroots = nz; // redundant information in the input!!
      }

      /* FXG modification suggested by Alex */
      int log_zroot;
      if (num_zroots > 1)
        log_zroot = (int)floor(log2(num_zroots - 1)) + 1;
      else
        log_zroot = 0;
      (*p2n)->log_zroots = log_zroot;

      p6est = p6est_new_ext(mpi->mpicomm, conn, 0, 0, log_zroot, num_zroots, 1,
                            3, init_fn_p6est, TEST_USER_POINTER);
      // printf("%s","p6est created\n");
    }
    else
    {
      p6est = stored_p6est;
      conn = p6est->connectivity;
    }

    // refinement levels
    refine_level = nref_levels;
    refine_zlevel = nref_zlevels;

    /*
     * Here goes the refinement after initial mesh creation
     */

    // printf("[p6estnuma] %s %d %d\n","------> refine horizontal, vertical",
    // 	nref_levels,nref_zlevels);

    // refine in horizontal
    for (i = 0; i < nref_levels; ++i)
      p6est_refine_columns(p6est, 0, refine_column_fn, init_fn_p6est);

    // refine in vertical
    for (i = 0; i < nref_zlevels; ++i)
      p6est_refine_layers(p6est, 0, refine_layer_fn, init_fn_p6est);

    // partition grid
    p6est_partition(p6est, NULL);

// Write vtk file
//#ifdef VTK_OUTPUT
    p6est_vtk_write_file(p6est, "p6est_test_initial"); /* Commented by FXG */
    //#endif

    /* Create the ghost layer to learn about parallel neighbors. */
    ghost_p6est = p6est_ghost_new(p6est, connect_type);

    /* Create a node numbering for continuous linear finite elements. */
    lnodes_p6est = p6est_lnodes_new(p6est, ghost_p6est, nop);

    p4est = p6est->columns;
    ghost_p4est = p4est_ghost_new(p4est, connect_type);

    mesh_p4est = p4est_mesh_new(p4est, ghost_p4est, connect_type);
    /* Get polynomial order */
    (*p2n)->nop = nop;
    (*p2n)->nopz = nopz;

    /* Get local (on rank) number of grid points */
    int npoin_p6est = lnodes_p6est->num_local_nodes;
    (*p2n)->npoin = npoin_p6est;

    /* Get nelem: number of elements on the processor */
    int nelem_p6est = (int)p6est->layers->elem_count;
    (*p2n)->nelem = nelem_p6est;
    (*p2n)->npoin_dg = (*p2n)->nelem * (nop + 1) * (nop + 1) * (nopz + 1);

    /* Get dg coordinates and intma */
    fill_coordinates_p6est(nop + 1, nopz + 1, xgl, xglz, p6est, lnodes_p6est,
                           *p2n);

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

    (*p2n)->node_column = malloc(sizeof(int) * (*p2n)->ncol * (*p2n)->nz);

    // loop over columns to fill node_column
    for (zz = 0; zz < (*p2n)->npoin; ++zz)
      (*p2n)->node_column[zz] =
          zz + 1; // simplification due to inherent p6est node numbering scheme

    /*
     * Count boundary faces
     */
    num_physical_boundary_faces = 0;
    num_processor_boundary_faces = 0;

    if (nopz > 0)
    {
      const p4est_topidx_t first_local_tree = p6est->columns->first_local_tree;
      const p4est_topidx_t last_local_tree = p6est->columns->last_local_tree;

      for (jt = first_local_tree; jt <= last_local_tree; ++jt)
      {
        tree = p4est_tree_array_index(p6est->columns->trees, jt);
        columns = &tree->quadrants;
        num_cols = (int)columns->elem_count;
        // each column has a boundary at the top and bottom only
        num_physical_boundary_faces += 2 * num_cols;
      }
    }
    /*
     * nbsido: number of (physical) domain boundary faces on the processor
     */

    // Count processor boundary faces
    for (q = 0; q < mesh_p4est->local_num_quadrants; ++q)
    {
      for (f = 0; f < P4EST_FACES; ++f)
      {
        nq = mesh_p4est->quad_to_quad[P4EST_FACES * q + f];
        nf = mesh_p4est->quad_to_face[P4EST_FACES * q + f];

        if (nq == q && nf == f && is_geometry_cube)
        {
          num_physical_boundary_faces += (*p2n)->nelz;
        }
        else if (nq >= mesh_p4est->local_num_quadrants)
        {
          num_processor_boundary_faces += (*p2n)->nelz;
        }
      }
    }

    (*p2n)->nbsido = num_physical_boundary_faces;

    /*
     * nface: number of unique faces on the processor (interior faces +
     * processor
     * boundary faces + domain boundary faces, see domain_decomp_metis.f90, line
     * 1438)
     */
    (*p2n)->nface = (P8EST_FACES * nelem_p6est - num_physical_boundary_faces -
                     num_processor_boundary_faces) /
                        2 +
                    num_physical_boundary_faces + num_processor_boundary_faces;

    /*
     * nboun: number of domain and processor boundary faces on the processor
     * (processor boundary faces + domain boundary faces, see
     * domain_decomp_metis.f90, line 1440)
     */
    (*p2n)->nboun =
        /*num_physical_boundary_faces +*/ num_processor_boundary_faces;

    /*/
      printf("%s\n","-----------------");
      printf("%s\n","-----------------");
      printf("%d %d\n",P8EST_FACES , nelem_p6est);
      printf("%s %d\n","proc_faces",num_processor_boundary_faces);
      printf("%s %d\n","phys_faces",num_physical_boundary_faces);
      printf("%s %d\n","nface",(*p2n)->nface);
      printf("%s %d\n","nboun",(*p2n)->nboun);
      printf("%s\n","-----------------");
      printf("%s\n","-----------------");
    */

    my_rank = mpi->mpirank;

    /**********************
     * create faces
     ***********************/
    p4est_locidx_t ncolz = (*p2n)->nelz;
    p4est_locidx_t nquads = mesh_p4est->local_num_quadrants;
    p4est_locidx_t nghosts = mesh_p4est->ghost_num_quadrants;
    p4est_locidx_t *intma_face =
        malloc(sizeof(p4est_locidx_t) * P8EST_FACES * ncolz * nquads);
    face_sort_t *mghost = malloc(sizeof(face_sort_t) * 8 * ncolz * nghosts);

    memset(intma_face, 0,
           sizeof(p4est_locidx_t) * P8EST_FACES * ncolz * nquads);

    /*face orientation transformation array p4est->numa*/
    static const int transform[6] = {4, 5, 2, 3, 0, 1};

    /*fill up intma_face*/
    int total_ghosts = 0, qq;
    for (q = 0, sk = 0; q < nquads; ++q)
    {
      for (f = 0; f < P4EST_FACES; ++f)
      {
        nq = mesh_p4est->quad_to_quad[P4EST_FACES * q + f];
        nf = mesh_p4est->quad_to_face[P4EST_FACES * q + f];

        if (nq == q && nf == f)
        {
          // should not be here if not a cube geometry
          for (zz = 0; zz < ncolz; zz++)
          {
            qq = q * ncolz + zz;
            intma_face[P8EST_FACES * qq + f] = sk;
            sk++;
          }
        }
        else if (nq >= nquads)
        {
          p4est_locidx_t ghostid = nq - nquads;
          p4est_locidx_t ghostp = mesh_p4est->ghost_to_proc[ghostid];
          p4est_quadrant_t *ghostquad =
              p4est_quadrant_array_index(&ghost_p4est->ghosts, (size_t)ghostid);
          p4est_locidx_t ghostk = ghostquad->p.piggy3.local_num;

          for (zz = 0; zz < ncolz; zz++)
          {
            qq = q * ncolz + zz;
            intma_face[P8EST_FACES * qq + f] = sk;

            mghost[total_ghosts].proc = ghostp;
            if (my_rank > ghostp)
            { // sort by neighbor's
              int nqq = ghostk * ncolz + zz;
              mghost[total_ghosts].ghostk = nqq;
              mghost[total_ghosts].nf = nf;
            }
            else
            { // sort by ours
              mghost[total_ghosts].ghostk = qq;
              mghost[total_ghosts].nf = f;
            }
            mghost[total_ghosts].face = sk;
            total_ghosts++;

            sk++;
          }
        }
        else if (q < nq)
        {
          for (zz = 0; zz < ncolz; zz++)
          {
            qq = q * ncolz + zz;
            intma_face[P8EST_FACES * qq + f] = sk;
            qq = nq * ncolz + zz;
            intma_face[P8EST_FACES * qq + nf] = sk;
            sk++;
          }
        }
      }

      // top-bottom of cubes
      if (nopz > 0)
      {
        for (zz = 0; zz < ncolz; zz++)
        {
          qq = q * ncolz + zz;
          if (zz == 0)
            intma_face[P8EST_FACES * qq + 4] = sk++;
          if (zz < ncolz - 1)
          {
            intma_face[P8EST_FACES * qq + 5] = sk;
            intma_face[P8EST_FACES * (qq + 1) + 4] = sk;
            sk++;
          }
          if (zz == ncolz - 1)
            intma_face[P8EST_FACES * qq + 5] = sk++;
        }
      }
    }

    /*
     * bsido(1:6,1:nbsido): domain boundary data:
     * bsido(1:4,i): point number (1...npoin) of the four corners of the face i
     *               (counter-clock wise when seen from outside of the domain)
     * bsido(5,i): element (1...nelem) to which face i belongs (only one because
     *             boundary face!)
     * bsido(6,i): boundary condition flag for face i (=4 for solid wall
     * boundary)
     */
    (*p2n)->bsido = malloc(sizeof(p4est_locidx_t) * 6 * (*p2n)->nbsido);
    memset((*p2n)->bsido, 0, sizeof(p4est_locidx_t) * 6 * (*p2n)->nbsido);
    /*
     * face(1:8, 1:nface): face data:
     * face(1:4,i): point number (1...npoin) of the four corners of the face i
     *               (counter-clock wise when seen from outside of the domain)
     * face(5,i): local face relative to left element
     * face(6,i): local face relative to right element
     * face(7,i): element to left
     * face(8,i): element to the right or boundary condition flag (-bsido(6,sk))
     */
    (*p2n)->face = malloc(sizeof(p4est_locidx_t) * FACE_LEN * (*p2n)->nface);
    (*p2n)->face_type = malloc(sizeof(p4est_locidx_t) * (*p2n)->nface);

    memset((*p2n)->face, 0, sizeof(p4est_locidx_t) * FACE_LEN * (*p2n)->nface);
    memset((*p2n)->face_type, 0, sizeof(p4est_locidx_t) * (*p2n)->nface);

    for (q = 0, skb = 0; q < nquads; ++q)
    {
      for (f = 0; f < P4EST_FACES; ++f)
      {
        nq = mesh_p4est->quad_to_quad[P4EST_FACES * q + f];
        nf = mesh_p4est->quad_to_face[P4EST_FACES * q + f];

        if (nq == q && nf == f)
        {
          // should not be here if not a cube geometry
          for (zz = 0; zz < ncolz; zz++)
          {
            qq = q * ncolz + zz;
            sk = intma_face[P8EST_FACES * qq + f];
            (*p2n)->face_type[sk] = 4;

            (*p2n)->face[FACE_LEN * sk + 4] = transform[f] + 1;
            (*p2n)->face[FACE_LEN * sk + 5] = 0;
            (*p2n)->face[FACE_LEN * sk + 6] = qq + 1;
            (*p2n)->face[FACE_LEN * sk + 7] = -iboundary[transform[f]];

            (*p2n)->bsido[6 * skb + 4] = qq + 1;
            (*p2n)->bsido[6 * skb + 5] = iboundary[transform[f]];
            skb++;
          }
        }
        else if (nq >= nquads)
        {
          for (zz = 0; zz < ncolz; zz++)
          {
            qq = q * ncolz + zz;
            sk = intma_face[P8EST_FACES * qq + f];
            (*p2n)->face_type[sk] = 2;

            (*p2n)->face[FACE_LEN * sk + 4] = transform[f] + 1;
            (*p2n)->face[FACE_LEN * sk + 5] = 0;
            (*p2n)->face[FACE_LEN * sk + 6] = qq + 1;
            (*p2n)->face[FACE_LEN * sk + 7] = 0;
          }
        }
        else if (q < nq || (q == nq && f < nf))
        {
          for (zz = 0; zz < ncolz; zz++)
          {
            qq = q * ncolz + zz;
            sk = intma_face[P8EST_FACES * qq + f];
            (*p2n)->face_type[sk] = 1;

            (*p2n)->face[FACE_LEN * sk + 4] = transform[f] + 1;
            (*p2n)->face[FACE_LEN * sk + 6] = qq + 1;
          }
        }
        else
        {
          for (zz = 0; zz < ncolz; zz++)
          {
            qq = q * ncolz + zz;
            sk = intma_face[P8EST_FACES * qq + f];
            (*p2n)->face_type[sk] = 1;

            (*p2n)->face[FACE_LEN * sk + 5] = transform[f] + 1;
            (*p2n)->face[FACE_LEN * sk + 7] = qq + 1;
          }
        }
      }

      // top-bottom of cubes
      if (nopz > 0)
      {
        for (zz = 0; zz < ncolz; zz++)
        {
          qq = q * ncolz + zz;
          if (zz == 0)
          {
            sk = intma_face[P8EST_FACES * qq + 4];
            (*p2n)->face_type[sk] = 4;

            (*p2n)->face[FACE_LEN * sk + 4] = transform[4] + 1;
            (*p2n)->face[FACE_LEN * sk + 6] = qq + 1;
            (*p2n)->face[FACE_LEN * sk + 5] = 0;
            (*p2n)->face[FACE_LEN * sk + 7] = -iboundary[transform[4]];

            (*p2n)->bsido[6 * skb + 4] = qq + 1;
            (*p2n)->bsido[6 * skb + 5] = iboundary[transform[4]];
            skb++;
          }
          if (zz < ncolz - 1)
          {
            sk = intma_face[P8EST_FACES * qq + 5];
            (*p2n)->face_type[sk] = 1;

            (*p2n)->face[FACE_LEN * sk + 4] = transform[5] + 1;
            (*p2n)->face[FACE_LEN * sk + 6] = qq + 1;
            (*p2n)->face[FACE_LEN * sk + 5] = transform[4] + 1;
            (*p2n)->face[FACE_LEN * sk + 7] = (qq + 1) + 1;
          }
          if (zz == ncolz - 1)
          {
            sk = intma_face[P8EST_FACES * qq + 5];
            (*p2n)->face_type[sk] = 4;

            (*p2n)->face[FACE_LEN * sk + 4] = transform[5] + 1;
            (*p2n)->face[FACE_LEN * sk + 6] = qq + 1;
            (*p2n)->face[FACE_LEN * sk + 5] = 0;
            (*p2n)->face[FACE_LEN * sk + 7] = -iboundary[transform[5]];

            (*p2n)->bsido[6 * skb + 4] = qq + 1;
            (*p2n)->bsido[6 * skb + 5] = iboundary[transform[5]];
            skb++;
          }
        }
      }
    }
    /*/
    //faces
    for(sk = 0; sk < (*p2n)->nface; ++sk)
    {
    printf("%4d. ", sk+1);
    for(f = 0;f < 8;f++)
    printf("%7d ",(*p2n)->face[sk*FACE_LEN+f]);
    printf("%7d ", (*p2n)->face_type[sk]);
    printf("\n");
    }
    //bsido
    for(sk = 0; sk < (*p2n)->nbsido; ++sk)
    {
    printf("%4d. ", sk+1);
    for(f = 0;f < 6;f++)
    printf("%7d ",(*p2n)->bsido[sk*6+f]);
    printf("\n");
    }
    */
    /**************************************
     * Parallel processing stuff goes here
     **************************************/
    num_nbh = (int)lnodes_p6est->sharers->elem_count;
    (*p2n)->num_send_recv = malloc(sizeof(int) * num_nbh);
    (*p2n)->nbh_proc = malloc(sizeof(int) * num_nbh);
    /*
     * Fill the sharing information for CG
     */
    if (!is_dg)
    {
      for (zz = 0, num_send_recv_total = 0, num_nbh = 0;
           zz < (int)lnodes_p6est->sharers->elem_count; ++zz)
      {
        p6est_lnodes_rank_t *lrank =
            p6est_lnodes_rank_array_index(lnodes_p6est->sharers, zz);

        if (lrank->rank != mpi->mpirank)
        {
          (*p2n)->nbh_proc[num_nbh] = lrank->rank + 1;
          (*p2n)->num_send_recv[num_nbh] = (int)lrank->shared_nodes.elem_count;

          num_send_recv_total += (int)lrank->shared_nodes.elem_count;
          num_nbh += 1;
        }
      }

      /* num_send_recv_total: total number of points to be communicated with all
       * neighboring processors
       */
      (*p2n)->num_send_recv_total = num_send_recv_total;

      /* num_nbh: number of neighboring processors */
      (*p2n)->num_nbh = num_nbh;

      /* nbh_send_recv(1:num_send_recv_total): lcoal number of the point that
       * must
       * be sent/received sorted in a global ordering (e.g.:
       * nbh_send_recv(1:num_send_recv(1)) are the points to be communicated
       * with
       * processor nbh_proc(1).
       * nbh_send_recv(num_send_recv(1)+1:num_send_recv(1)+num_send_recv(2)) are
       * the points to be communicated with processor nbh_proc(2)...)
       */

      (*p2n)->nbh_send_recv =
          malloc(sizeof(p4est_locidx_t) * num_send_recv_total);
      for (zz = 0, num_send_recv_offset = 0;
           zz < (int)lnodes_p6est->sharers->elem_count; ++zz)
      {
        p6est_lnodes_rank_t *lrank =
            p6est_lnodes_rank_array_index(lnodes_p6est->sharers, zz);
        const int nshared = (int)lrank->shared_nodes.elem_count;
        int jn;

        if (lrank->rank != mpi->mpirank)
        {
          for (jn = 0; jn < nshared; ++jn)
          {
            int gl = (int)*(p4est_locidx_t *)sc_array_index_int(
                &lrank->shared_nodes, jn);
            (*p2n)->nbh_send_recv[num_send_recv_offset + jn] = gl + 1;
          }
          num_send_recv_offset += nshared;
        }
      }
    }
    else
    {
      /*
       * Fill the sharing information for DG
       */
      (*p2n)->nbh_send_recv = malloc(sizeof(p4est_locidx_t) * total_ghosts);

      qsort(mghost, total_ghosts, sizeof(face_sort_t), cmpfunc);

      int p = -1;
      for (zz = 0, num_nbh = 0; zz < total_ghosts; zz++)
      {
        int proc = mghost[zz].proc;
        int face = mghost[zz].face;
        if (proc != p)
        {
          num_nbh++;
          (*p2n)->nbh_proc[num_nbh - 1] = proc + 1;
          (*p2n)->num_send_recv[num_nbh - 1] = 1;
          p = proc;
        }
        else
        {
          (*p2n)->num_send_recv[num_nbh - 1]++;
        }
        (*p2n)->nbh_send_recv[zz] = face + 1;
      }
      (*p2n)->num_send_recv_total = total_ghosts;
      (*p2n)->num_nbh = num_nbh;
    }
    /*
       printf("NumSend %d\n",(*p2n)->num_send_recv_total);
       for(i = 0; i < (*p2n)->num_send_recv_total; ++i)
       printf("%d ",(*p2n)->nbh_send_recv[i]);
       printf("\n");
     *
     printf("===========================\n");
     printf(" Processor %d\n",my_rank);
     printf("===========================\n");
     for(sk = 0; sk < total_ghosts; ++sk)
     {
     printf("%4d. ", sk+1);
     printf("%7d %7d %7d %7d\n",mghost[sk].proc,mghost[sk].ghostk,
     mghost[sk].nf,mghost[sk].face);
     }
    */

    /*free local vars*/
    free(intma_face);
    free(mghost);

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
    p6est_lnodes_destroy(lnodes_p6est);
    p6est_ghost_destroy(ghost_p6est);

    stored_p6est = p6est;
  }
  /* End of p6est section */
}

#ifdef ibmcompiler
void p6esttonuma_get_mesh_scalars
#else
void p6esttonuma_get_mesh_scalars_
#endif
    (p4esttonuma_t **p2n, int *npoin, int *nelem, int *log_zroots, int *num_nbh,
     int *num_send_recv_total, int *nbsido, int *nface, int *nboun, int *ncol,
     int *nz, int *nelz)
{
  *npoin = (*p2n)->npoin;
  *nelem = (*p2n)->nelem;
  *num_nbh = (*p2n)->num_nbh;
  *num_send_recv_total = (*p2n)->num_send_recv_total;
  *nbsido = (*p2n)->nbsido;
  *nface = (*p2n)->nface;
  *nboun = (*p2n)->nboun;
  *ncol = (*p2n)->ncol;
  *nz = (*p2n)->nz;
  *nelz = (*p2n)->nelz;
  *log_zroots = (*p2n)->log_zroots;
}

#ifdef ibmcompiler
void p6esttonuma_get_mesh_arrays
#else
void p6esttonuma_get_mesh_arrays_
#endif
    (p4esttonuma_t **p2n, real *coord, int *intma, int *face, int *face_type,
     int *num_send_recv, int *nbh_proc, int *nbh_send_recv, int *bsido,
     int *node_column, int *nis_cgc)
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
    intma[i] = (*p2n)->intma[i];

  for (i = 0; i < (*p2n)->num_nbh; ++i)
    num_send_recv[i] = (*p2n)->num_send_recv[i];

  for (i = 0; i < (*p2n)->num_nbh; ++i)
    nbh_proc[i] = (*p2n)->nbh_proc[i];

  for (i = 0; i < (*p2n)->num_send_recv_total; ++i)
    nbh_send_recv[i] = (*p2n)->nbh_send_recv[i];

  for (i = 0; i < 6 * (*p2n)->nbsido; ++i)
    bsido[i] = (*p2n)->bsido[i];

  for (i = 0; i < FACE_LEN * (*p2n)->nface; ++i)
    face[i] = (*p2n)->face[i];

  for (i = 0; i < (*p2n)->nface; ++i)
    face_type[i] = (*p2n)->face_type[i];

  for (i = 0; i < ((*p2n)->nz) * ((*p2n)->ncol); ++i)
    node_column[i] = (*p2n)->node_column[i];
}

#ifdef ibmcompiler
void p6esttonuma_free
#else
void p6esttonuma_free_
#endif
    (p4esttonuma_t **p2n)
{
  free((*p2n)->coord_dg);
  free((*p2n)->coord_cg);
  free((*p2n)->intma);
  free((*p2n)->num_send_recv);
  free((*p2n)->nbh_proc);
  free((*p2n)->nbh_send_recv);
  free((*p2n)->bsido);
  free((*p2n)->face);
  free((*p2n)->face_type);
  free((*p2n)->node_column);
  free(*p2n);
}

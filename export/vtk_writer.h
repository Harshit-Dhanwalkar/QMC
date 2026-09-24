#ifndef QMC_VTK_WRITER_H
#define QMC_VTK_WRITER_H

#include <stddef.h>

// NOTE: Write a scalar field on a regular structured grid as a legacy
// VTK ASCII file (STRUCTURED_POINTS dataset) - openable directly in
// ParaView or VisIt (no conversion step needed)
//
// `data` is C-order (row-major) with x fastest-varying index:
//   data[k * ny * nx + j * nx + i] is value at grid point (i, j, k)
// For a 2D field, pass nz = 1
//
// spacing_{x,y,z} is grid spacing (dx, dy, dz); origin_{x,y,z} is
// physical-space coordinate of grid point (0, 0, 0)
//
// Returns 0 on success, -1 if the output file could not be opened
int vtk_write_structured_points(const char *filename, const char *varname,
                                const double *data, size_t nx, size_t ny,
                                size_t nz, double spacing_x, double spacing_y,
                                double spacing_z, double origin_x,
                                double origin_y, double origin_z);

#endif // QMC_VTK_WRITER_H

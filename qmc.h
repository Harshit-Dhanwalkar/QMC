#ifndef QMC_H
#define QMC_H

/*
 * QMC umbrella header: one include for downstream projects
 *
 *   #include <qmc/qmc.h>
 *
 * Pulls in version macros, every physics module (via physics/physics.h, which
 * also brings in core/), always-available exporters (CSV, JSON, VTK) and LaTeX
 * generator
 *
 * Deliberately NOT included, because they need optional third-party headers:
 *   <qmc/export/hdf5_writer.h>     (build with USE_HDF5=1)
 *   <qmc/export/netcdf_writer.h>   (build with USE_NETCDF=1)
 *   <qmc/export/plot.h>            (plotting backend)
 * Include those directly where needed
 */

// IWYU pragma: begin_exports
#include "core/qmc_version.h"
#include "export/csv_writer.h"
#include "export/export_schema.h"
#include "export/json_writer.h"
#include "export/vtk_writer.h"
#include "latex/latex_gen.h"
#include "physics/physics.h"
// IWYU pragma: end_exports

#endif /* QMC_H */

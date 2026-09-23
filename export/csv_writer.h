#ifndef QMC_CSV_WRITER_H
#define QMC_CSV_WRITER_H

// Write 1D data (x, y) to CSV with a header row 
// Returns 0 on success, -1 if output file could not be opened
// Path is QMC_OUTPUT_DIR/filename
int csv_write_1d(const char *filename, const double *x, const double *y, int n,
                 const char *xlabel, const char *ylabel);

// Write a row-major 2D matrix to CSV. col_headers may be NULL, in which case
// columns are labeled col0, col1, ...
// Returns 0 on success, -1 on failure
int csv_write_matrix(const char *filename, const double *data, int rows,
                     int cols, const char **col_headers);

#endif // QMC_CSV_WRITER_H

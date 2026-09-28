/* Bounded nearest-PD and simplex QP for aggregated ancestry projections. */
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#define ANCESTRY_GROUPS 30
#define ANCESTRY_PCS 64
#define ANCESTRY_DIM (ANCESTRY_GROUPS + 1)

static double matrix_norm(const double *a, int n) {
    double maximum = 0;
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < n; j++) sum += fabs(a[i * n + j]);
        if (sum > maximum) maximum = sum;
    }
    return maximum;
}

/* Orthogonal Jacobi rotations for the small symmetric matrices used here. */
static bool eigen_symmetric(const double *a, int n, double *values, double *vectors) {
    double work[ANCESTRY_GROUPS * ANCESTRY_GROUPS];
    memcpy(work, a, (size_t)n * n * sizeof(double));
    memset(vectors, 0, (size_t)n * n * sizeof(double));
    for (int i = 0; i < n; i++) vectors[i * n + i] = 1;
    for (int sweep = 0; sweep < 100 * n * n; sweep++) {
        int p = 0, q = 1;
        double largest = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
            double v = fabs(work[i * n + j]);
            if (v > largest) { largest = v; p = i; q = j; }
        }
        if (largest <= 1e-14 * matrix_norm(work, n)) {
            for (int i = 0; i < n; i++) values[i] = work[i * n + i];
            return true;
        }
        double phi = 0.5 * atan2(2 * work[p * n + q], work[q * n + q] - work[p * n + p]);
        double c = cos(phi), s = sin(phi);
        for (int i = 0; i < n; i++) if (i != p && i != q) {
            double x = work[i * n + p], y = work[i * n + q];
            work[i * n + p] = work[p * n + i] = c * x - s * y;
            work[i * n + q] = work[q * n + i] = s * x + c * y;
        }
        double pp = work[p * n + p], qq = work[q * n + q], pq = work[p * n + q];
        work[p * n + p] = c * c * pp - 2 * c * s * pq + s * s * qq;
        work[q * n + q] = s * s * pp + 2 * c * s * pq + c * c * qq;
        work[p * n + q] = work[q * n + p] = 0;
        for (int i = 0; i < n; i++) {
            double x = vectors[i * n + p], y = vectors[i * n + q];
            vectors[i * n + p] = c * x - s * y;
            vectors[i * n + q] = s * x + c * y;
        }
    }
    return false;
}

static void reconstruct(double *out, const double *values, const double *vectors, int n) {
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
        double sum = 0;
        for (int k = 0; k < n; k++) sum += vectors[i * n + k] * values[k] * vectors[j * n + k];
        out[i * n + j] = sum;
    }
}

/* Matrix::nearPD defaults: Dykstra, eig.tol=1e-6, conv.tol=1e-7,
 * maxit=100; final posd.tol=1e-8 with diagonal restoration. */
static bool nearest_pd(double *a, int n) {
    double previous[ANCESTRY_GROUPS * ANCESTRY_GROUPS];
    double residual[ANCESTRY_GROUPS * ANCESTRY_GROUPS] = {0};
    double vectors[ANCESTRY_GROUPS * ANCESTRY_GROUPS], values[ANCESTRY_GROUPS];
    bool converged = false;
    for (int iteration = 0; iteration < 100; iteration++) {
        memcpy(previous, a, (size_t)n * n * sizeof(double));
        for (int i = 0; i < n * n; i++) a[i] -= residual[i];
        if (!eigen_symmetric(a, n, values, vectors)) return false;
        double largest = values[0];
        for (int i = 1; i < n; i++) if (values[i] > largest) largest = values[i];
        if (largest <= 0) return false;
        for (int i = 0; i < n; i++) if (values[i] <= 1e-6 * largest) values[i] = 0;
        double projected[ANCESTRY_GROUPS * ANCESTRY_GROUPS];
        reconstruct(projected, values, vectors, n);
        for (int i = 0; i < n * n; i++) {
            residual[i] = projected[i] - a[i];
            a[i] = projected[i];
        }
        double delta[ANCESTRY_GROUPS * ANCESTRY_GROUPS];
        for (int i = 0; i < n * n; i++) delta[i] = a[i] - previous[i];
        if (matrix_norm(delta, n) <= 1e-7 * matrix_norm(previous, n)) {
            converged = true;
            break;
        }
    }
    if (!converged || !eigen_symmetric(a, n, values, vectors)) return false;
    double largest = values[0];
    for (int i = 1; i < n; i++) if (values[i] > largest) largest = values[i];
    if (largest <= 0) return false;
    double floor_value = 1e-8 * largest;
    double original_diag[ANCESTRY_GROUPS];
    for (int i = 0; i < n; i++) {
        original_diag[i] = a[i * n + i];
        if (values[i] < floor_value) values[i] = floor_value;
    }
    reconstruct(a, values, vectors, n);
    double scale[ANCESTRY_GROUPS];
    for (int i = 0; i < n; i++) scale[i] = sqrt(fmax(floor_value, original_diag[i]) / a[i * n + i]);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) a[i * n + j] *= scale[i] * scale[j];
    return true;
}

/* Gaussian elimination with partial pivoting on an active face's KKT system. */
static bool linear_solve(double a[ANCESTRY_DIM + 1][ANCESTRY_DIM + 2], int n, double *out) {
    for (int col = 0; col < n; col++) {
        int pivot = col;
        for (int row = col + 1; row < n; row++)
            if (fabs(a[row][col]) > fabs(a[pivot][col])) pivot = row;
        if (fabs(a[pivot][col]) < 1e-14) return false;
        if (pivot != col) for (int j = col; j <= n; j++) {
            double tmp = a[col][j]; a[col][j] = a[pivot][j]; a[pivot][j] = tmp;
        }
        for (int row = col + 1; row < n; row++) {
            double ratio = a[row][col] / a[col][col];
            for (int j = col; j <= n; j++) a[row][j] -= ratio * a[col][j];
        }
    }
    for (int i = n - 1; i >= 0; i--) {
        double value = a[i][n];
        for (int j = i + 1; j < n; j++) value -= a[i][j] * out[j];
        out[i] = value / a[i][i];
    }
    return true;
}

/* Slack is a zero-cost extra coordinate for sum(q)<=1. */
static bool simplex_qp(const double *h, const double *d, int k, bool equality, double *q) {
    bool active[ANCESTRY_DIM] = {0};
    int n = k + !equality;
    int initial = equality ? 0 : k;
    if (equality) {
        for (int i = 1; i < k; i++) if (d[i] - h[i * k + i] / 2 >
            d[initial] - h[initial * k + initial] / 2) initial = i;
    }
    q[initial] = 1;
    active[initial] = true;
    for (int iteration = 0; iteration < 10000; iteration++) {
        int indices[ANCESTRY_DIM], m = 0;
        for (int i = 0; i < n; i++) if (active[i]) indices[m++] = i;
        double system[ANCESTRY_DIM + 1][ANCESTRY_DIM + 2] = {{0}};
        double target[ANCESTRY_DIM + 1] = {0};
        for (int i = 0; i < m; i++) {
            int row = indices[i];
            for (int j = 0; j < m; j++) {
                int col = indices[j];
                system[i][j] = row == k || col == k ? 0 : h[row * k + col];
            }
            system[i][m] = 1;
            system[i][m + 1] = row == k ? 0 : d[row];
        }
        for (int j = 0; j < m; j++) system[m][j] = 1;
        system[m][m + 1] = 1;
        if (!linear_solve(system, m + 1, target)) return false;
        double step = 1;
        for (int i = 0; i < m; i++) {
            int idx = indices[i];
            if (target[i] < 0 && q[idx] > 0) {
                double limit = q[idx] / (q[idx] - target[i]);
                if (limit < step) step = limit;
            }
        }
        for (int i = 0; i < m; i++) {
            int idx = indices[i];
            q[idx] += step * (target[i] - q[idx]);
            if (fabs(q[idx]) < 1e-12) q[idx] = 0;
        }
        if (step < 1 - 1e-12) {
            for (int i = 0; i < m; i++) if (q[indices[i]] == 0) active[indices[i]] = false;
            continue;
        }
        double multiplier = target[m], worst = -1e-9;
        int add = -1;
        for (int i = 0; i < n; i++) if (!active[i]) {
            double gradient = i == k ? 0 : -d[i];
            if (i != k) for (int j = 0; j < k; j++) gradient += h[i * k + j] * q[j];
            double reduced = gradient + multiplier;
            if (reduced < worst) { worst = reduced; add = i; }
        }
        if (add < 0) return true;
        active[add] = true;
    }
    return false;
}

static bool valid_at(duckdb_vector vector, idx_t row) {
    uint64_t *mask = duckdb_vector_get_validity(vector);
    return !mask || duckdb_validity_row_is_valid(mask, row);
}

static void ancestry_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector xvec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector yvec = duckdb_data_chunk_get_vector(input, 1);
    duckdb_vector kvec = duckdb_data_chunk_get_vector(input, 2);
    bool has_equality = duckdb_data_chunk_get_column_count(input) == 4;
    duckdb_vector eqvec = has_equality ? duckdb_data_chunk_get_vector(input, 3) : NULL;
    duckdb_vector xc = duckdb_list_vector_get_child(xvec);
    duckdb_vector yc = duckdb_list_vector_get_child(yvec);
    duckdb_list_entry *xe = duckdb_vector_get_data(xvec), *ye = duckdb_vector_get_data(yvec);
    int32_t *ks = duckdb_vector_get_data(kvec);
    bool *equalities = has_equality ? duckdb_vector_get_data(eqvec) : NULL;
    double *xs = duckdb_vector_get_data(xc), *ys = duckdb_vector_get_data(yc);
    duckdb_list_vector_set_size(output, 0);
    for (idx_t row = 0; row < duckdb_data_chunk_get_size(input); row++) {
        if (!valid_at(xvec, row) || !valid_at(yvec, row) || !valid_at(kvec, row) ||
            (has_equality && !valid_at(eqvec, row))) {
            duckdb_vector_ensure_validity_writable(output);
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(output), row);
            continue;
        }
        int k = ks[row];
        idx_t pcs = ye[row].length;
        if (k < 1 || k > ANCESTRY_GROUPS || pcs < 1 || pcs > ANCESTRY_PCS ||
            xe[row].length != pcs * (idx_t)k) {
            duckdb_scalar_function_set_error(info, "duckhts_ancestry_proportions: expected 1..30 groups and 1..64 PCs in PC-major X");
            return;
        }
        double h[ANCESTRY_GROUPS * ANCESTRY_GROUPS] = {0};
        double d[ANCESTRY_GROUPS] = {0}, q[ANCESTRY_DIM] = {0};
        double scale = 0;
        for (idx_t pc = 0; pc < pcs; pc++) {
            if (!valid_at(yc, ye[row].offset + pc) || !isfinite(ys[ye[row].offset + pc])) goto invalid;
            for (int i = 0; i < k; i++) {
                idx_t at = xe[row].offset + pc * (idx_t)k + (idx_t)i;
                if (!valid_at(xc, at) || !isfinite(xs[at])) goto invalid;
                scale = fmax(scale, fabs(xs[at]));
            }
        }
        if (scale == 0) {
            duckdb_scalar_function_set_error(info, "duckhts_ancestry_proportions: singular X");
            return;
        }
        for (idx_t pc = 0; pc < pcs; pc++) {
            double yi = ys[ye[row].offset + pc] / scale;
            if (!isfinite(yi)) goto invalid_scale;
            for (int i = 0; i < k; i++) {
                double xi = xs[xe[row].offset + pc * (idx_t)k + (idx_t)i] / scale;
                d[i] += xi * yi;
                for (int j = 0; j < k; j++) {
                    double xj = xs[xe[row].offset + pc * (idx_t)k + (idx_t)j] / scale;
                    h[i * k + j] += xi * xj;
                }
            }
        }
        for (int i = 0; i < k; i++) if (!isfinite(d[i])) goto invalid_scale;
        if (!nearest_pd(h, k) || !simplex_qp(h, d, k,
                                            has_equality ? equalities[row] : true, q)) {
            duckdb_scalar_function_set_error(info, "duckhts_ancestry_proportions: nearest-PD or QP did not converge");
            return;
        }
        {
            idx_t offset = duckdb_list_vector_get_size(output);
            if (duckdb_list_vector_reserve(output, offset + (idx_t)k) != DuckDBSuccess ||
                duckdb_list_vector_set_size(output, offset + (idx_t)k) != DuckDBSuccess) {
                duckdb_scalar_function_set_error(info, "duckhts_ancestry_proportions: output allocation failed");
                return;
            }
            duckdb_vector child = duckdb_list_vector_get_child(output);
            double *values = duckdb_vector_get_data(child);
            for (int i = 0; i < k; i++) values[offset + (idx_t)i] = q[i];
            ((duckdb_list_entry *)duckdb_vector_get_data(output))[row] = (duckdb_list_entry){offset, (idx_t)k};
        }
        continue;
invalid:
        duckdb_scalar_function_set_error(info, "duckhts_ancestry_proportions: X and y must be finite without NULLs");
        return;
invalid_scale:
        duckdb_scalar_function_set_error(info, "duckhts_ancestry_proportions: scaled X and y exceed finite range");
        return;
    }
}

void register_duckhts_ancestry_functions(duckdb_connection connection) {
    duckdb_scalar_function function = duckdb_create_scalar_function();
    duckdb_logical_type real = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_logical_type values = duckdb_create_list_type(real);
    duckdb_logical_type integer = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_logical_type boolean = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_scalar_function_set_name(function, "duckhts_ancestry_proportions");
    duckdb_scalar_function_add_parameter(function, values);
    duckdb_scalar_function_add_parameter(function, values);
    duckdb_scalar_function_add_parameter(function, integer);
    duckdb_scalar_function_add_parameter(function, boolean);
    duckdb_scalar_function_set_return_type(function, values);
    duckdb_scalar_function_set_special_handling(function);
    duckdb_scalar_function_set_function(function, ancestry_scalar);
    duckdb_register_scalar_function(connection, function);
    duckdb_destroy_scalar_function(&function);
    function = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(function, "duckhts_ancestry_proportions");
    duckdb_scalar_function_add_parameter(function, values);
    duckdb_scalar_function_add_parameter(function, values);
    duckdb_scalar_function_add_parameter(function, integer);
    duckdb_scalar_function_set_return_type(function, values);
    duckdb_scalar_function_set_special_handling(function);
    duckdb_scalar_function_set_function(function, ancestry_scalar);
    duckdb_register_scalar_function(connection, function);
    duckdb_destroy_scalar_function(&function);
    duckdb_destroy_logical_type(&real);
    duckdb_destroy_logical_type(&values);
    duckdb_destroy_logical_type(&integer);
    duckdb_destroy_logical_type(&boolean);
}

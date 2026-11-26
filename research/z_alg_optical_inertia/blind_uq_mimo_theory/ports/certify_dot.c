/* Explicit binary64 matrix product for a backward-error enclosure.
   Clang must compile with -fno-fast-math -ffp-contract=off. */
#include <pthread.h>
#include <stddef.h>

struct job {
    size_t row0, row1, inner, cols;
    const double *ar, *ai, *br, *bi;
    double *cr, *ci;
};

static void *multiply_rows(void *input) {
    const struct job *j = input;
    for (size_t k0 = 0; k0 < j->inner; k0 += 64) {
        const size_t k1 = k0 + 64 < j->inner ? k0 + 64 : j->inner;
        for (size_t c0 = 0; c0 < j->cols; c0 += 128) {
            const size_t c1 = c0 + 128 < j->cols ? c0 + 128 : j->cols;
            for (size_t i = j->row0; i < j->row1; i++) {
                double *out_re = j->cr + i * j->cols;
                double *out_im = j->ci + i * j->cols;
                for (size_t k = k0; k < k1; k++) {
                    const double a_r = j->ar[i*j->inner+k];
                    const double a_i = j->ai[i*j->inner+k];
                    const double *row_re = j->br + k * j->cols;
                    const double *row_im = j->bi + k * j->cols;
                    for (size_t c = c0; c < c1; c++) {
                        out_re[c] += a_r*row_re[c] - a_i*row_im[c];
                        out_im[c] += a_r*row_im[c] + a_i*row_re[c];
                    }
                }
            }
        }
    }
    return NULL;
}

int complex_product(size_t rows, size_t inner, size_t cols,
                    const double *ar, const double *ai, const double *br, const double *bi,
                    double *cr, double *ci, size_t threads) {
    if (!threads || threads > 16) {
        return 1;
    }
    pthread_t pool[16];
    struct job args[16];
    size_t launched = 0;
    for (size_t t=0; t<threads; t++) {
        args[t] = (struct job){rows*t/threads, rows*(t+1)/threads, inner, cols, ar, ai, br, bi, cr, ci};
        if (pthread_create(&pool[t], NULL, multiply_rows, &args[t])) {
            break;
        }
        launched++;
    }
    for (size_t t=0; t<launched; t++) {
        pthread_join(pool[t], NULL);
    }
    return launched != threads;
}

double binary64_sum_squares(size_t n, const double *real, const double *imag) {
    double value = 0.0;
    for (size_t i=0; i<n; i++) {
        value += real[i]*real[i] + imag[i]*imag[i];
    }
    return value;
}

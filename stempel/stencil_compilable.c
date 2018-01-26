#include <stdlib.h>
#include <math.h>

#include "timing.h"
#include "kerncraft.h"
#ifdef LIKWID_PERFMON
#include <likwid.h>
#endif

#define N 1767L
#define P 500L
#define M 1767L

void* aligned_malloc(size_t, size_t);
extern int var_false;
void dummy(double *);
void kernel_loop(double *a, double *b, double *W)
{
  #pragma omp parallel for schedule(runtime)
  for (int k = 2; k < (M - 2); k++)
  {
    for (int j = 2; j < (N - 2); j++)
    {
      for (int i = 2; i < (P - 2); i++)
      {
        b[(i + (j * P)) + (k * (P * N))] = ((W[((i + (j * P)) + (0 * (P * N))) + (k * ((P * N) * M))] * a[(i + (j * P)) + (k * (P * N))]) + (W[((i + (j * P)) + (1 * (P * N))) + (k * ((P * N) * M))] * (((a[((i - 1) + (j * P)) + (k * (P * N))] + a[((i + 1) + (j * P)) + (k * (P * N))]) + (a[(i + (j * P)) + ((k - 1) * (P * N))] + a[(i + (j * P)) + ((k + 1) * (P * N))])) + (a[(i + ((j - 1) * P)) + (k * (P * N))] + a[(i + ((j + 1) * P)) + (k * (P * N))])))) + (W[((i + (j * P)) + (2 * (P * N))) + (k * ((P * N) * M))] * (((a[((i - 2) + (j * P)) + (k * (P * N))] + a[((i + 2) + (j * P)) + (k * (P * N))]) + (a[(i + (j * P)) + ((k - 2) * (P * N))] + a[(i + (j * P)) + ((k + 2) * (P * N))])) + (a[(i + ((j - 2) * P)) + (k * (P * N))] + a[(i + ((j + 2) * P)) + (k * (P * N))])));
      }

    }

  }

  if (var_false)
  {
    dummy(a);
  }

  if (var_false)
  {
    dummy(b);
  }

  if (var_false)
  {
    dummy(W);
  }

}

int main(int argc, char **argv)
{
  
  #ifdef LIKWID_PERFMON
  LIKWID_MARKER_INIT;
  #endif

  
  #ifdef LIKWID_PERFMON
  #pragma omp parallel
  {
    LIKWID_MARKER_THREADINIT;
  }
  #endif

  double *a = aligned_malloc((sizeof(double)) * ((M * N) * P), 32);
  for (int i = 0; i < ((M * N) * P); ++i)
    a[i] = rand() / ((double ) RAND_MAX);

  if (var_false)
  {
    dummy(a);
  }

  double *b = aligned_malloc((sizeof(double)) * ((M * N) * P), 32);
  for (int i = 0; i < ((M * N) * P); ++i)
    b[i] = rand() / ((double ) RAND_MAX);

  if (var_false)
  {
    dummy(b);
  }

  double *W = aligned_malloc((sizeof(double)) * (((3 * M) * N) * P), 32);
  for (int i = 0; i < (((3 * M) * N) * P); ++i)
    W[i] = (rand() / ((double ) RAND_MAX)) / 5.0;

  if (var_false)
  {
    dummy(W);
  }

  int repeat = 1;
  double runtime = 0.0;
  double wct_start;
  double wct_end;
  double cput_start;
  double cput_end;
  double *tmp;
  
  #ifdef LIKWID_PERFMON
  #pragma omp parallel
  {
    LIKWID_MARKER_START("Sweep");
  }
  #endif

  double total = 0.0;
  while (runtime < 0.5)
  {
    timing(&wct_start, &cput_start);
    for (int n = 0; n < repeat; ++n)
    {
      kernel_loop(a, b, W);
      tmp = a;
      a = b;
      b = tmp;
    }

    timing(&wct_end, &cput_end);
    runtime = wct_end - wct_start;
    repeat *= 2;
  }

  
  #ifdef LIKWID_PERFMON
  #pragma omp parallel
  {
    LIKWID_MARKER_STOP("Sweep");
  }
  #endif

  repeat /= 2;
  printf("b'Performance in mlup/s: %lf \\n'", (((double ) repeat) * ((((double ) (M - 4)) * ((double ) (N - 4))) * ((double ) (P - 4)))) / (runtime * 1000000.));
  printf("b'size: %d    time: %lf    iter: %d    mlup/s: %lf\\n'", (M * P * N), runtime, repeat, (((double ) repeat) * ((((double ) (M - 4)) * ((double ) (N - 4))) * ((double ) (P - 4)))) / (runtime * 1000000.));
  for (int k = 2; k < (M - 2); k++)
  {
    for (int j = 2; j < (N - 2); j++)
    {
      for (int i = 2; i < (P - 2); i++)
      {
        total = total + (a[(i + (j * P)) + (k * (P * N))] * a[(i + (j * P)) + (k * (P * N))]);
      }

    }

  }

  printf("b'norm(a): %lf\\n'", sqrt(total));
  
  #ifdef LIKWID_PERFMON
  LIKWID_MARKER_CLOSE;
  #endif

}
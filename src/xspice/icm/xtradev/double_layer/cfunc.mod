#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define INPUT_SAMPLES 64

// Code model memory.
#define INP 0
#define IMPULSE 1

typedef struct {
  double x;
  double y;
} pair_t;

typedef struct {
  char name[256];
  pair_t *p_values;
  int n_values;
  int max_values;
  double x_min, x_max;
  double y_min, y_max;
} plot_t;

#define MAX(a, b) (a > b ? a : b)
#define MIN(a, b) (a < b ? a : b)

void add_to_plot(plot_t *p, double x, double y) {

  if (p->n_values == p->max_values) {
    p->max_values *= 2;
    pair_t *p_values = malloc(sizeof(pair_t) * p->max_values);
    memcpy(p_values, p->p_values, sizeof(pair_t) * p->n_values);
    free(p->p_values);
    p->p_values = p_values;
  }

  p->p_values[p->n_values].x = x;
  p->p_values[p->n_values].y = y;
  p->n_values++;
  if (p->n_values == 1) {
    p->x_min = p->x_max = x;
    p->y_min = p->y_max = y;
  } else {
    p->x_min = MIN(p->x_min, x);
    p->x_max = MAX(p->x_max, x);
    p->y_min = MIN(p->y_min, y);
    p->y_max = MAX(p->y_max, y);
  }
}

plot_t *new_plot(char *p_name) {
  plot_t *p = malloc(sizeof(plot_t));
  strncpy(p->name, p_name, sizeof(p->name));
  p->p_values = malloc(sizeof(pair_t) * 16);
  p->n_values = 0;
  p->max_values = 16;
  return p;
}

void end_plot(plot_t *p, int do_draw) {

  if (do_draw) {
    printf("\"%s\", %i points, x_{min,max} = (%.3e, %.3e), y_{min,max} = (%.3e, %.3e)\n", p->name, p->n_values, p->x_min, p->x_max, p->y_min, p->y_max);
    const int rows = 16, cols = 64;
    int r, c;
    for (r = rows; r > 0; r--) {

      double row_min = p->y_min + ((p->y_max - p->y_min) * ((r - 1) / (double)rows));
      double row_max = p->y_min + ((p->y_max - p->y_min) * (r / (double)rows));
      printf("%s%.3e, %s%.3e | ",
        row_min < 0 ? "" : " ", row_min,
        row_max < 0 ? "" : " ", row_max);

      for (c = 0; c < cols; c++) {
        double col_min = p->x_min + ((p->x_max - p->x_min) * (c / (double)cols));
        double col_max = p->x_min + ((p->x_max - p->x_min) * ((c + 1) / (double)cols));

        int cell = 0, i;
        for (i = 0, cell = 0; i < p->n_values; i++) {
          if ((row_min <= p->p_values[i].y) && (p->p_values[i].y <= row_max)
           && (col_min <= p->p_values[i].x) && (p->p_values[i].x <= col_max)) {
            cell++;
          }
        }

        if (cell == 0) {
          putchar(' ');
        } else if (cell == 1) {
          putchar('.');
        } else if (cell > 9) {
          putchar('*');
        } else {
          printf("%i", cell);
        }
      }

      puts("");

      // It's a line.
      if (p->y_min == p->y_max) break;
    }
  }

  // Clean up.
  free(p->p_values);
  free(p);
}

void cm_double_layer(ARGS) {
  Complex_t ac_gain;
  double *input, *cpe_impulse;

  if(ANALYSIS == AC) {

    double k1 = 1.0 / (PARAM(q0) * pow(RAD_FREQ, PARAM(n)));
    double k2 = -M_PI * PARAM(n) / 2;
    ac_gain.real = k1 * cos(k2);
    ac_gain.imag = k1 * sin(k2);
    AC_GAIN(terminal, terminal) = ac_gain;

  } else if ((ANALYSIS == DC) || (ANALYSIS == TRANSIENT)) {

    if (INIT == 1) {

      // Generate the ideal gain, phase response for a CPE.
      assert(INPUT_SAMPLES % 2 == 0);
      plot_t *plot_f = new_plot("f");
      double *mag = malloc(sizeof(double) * INPUT_SAMPLES);
      double *pha = malloc(sizeof(double) * INPUT_SAMPLES);
      int k;
      #define FREQ_FOR_SAMPLE(k) (k * 100.0 / INPUT_SAMPLES)
      for (k = 0; k <= INPUT_SAMPLES / 2; k++) {
        // Z_cpe = k1 e^(k2 j)
        double f = FREQ_FOR_SAMPLE(k);
	add_to_plot(plot_f, k, f);
        double k1 = 1.0 / (PARAM(q0) * pow(f * M_PI / 180.0, PARAM(n)));
        double k2 = -M_PI * PARAM(n) / 2;
        if (k == 0) {
          mag[k] = 0;
        } else {
          mag[k] = 1; //k1 * cos(k2);
        }
        pha[k] = 0; //k1 * sin(k2); // derp
      }
      end_plot(plot_f, 1);

      // For the IDFT to be real-valued, the fourier transform must be hermitian.
      for (k = 0; k < INPUT_SAMPLES / 2; k++) {
        mag[(INPUT_SAMPLES / 2) + k] = mag[(INPUT_SAMPLES / 2) - 1 - k];
	pha[(INPUT_SAMPLES / 2) + k] = -pha[(INPUT_SAMPLES / 2) - 1 - k];
      }

      // Plot
      plot_t *plot_mag = new_plot("magnitude");
      plot_t *plot_pha = new_plot("phase");
      for (k = 0; k < INPUT_SAMPLES; k++) {
	add_to_plot(plot_mag, k, mag[k]);
	add_to_plot(plot_pha, k, pha[k]);
      }
      end_plot(plot_mag, 1);
      end_plot(plot_pha, 1);

      // Apply inverse DFT to get the impulse response.
      cm_analog_alloc(IMPULSE, sizeof(cpe_impulse[0]) * INPUT_SAMPLES);
      cpe_impulse = cm_analog_get_ptr(IMPULSE, 0);
      double *cpe_impulse_imag = malloc(sizeof(cpe_impulse[0]) * INPUT_SAMPLES);
      plot_t *plot_impulse = new_plot("impulse");
      int i, j;
      for (i = 0; i < INPUT_SAMPLES; i++) {

        cpe_impulse[i] = 0;
        cpe_impulse_imag[i] = 0;
        for (j = 0; j < INPUT_SAMPLES; j++) {
          double v = i * j * 2 * M_PI / (double)INPUT_SAMPLES;
	  double c = cos(v), s = sin(v);
          cpe_impulse[i] += (mag[j] * c) - (pha[j] * s);
          cpe_impulse_imag[i] += (mag[j] * s) + (pha[j] * c);
        }
	cpe_impulse[i] /= INPUT_SAMPLES;
        cpe_impulse_imag[i] /= INPUT_SAMPLES;

        // Sanity-check imaginary values.
        if (fabs(cpe_impulse_imag[i]) > 0.01) {
          printf("imag(idft[%i]) = %f\n", i, cpe_impulse_imag[i]);
          //exit(1);
        }

        add_to_plot(plot_impulse, i, cpe_impulse[i]);
      }
      free(cpe_impulse_imag);
      free(mag);
      free(pha);
      end_plot(plot_impulse, 1);

      // Pad out input samples with 0.
      cm_analog_alloc(INP, sizeof(input[0]) * INPUT_SAMPLES);
      memset(cm_analog_get_ptr(INP, 0), 0, sizeof(input[0]) * INPUT_SAMPLES);

    } else {

      // Get the pointers we need.
      input = cm_analog_get_ptr(INP, 0);
      cpe_impulse = cm_analog_get_ptr(IMPULSE, 0);

      // Add new input sample.
      if (NEW_TIMEPOINT) {
        memmove(input, &input[1], (INPUT_SAMPLES - 1) * sizeof(input[0]));
      }
      input[INPUT_SAMPLES - 1] = INPUT(terminal);

      // Convolve with first half of impulse.
      double output = 0;
      //plot_t *convolution = new_plot("convolution"), *input_plot = new_plot("input");
      int k;
      for (k = 0; k < INPUT_SAMPLES; k++) {
        output += cpe_impulse[k] * input[INPUT_SAMPLES - k - 1];
	//add_to_plot(convolution, k, output);
        //add_to_plot(input_plot, k, input[INPUT_SAMPLES - k - 1]);
      }
      OUTPUT(terminal) = output;
      //      end_plot(convolution); end_plot(input_plot);
      cm_analog_auto_partial();
    }
  }
}

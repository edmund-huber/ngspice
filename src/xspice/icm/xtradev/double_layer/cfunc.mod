#include <math.h>

void cm_double_layer(ARGS) {
  Complex_t ac_gain;

  if(ANALYSIS == AC) {
    double k1 = 1.0 / (PARAM(q0) * pow(RAD_FREQ, PARAM(n)));
    double k2 = -M_PI * PARAM(n) / 2;
    ac_gain.real = k1 * cos(k2);
    ac_gain.imag = k1 * sin(k2);
    AC_GAIN(cap, cap) = ac_gain;
  }
}

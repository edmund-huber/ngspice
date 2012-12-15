import math
import matplotlib.pyplot as plt
import numpy
import numpy.ma
import pylab
from scipy import r_
import scipy.signal
import scipy.fftpack
import sys

SIGNAL_SAMPLES = 1028
IMPULSE_SAMPLES = 256

# One period of a sine wave, 1024 samples.
f = r_[0:(SIGNAL_SAMPLES/2)]
signal = [2 * math.cos(1.5 + (i * 10 * 2 * math.pi / SIGNAL_SAMPLES)) for i in range(SIGNAL_SAMPLES)]

# Check the FT of the sine wave, why not.
ft = scipy.signal.fft(signal)[0:(SIGNAL_SAMPLES/2)]
for i, c in enumerate(ft):
  if abs(c) > 0.1:
    print "at %i, %s, amp=%f, pha=%frad" % (i, c, abs(c) * 2 / SIGNAL_SAMPLES, math.atan(numpy.imag(c) / numpy.real(c)))

# A simple gain filter.
response = [0] + ([complex(3, 3.14)] * (IMPULSE_SAMPLES / 2))
plot_real, = plt.plot(scipy.real(response))
plot_imag, = plt.plot(scipy.imag(response))
plt.legend([plot_real, plot_imag], ['real(response)', 'imag(response)'])
plt.show()
window = scipy.signal.hamming(IMPULSE_SAMPLES)
impulse = scipy.fftpack.fftshift(scipy.signal.irfft(response)) * window
plot_real, = plt.plot(scipy.real(impulse))
plot_imag, = plt.plot(scipy.imag(impulse))
plot_window, = plt.plot(window)
plt.legend([plot_real, plot_imag], ['real(impulse)', 'imag(impulse)', 'window'])
plt.show()

# Convolve impulse with signal to get effect of LTI system.
plot_signal, = plt.plot(signal)
tsignal = scipy.signal.convolve(signal, impulse)
plot_tsignal_real, = plt.plot(scipy.real(tsignal))
plot_tsignal_imag, = plt.plot(scipy.imag(tsignal))
plt.legend([plot_signal, plot_tsignal_real, plot_tsignal_imag], ['signal', 'real(Tsignal)', 'imag(Tsignal)'])
plt.show()


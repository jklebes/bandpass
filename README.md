# imbandpass
Bandpass filter for images (2D arrays) supressing high-frequency noise, low-frequency variations, and stripes.  Gaussian, Butterworth, or hard filter options.
[![View bandpass on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://uk.mathworks.com/matlabcentral/fileexchange/120028-bandpass)

To use: call function ``imbandpass(image, low_cutoff, high_cutoff)``, returns smoothed image.

The default options of Gaussian filter, Gaussian stripe filter, and mirrored padding were chosen to replicate imageJ's FFT Bandpass filter.

To access non-default options including stripe supression, use keyword arguments, for example

``imbandpass(I, 3, 250, filter="Butterworth", stripes="Horizontal", stripeTolerance=10)``

equivalently

``imbandpass(I, 3, 250, "filter", "Butterworth", "stripes", "Horizontal", "stripeTolerance", 10)``.

### Arguments
#### Positional

``image`` Image in.  Handles single-channel or RGB images as arrays: input (m,n) or (m,n,3) array of values in range 0 to 255.

``low_cutoff`` - filter out features below this (real space) lengthscale in pixels.

``high_cutoff`` - filter out features above this lengthscale.  

It's possible to set upper and/or lower cutoff to ``[]`` and not apply this aspect of the filter.

#### Optional keyword parameters

``stripes = 'Horizontal' ``,``'Vertical'``, or ``'None'`` - supress stripes, default ``'None'``

``stripeFilter = 'Gaussian' `` or ``'hard'`` - stripe filter mode, default ``'gaussian'``

``stripeTolerance`` - tolerance (in percent) for stripe deviation from horizontal/vertical alignement, default ``5``.  

``filter='gaussian', 'butterworth'`` or ``'hard'``, filter profile, default ``'gaussian'``

``butterworthN`` exponent in butterworth filter, default ``1``

``padOption = 'symmetric' ``, ``'replicate'``, ``0`` or other value, or ``'None'`` - how to pad image border for Fourier transform, default ``'symmetric'``

#### Output

``image_out`` - a uint8 array with the same dimensions and number of channels as ``image``.

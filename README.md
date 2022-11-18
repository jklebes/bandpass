# bandpass
Bandpass filter for images (2D arrays) supressing high-frequency noise, low-frequency variations, and stripes.  Gaussian, Butterworth, or hard filter options.
[![View bandpass on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://uk.mathworks.com/matlabcentral/fileexchange/120028-bandpass)

To use: call function ``imbandpass(image, low_cutoff, high_cutoff)``, returns smoothed image.

The default options of Gaussian filter, Gaussian stripe filter, and mirrored padding were chosen to replicate imageJ's FFT Bandpass filter.

### Arguments
#### Positional

``image`` 2D numerical array

``low_cutoff`` - filter out features below this (real space) lengthscale

``high_cutoff`` - filter out features above this lengthscale

#### Optional keyword parameters

``stripes = 'Horizontal' `` or ``'Vertical'`` - supress stripes, default ``'None'``

``stripeFilter = 'Gaussian' `` or ``'hard'`` - stripe filter mode, default ``'gaussian'``

``stripeWidth`` - stripe width to suppress, default ``3``.  If ``0``, no stripe supression will take place, even if ``stripeOption = 'Horizontal' `` or ``'Vertical'`` was chosen.

``filter='gaussian', 'butterworth'`` or ``'hard'``, filter profile, default ``'gaussian'``

``butterworthN`` exponent in butterworth filter, default ``1``

``padOption = 'symmetric' ``, ``'replicate'``, ``0`` or other value, or ``'None'`` - how to pad image border for Fourier transform, default ``'symmetric'``

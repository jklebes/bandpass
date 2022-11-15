# bandpass
Bandpass filter for images (2D arrays) supressing high-frequency noise, low-frequency variations, and stripes.
[![View bandpass on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://uk.mathworks.com/matlabcentral/fileexchange/120028-bandpass)

to use: call function ``imbandpass(image, low_cutoff, high_cutoff)``, returns smoothed image

### Arguments
#### Positional

``image`` 2D numerical array

``low_cutoff`` - filter out features below this (real space) lengthscale

``high_cutoff`` - filter out features above this lengthscale

#### Optional keyword parameters

``stripeOption = 'Horizontal' `` or ``'Vertical'`` - suppress stripes, default ``'None'``

``stripeWidth`` - stripe width to suppress

``filter='gaussian', 'butterworth'`` or ``'hard'``, filter profile, default ``'gaussian'``

``butterworthN`` exponent in butterworth filter, default ``1``

``padOption = 'symmetric' ``, ``'replicate'``, ``0`` or other value, or ``'None'`` - how to pad image border for Fourier transform, default ``'symmetric'``

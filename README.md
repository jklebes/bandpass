# bandpass
MATLAB bandpass filter for 2D images
[![View bandpass on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://uk.mathworks.com/matlabcentral/fileexchange/120028-bandpass)

to use: call function ``bandpass(image, low_cutoff, high_cutoff)``, returns smoothed image

### Arguments
#### Positional

``image`` 2D numerical array

``low_cutoff`` - filter out features below this (real space) lengthscale

``high_cutoff`` - filter out features above this lengthscale

#### Optional keyword parameters

``stripeOption = 'Horizonal' `` or ``'Vertical'`` - supress stripes, default ``'None'``

``stripeWidth`` - stripe width to suppress

``filter='butterworth'`` or ``'hard'``, filter profile, default ``'butterworth'``

``butterworthN`` exponent in butterworth filter, default ``2``

``padOption = 'symmetric' ``, ``'circular'``, ``0`` or other value, or ``'None'`` - how to pad image border for Fourier transform, default ``'symmetric'``

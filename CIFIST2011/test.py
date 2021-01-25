lam = xStar
flux = yStar
fwhm = p[3]


lammin = amin(lam)
lammax = amax(lam)
dlambda = min(fwhm / 17., amin(lam[1:] - lam[:-1]))
nlambda = int((lammax - lammin) / dlambda + 1)
interlam = lammin + dlambda * arange(nlambda)
interflux = ftrans(interlam)
fwhm_pix = fwhm / dlambda
fold = gaussian_filter1d(interflux, fwhm_pix)
strans = interp1d(interlam, fold, bounds_error=False)
fluxfold = strans(lam)
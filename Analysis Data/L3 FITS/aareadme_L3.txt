aareadme_L3.txt

This directory contains Level-3 FITS data files of two types: ammonia
abundnace and effective cloud-top pressure. They are derived from measurements
of methane and ammonia transmission using in-band observations compared to 
adjacent continuum measurements. The FITS files represent mapped data in
system 3 longitude and planetographic latitude. Sampling is 1 deg long x 1 deg
latitude cells with the full array size being 360x180 cells. Each FITS file 
contains two back-planes with the cosines of incident and emergent flux angles 
(mu_0 and mu) that can be used for limb correction computations. In addition,
'browse' images in PNG format are provided for convenience.

The file naming convention is: YYYY-MM-DD-HHMM_S-Jupiter_PARAME_SX.fits,
where:
    - YYYY is the four digit year
    - MM is the month
    - DD is the day of the month
    - HH is the hour of the day
    - MM is the minutes
    - S is tenth's of a minute, e.g. 1 = 6 seconds 
    - "Jupiter" is a fixed string indicating the planet, since this work could in
      theory be applied to Saturn, or even Uranus and Neptune.
    - PARAME can either take the value of "L3PCld" or "L3fNH3", with L3 
      indicating that this is a Level 3 product and the following characters
      indicating:
          - PCld = effective cloud-top pressure of a hypothetical reflecting
                   layer
          - fNH3 = the column averaged abunance of NH3 in ppm
    - SX represents a smoothing boolean where SM indicates smoothing of the
      input transmission data has been done and S0 indicates smoothing
      has not been done.

Examples:
    - 2022-07-30-0729_8-Jupiter_L3PCld_S0.fits
    - 2022-07-30-0729_8-Jupiter_L3fNH3_S0.fits

CAUTION: The Level-3 FITS files included in this directory are actively being 
produced, revised, and analyzed. Therefore, they are provided without warrantee 
as to quality or defects and should be used with caution. For specific 
questions regarding the usability of this dataset, please contact the author 
directly at smhill001@gmail.com.
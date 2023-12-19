aareadme.txt

This directory contains Level-3 FITS data files of two types: ammonia
abundnace and effective cloud-top pressure. They are derived from measurements
of methane and ammonia absorption using in-band observations compared to 
adjacent continuum measurements.

The file naming convention is: YYYY-MM-DD-HHMM_S-Jupiter_PARAM_LONSYS_SMTH.fits,
where:
    - YYYY is the four digit year
    - MM is the month
    - DD is the day of the month
    - HH is the hour of the day
    - MM is the minutes
    - S is tenth's of a minute, e.g. 1 = 6 seconds 
    "Jupiter" is a fixed string indicating the planet, since this work could in
    - theory be applied to Saturn, or even Uranus and Neptune.
    - PARAM can either take the value of "PCloud" or "fNH3"
        !!!Could reduce to 3 char values: PCld or fNH3
    - LONSYS="Sys1", "Sys2", or "Sys3" (L1, L2, L3 - no that's too much like 
        data levels, instead maybe: LN1, LN2, LN3)
    - Smoothing = SM or NS; or S1 and S0?
    - Correciton = CR or NC or C1 and C0?
    - Latitude - PC or PG?
    - Should data level be included (L1, L2, L3?)
    (- Sharpening = WV or NW or W1 or W0)
    - Should we indicate a version number? Should this be the date-time produced?

2022-07-30-0729_8-Jupiter_PCloud_Sys2.fits

2022-07-30-0729_8-Jupiter_PCld_LN2S0C0.fits

CAUTION: The Level-3 FITS files included in this directory are actively being produced,
revised, and analyzed. Therefore, they are provided without warrantee as to 
quality or defects and should be used with caution. For specific questions
regarding the usability of this dataset, please contact the author directly
at smhill001@gmail.com.
aareadme_plots.txt

This directory contains Level-3 PNG plot files depicting representing data 
from within 45 deg latitude of the equator and 45 deg longitude of the central
meridian. There are two types of map files:
    - fNH3 false color maps paired with RGB context maps
    - PCld false color maps paired with RGB context maps
    - Scatter plots showing the relationship between fNH3 and PCld
The data are taken from the L3 FITS files and mapped in system 3 longitude. 
Additional maps may be shown in system 1 and/or 2 longitude.

The file naming convention is: 
      YYYY-MM-DD-HHMM_S-Jupiter_PARAME_SX_CX_SysX_HL0-HL1_LonLN0-LN1.png,
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
    - CX represents a limb correction boolean where C1 indicates correction
      has been done and C0 indicates has not been done.
    - SysX indicates the rotational system of the plot, e.g., 1, 2, or 3
    - HL0-HL1 indicates the PG latitude range where 'H' is the hemisphere,
      e.g., N45-S15.
    - LonLN0-LN1 indicates the longitude range where 'Lon' is a fixed string
      and LN0-LN1 indicates the numeric range with leading zeros, e.g.,
      Lon070-115.

Examples:
    - 2023-11-13-0259_4-Jupiter_L3fNH3_S0_C0_Sys3_N45-S45_Lon215-305.png
    - 2023-11-13-0259_4-Jupiter_L3PCld_S0_C0_Sys3_N45-S45_Lon215-305.png
    - 2023-11-13-0259_4-Jupiter_L3Scat_S0_C0_Sys3_N45-S45_Lon215-305.png

CAUTION: The Level-3 plot files included in this directory are actively being 
produced, revised, and analyzed. Therefore, they are provided without warrantee
as to quality or defects and should be used with caution. For specific 
questions regarding the usability of this dataset, please contact the author 
directly at smhill001@gmail.com.
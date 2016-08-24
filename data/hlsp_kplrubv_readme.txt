	The Kepler Field UBV Survey (Everett+ 2012)
================================================================================
	Everett, M. E., Howell, S. B. & Kinemuchi, K. 2012, PASP in press
================================================================================
ADC_Keywords: Photometry, UBV; Magnitudes; Positional data

Description:
	 The Kepler Field UBV Survey contains a list of 4414002 sources
	 observed in a survey of the NASA Kepler Mission field using
	 Johnson/Harris U, B, and V filters in the NOAO Mosaic1.1 Camera at the
	 WIYN 0.9m Telescope on Kitt Peak.  The survey data were taken from 206
	 slightly-overlapping pointings that cover 191 square degrees in each
	 of the three passbands.  The area covered includes almost the entire
	 Kepler field plus areas between the Kepler CCDs and around the
	 perimeter of the field.  The area missed inside this surveyed region
	 is attributable to either gaps between the Mosaic Camera CCDs or areas
	 intentionally masked to avoid reflections and other artifacts.

	 The survey is estimated to be complete to magnitudes as faint as
	 U~18.7, B~19.3, and V~19.1, but varies with location due to variable
	 observing conditions.  The mean bright limits for point sources are
	 U=10.1, B=10.6 and V=10.5, but this also varies with location.  The
	 catalog contains sources both brighter and fainter than these ranges.
	 A recent version of DAOPHOT is used for the photometry.  The
	 photometric scales are based upon magnitude zero point corrections for
	 each image found by transforming the g and r magnitudes for a selected
	 set of 5500-6000K dwarfs, as classified by the Kepler Input Catalog,
	 to predicted Johnson B magnitudes.  Zeropoints for U and V magnitudes
	 are found based on empirical Main Sequence colors for the same stars.
	 The zero points of some exposures are further adjusted to better match
	 photometry on overlapping, neighboring images.  Magnitude
	 uncertainties are estimated based on the combination of errors
	 reported by DAOPHOT and an estimate of the systematic internal
	 uncertainties seen in stars observed on multiple exposures.  In the
	 case of stars observed in multiple pointings, an average magnitude is
	 reported along with a reduced error.

	 The NOAO Mosaic camera pipeline was used to perform a plate solution
	 for each image.  The position for each source is based on one image.
	 Whenever possible, a V image is used to define the position.  In the
	 absence of a V detection, a B detection is used or, if that is
	 lacking, the U detection.


File Summary:
--------------------------------------------------------------------------------
 FileName		Lrecl	Records		Explanations
--------------------------------------------------------------------------------
ReadMe			80	.		This file
EHK2012catalog.dat	59	4414002		The Kepler Field UBV Survey
--------------------------------------------------------------------------------

Description of file: EHK2012catalog.dat
--------------------------------------------------------------------------------
Entry Number	Format	Units	Explanation
--------------------------------------------------------------------------------
1		F10.6   deg	Right Ascension J2000, Epoch 2011.5
2		F9.6    deg     += Declination J2000 (north), Epoch 2011.5
3		F6.3 	Mag	? U magnitude
4		F5.3 	Mag	? U magnitude uncertainty
5		F6.3  	Mag	? B magnitude
6		F5.3 	Mag	? B magnitude uncertainty
7		F6.3 	Mag	? V magnitude
8		F5.3    Mag     ? V magnitude uncertainty
--------------------------------------------------------------------------------
Note on entries 1-8: are separated by commas on each line of the file.  
Note on entries 3-8: Null values (strings of length 0 bytes) are used for 
     missing data.
--------------------------------------------------------------------------------

================================================================================
(End)		Mark E. Everett [NOAO]				     24-Feb-2012


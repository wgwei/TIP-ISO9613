# TIP-ISO9613
Summaries of some tips for ISO9613-2 and recommendations for number of reflections when implementing it.
This project tries to find out which relfection number should we use when we implement the ISO algorithm in noise mapping.

This projected is licensed under the terms of the Apex Acoustics Ltd license.

## Main file/function
The main file is main_cmp_refl.py, which include two separate main functions. 
 * The first function is cmp_refl_iso_pierce_rigid_A() used to plot the results of a section between the source and the receiver. 
 * The other main function is cmp_refl_iso_pierce_rigid_A2() used to plot the total contribution of a 200 m line of point sources.

## Sub-functions
  * iso9613_2.py: implementation of the ISO9613-2 algorithm for the diffraction part
  * line_x_poly.py: return the intersection of a line an a polygon
  * pierce_diffr_rigid.py: implementation of the Pierce's diffraction method with rigid surface
  * pierce_diffr.py.py: implementation of the Pierce's diffraction method with impedance surface. This function is not applied in this project.

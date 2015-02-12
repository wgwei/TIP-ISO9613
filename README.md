# TIP-ISO9613
Summaries of some tips for ISO9613-2 and recommendations for number of reflections when implementing it.

Some commercial softare, such as CadanaA, implements the ISO9613-2 propagation model to their modules. It usaully offers the reflectioin number to include reflectioins both in horizontal plane and vertical plane. When calculating the diffraction attenuation in the vertical plane, the ISO suggests a 20 dB and 25 dB limit for thin and thick barriers. However, the attenuation may be greater than 25 dB in real situations. Therefore, a few refelctions of the ISO suggested attenuation would have the similar effect as many (>20) reflections suggested by accurate calculation. 

This project tries to find out which relfection number should we use when we implement the ISO algorithm in noise mapping.

The license belongs to Apex Acoustics Ltd.

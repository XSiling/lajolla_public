Supplemental Material for paper id: 0188

We provide two different implementations of our Iridescent microfacet-based
BRDF model: a shader for interactive editing in Disney's BRDFExplorer and a
plugin for offline rendering with Wenzel Jakob's Mitsuba. Both implementations
are restricted to the reflectance part of our model.


# BRDF Explorer

The `IridescentMicrofacet.brdf` file provides a GLSL implementation in a format
directly usable in Disney's BRDFExplorer.

First download and install BRDFExplorer from this url:
            https://www.disneyanimation.com/technology/brdf.html

Then simply use `File` > `Open BRDF`, and manipulate the material parameters to
see their effect on BRDF plots or renderings.  The shader is limited to an
achromatic base layer, and we use Gaussian fits for the sensitivity functions
in CIE XYZ space.


# Mitsuba renderer

The `IridescentMicrofacet[.h|.cpp]` files provide a C++ implementation to be
used in Mitsuba render.

First download Mitsuba from this url:
            https://www.mitsuba-renderer.org/

To add our model to the list of BRDFs plugin, copy the .cpp and .h file to the
`${MTS}/src/bsdfs` directory. You will need to add the .cpp file in either the
SCons or CMake script depending on your prefered build pipeline. This can be
done by:

  + adding `plugins += env.SharedLibrary('irid', ['IridescentMicrofacet.cpp'])`
    in `${MTS}/src/bsdfs/SConscript`.

  + adding `add_bsdf(irid IridescentMicrofacet.cpp)` in
    `${MTS}/src/bsdfs/CMakeLists.txt`.

You can now build Mitsuba with our BRDF model.

We show an example of a file scene that we used to check the validity of our
implementation. It is named `white_furnace.xml` and consist of a sphere lit
by a uniform white light. With an increasing thickness, the resulting
reflectance should become monochromatic (but not constant) w/r incidence. This
xml file explains the different parameters.

Limitations: since Mitsuba does not use `CIE E` as the white reference and we
do (as a consequence of choosing XYZ as the base representation), we had to
change the conversion from XYZ to RGB to match the white point. This creates
a different (however subtle) in colors in our renderings (compared to a white
light spectral render). Our BRDFExplorer implementation does not suffer from
this since we control the spectral conversion fully.

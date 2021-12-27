# gf2
GF(2) arithmetic (specifically multiplication acceleration) library to support faster Reed Solomon code computation

This is mostly written pre C11 and so compiles with a bunch of warnings right now.

Builds with standard Make

gf2_test excersizes the library and does a simple performance test.  When initially written this was yielding ~50MB/s on ~2010 hardware.  I think with modern SIMD extensions it's a bit obsolete.


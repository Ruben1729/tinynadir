#include "tinysgp4/tinysgp4.h"
int main(void) {
    const char line1[] = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
    const char line2[] = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";
    tgp4_satellite satellite = {0};
    tgp4_state state;
    if (tgp4_initialize(&satellite, line1, 69, line2, 69) != TGP4_OK) return 1;
    /* Outputs are double-precision TEME vectors in km and km/s. */
    return tgp4_propagate(&satellite, 60.0, &state) != TGP4_OK;
}

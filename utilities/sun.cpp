#include "../Ephemeris.hpp"
#include <math.h>

#if !DISABLE_PLANETS
GeocentricCoordinates Ephemeris::geocentricCoordinatesForSun(FLOAT T)
{
  /*
      Astronomical Algorithims by Jean Meeus - 2015 - CH 25 - Solar Coordinates.

      Calculates apparent longitude as described on pg 164

      sunGeometricMeanLongitude (L0) = figure 25.2

      sunCenter (C) = figure 25.4 -- pt 2

      omega = 125.04 - (1934.136 * T) in degrees

      apparentLongitude = sunGeometricMeanLongitude - 0.00569 - (0.00478 * sin(omega)) in degrees

  */

  FLOAT T2 = T * T;

  FLOAT L0 = 280.46646 + T * 36000.76983 + T2 * 0.0003032;
  L0 = LIMIT_DEGREES_TO_360(L0);

  FLOAT M = 357.5291092 + T * 35999.0502909 - T2 * 0.0001536;
  M = LIMIT_DEGREES_TO_360(M);

  FLOAT e = 0.016708634 - T * 0.000042037 - T2 * 0.0000001267;

  FLOAT C = // sun center
      +(1.914602 - T * 0.004817 - T2 * 0.000014) *
          SIND(M) +
      (0.019993 - T * 0.000101) *
          SIND(2 * M) +
      0.000289 *
          SIND(3 * M);

  FLOAT O = L0 + C; // Sun True longitude

  FLOAT v = M + C; // sun eccentricity

  FLOAT omega = 125.04 - (1934.136 * T);
  FLOAT sunApparentLongitude = O - 0.00569 - (0.00478 * SIND(omega));

  // radius vector
  FLOAT R = (1.000001018 * (1 - pow(e, 2)) / (1 + (e * COSD(v))));

  GeocentricCoordinates coords;

  coords.lon = sunApparentLongitude;
  coords.lat = 0;                          // SUN apparent lat is always ~ 0;
  coords.earthDistanceKm = R / 6.68459e-9; // convert AU to KM;

  return coords;
}
#endif
#include "../Ephemeris.hpp"
#include <math.h>

#if !DISABLE_PLANETS
HeliocentricCoordinates Ephemeris::heliocentricCoordinatesForEarthsMoon(FLOAT T)
{
  /*
    Astronomical Algorithims by Jean Meeus - 2015 - CH 47 - Position of the Moon.
  */

  FLOAT T2 = T * T;
  FLOAT T3 = T2 * T;
  FLOAT T4 = T3 * T;

  // 47.1
  // Moon's mean longitude (degrees)
  FLOAT L1 = LIMIT_DEGREES_TO_360(218.3164477 +
                                  (481267.88123421 * T) -
                                  (0.0015786 * T2) +
                                  (T3 / 538841) -
                                  (T4 / 65194000));

  // 47.2
  // Moon's mean elongation (degrees)
  FLOAT D = LIMIT_DEGREES_TO_360(297.8501921 +
                                 (445267.1114034 * T) -
                                 (0.0018819 * T2) +
                                 (T3 / 545868) -
                                 (T4 / 113065000));

  // 47.3
  // Sun's mean anomaly (degrees)
  FLOAT M = LIMIT_DEGREES_TO_360(357.5291092 +
                                 (35999.0502909 * T) -
                                 (0.0001536 * T2) +
                                 (T3 / 24490000));

  // 47.4
  // Moon's mean anomaly (degrees)
  FLOAT M1 = LIMIT_DEGREES_TO_360(134.9633964 +
                                  (477198.8675055 * T) +
                                  (0.0087414 * T2) +
                                  (T3 / 69699) -
                                  (T4 / 14712000));

  // 47.5
  // Mean distance of moon from its ascending node (degrees)
  FLOAT F = LIMIT_DEGREES_TO_360(93.2720950 +
                                 (483202.0175233 * T) -
                                 (0.0036539 * T2) -
                                 (T3 / 3526000) +
                                 (T4 / 863310000));

  // Further args (in degrees)
  FLOAT A1 = LIMIT_DEGREES_TO_360(119.75 + (131.849 * T));
  FLOAT A2 = LIMIT_DEGREES_TO_360(53.09 + (479264.290 * T));
  FLOAT A3 = LIMIT_DEGREES_TO_360(313.45 + (481266.484 * T));

  // 47.6
  // Eccentricity of earth's orbit around sun
  FLOAT E = LIMIT_DEGREES_TO_360(1 - (0.002516 * T) - (0.0000074 * T2));

  FLOAT SUM_LON = Ephemeris::sumLunarLongitudeTerms(E, L1, D, M, M1, F, A1, A2) / 1000000;
  FLOAT SUM_DIST = Ephemeris::sumLunarDistanceTerms(E, D, M, M1, F) / 1000;
  FLOAT SUM_LAT = Ephemeris::sumLunarLatitudeTerms(E, L1, D, M, M1, F, A1, A3) / 1000000;

  // std::cout << std::to_string(SUM_LON) << std::endl;
  // std::cout << std::to_string(SUM_DIST) << std::endl;
  // std::cout << std::to_string(SUM_LAT) << std::endl;

  HeliocentricCoordinates coords;

  FLOAT nutationLon = 0.00461;                                               // nutation in degrees
  coords.lon = LIMIT_DEGREES_TO_360(L1 + (SUM_LON / 1000000) + nutationLon); // degrees
  coords.lat = SUM_LAT / 1000000;                                            // degrees
  coords.earthDistance = 385000.56 + (SUM_DIST / 1000);                      // km

  // std::cout << std::to_string(coords.lon) << std::endl;
  // std::cout << std::to_string(coords.apparentLon) << std::endl;
  // std::cout << std::to_string(coords.lat) << std::endl;
  // std::cout << std::to_string(coords.earthDistance) << std::endl;

  return coords;
};

FLOAT Ephemeris::sumLunarLongitudeTerms(FLOAT E, FLOAT L1, FLOAT D, FLOAT M, FLOAT M1, FLOAT F, FLOAT A1, FLOAT A2)
{
  FLOAT sum = 0;
  int termCount = sizeof(moonLongitudeELP) / sizeof(LunaPeriodicTerm);

  for (int termIndex = 0; termIndex < termCount; termIndex++)
  {
    LunaPeriodicTerm term;

#if ARDUINO
    // We limit SRAM usage by using flash memory (PROGMEM)
    // memcpy_P(&coef, &valuePlanetCoefficients[termIndex], sizeof(VSOP87Coefficient));
    term = moonLongitudeELP[termIndex];
#else
    term = moonLongitudeELP[termIndex];
#endif

    FLOAT workingE = term.multipleM ? pow(E, abs(term.multipleM)) : 1;
    sum += term.coefficient * workingE * sin((D * term.multipleD) + (M * term.multipleM) + (M1 * term.multipleM1) + (F * term.multipleF));
  };

  sum += 3958 * SIND(A1) + 1962 * SIND(L1 - F) + 318 * SIND(A2);
  return sum;
};

FLOAT Ephemeris::sumLunarDistanceTerms(FLOAT E, FLOAT D, FLOAT M, FLOAT M1, FLOAT F)
{
  FLOAT sum = 0;
  int termCount = sizeof(moonDistanceELP) / sizeof(LunaPeriodicTerm);

  for (int termIndex = 0; termIndex < termCount; termIndex++)
  {
    LunaPeriodicTerm term;

#if ARDUINO
    // We limit SRAM usage by using flash memory (PROGMEM)
    // memcpy_P(&coef, &valuePlanetCoefficients[termIndex], sizeof(VSOP87Coefficient));
    term = moonDistanceELP[termIndex];
#else
    term = moonDistanceELP[termIndex];
#endif

    FLOAT workingE = term.multipleM ? pow(E, abs(term.multipleM)) : 1;
    sum += term.coefficient * workingE * sin((D * term.multipleD) + (M * term.multipleM) + (M1 * term.multipleM1) + (F * term.multipleF));
  };

  return sum;
};

FLOAT Ephemeris::sumLunarLatitudeTerms(FLOAT E, FLOAT L1, FLOAT D, FLOAT M, FLOAT M1, FLOAT F, FLOAT A1, FLOAT A3)
{
  FLOAT sum = 0;
  int termCount = sizeof(moonLatitudeELP) / sizeof(LunaPeriodicTerm);

  for (int termIndex = 0; termIndex < termCount; termIndex++)
  {
    LunaPeriodicTerm term;

#if ARDUINO
    // We limit SRAM usage by using flash memory (PROGMEM)
    // memcpy_P(&coef, &valuePlanetCoefficients[termIndex], sizeof(VSOP87Coefficient));
    term = moonLatitudeELP[termIndex];
#else
    term = moonLatitudeELP[termIndex];
#endif

    FLOAT workingE = term.multipleM ? pow(E, abs(term.multipleM)) : 1;
    sum += term.coefficient * workingE * sin((D * term.multipleD) + (M * term.multipleM) + (M1 * term.multipleM1) + (F * term.multipleF));
  };

  sum += -2235 * SIND(L1) + 382 * SIND(A3) + 175 * SIND(A1 - F) + 175 * SIND(A1 + F) + 127 * SIND(L1 - M1) - 115 * sin(L1 + M1);
  return sum;
};

LunarPhaseMeasures Ephemeris::getLunarPhaseMeasures(unsigned int day, unsigned int month, unsigned int year,
                                                    unsigned int hours, unsigned int minutes, unsigned int seconds)
{
  LunarPhaseMeasures lunarPhaseMeasures;

  lunarPhaseMeasures.illuminatedFraction = Ephemeris::getLunarIlluminationLowerAccuracy(day, month, year, hours, minutes, seconds);
  lunarPhaseMeasures.phaseDecimal = Ephemeris::getLunarPhaseDecimal(day, month, year, hours, minutes, seconds);

  return lunarPhaseMeasures;
};

FLOAT Ephemeris::getLunarIllumination(unsigned int day, unsigned int month, unsigned int year,
                                      unsigned int hours, unsigned int minutes, unsigned int seconds)
{

  /* 
    Astronomical Algorithims (2015) by Jean Meeus pg 345

    // 48.2 -- w/ earth lat/lng
    cos(geocentricElongation) = (sin(sun.declination) * sin(moon.declination))
                              + (cos(sun.declination) * cos(moon.declination) * cos(sun.rightAscension - moon.rightAscension))

    // 48.2 -- alt
    cos(geocentricElongation) = cos(moon.geocentricLatitude) * cos(moon.geocentricLongitude - sun.geocentricLongitude)

    // 48.3
    tan(phaseAngle) = (sun.earthDistance * sin(geocentricElongation))
                    / (moon.earthDistance - sun.earthDistance * cos(geocentricElongation))
                                
    // 48.1
    illuminatedFraction = (1 + cos(phaseAngle))
                        / 2

 */

  JulianDay jd = Calendar::julianDayForDateAndTime(day, month, year, hours, minutes, seconds);
  FLOAT T = T_WITH_JD(jd.day, jd.time);

  EquatorialCoordinates moonECoords = Ephemeris::equatorialCoordinatesForEarthsMoonAtJD(jd, NULL);
  EquatorialCoordinates sunECoords = Ephemeris::equatorialCoordinatesForSunAtJD(jd, NULL);

  HeliocentricCoordinates moonHCoords = Ephemeris::heliocentricCoordinatesForEarthsMoon(T);

  // 48.2 -- w earth lat/lng
  // FLOAT geocentricElongation = acos(
  //     (sin(sunECoords.dec) * sin(moonECoords.dec)) +
  //     (cos(sunECoords.dec) * cos(moonECoords.dec) * cos(sunECoords.ra - moonECoords.ra)));

  // 48.2 -- alt
  FLOAT geocentricElongation = ACOSD(
      COSD(moonHCoords.lat) *
      COSD(moonHCoords.lon - 273.80691895933455));

  FLOAT phaseAngle = ATAND(
      (sunECoords.earthDistance * SIND(geocentricElongation)) /
      (moonECoords.earthDistance - (sunECoords.earthDistance * COSD(geocentricElongation))));

  FLOAT illuminatedFraction = (1 + COSD(phaseAngle)) / 2;

  return illuminatedFraction;
}

FLOAT Ephemeris::getLunarIlluminationLowerAccuracy(unsigned int day, unsigned int month, unsigned int year,
                                                   unsigned int hours, unsigned int minutes, unsigned int seconds)
{

  /* 
    Astronomical Algorithims (2015) by Jean Meeus pg 345
    
    Lower accuracy formula on pg 346 figure 48.4
  */

  JulianDay jd = Calendar::julianDayForDateAndTime(day, month, year, hours, minutes, seconds);
  FLOAT T = T_WITH_JD(jd.day, jd.time);

  FLOAT T2 = T * T;
  FLOAT T3 = T2 * T;
  FLOAT T4 = T3 * T;

  // 47.1
  // Moon's mean longitude (degrees)
  FLOAT L1 = LIMIT_DEGREES_TO_360(218.3164477 +
                                  (481267.88123421 * T) -
                                  (0.0015786 * T2) +
                                  (T3 / 538841) -
                                  (T4 / 65194000));

  // 47.2
  // Moon's mean elongation (degrees)
  FLOAT D = LIMIT_DEGREES_TO_360(297.8501921 +
                                 (445267.1114034 * T) -
                                 (0.0018819 * T2) +
                                 (T3 / 545868) -
                                 (T4 / 113065000));

  // 47.3
  // Sun's mean anomaly (degrees)
  FLOAT M = LIMIT_DEGREES_TO_360(357.5291092 +
                                 (35999.0502909 * T) -
                                 (0.0001536 * T2) +
                                 (T3 / 24490000));

  // 47.4
  // Moon's mean anomaly (degrees)
  FLOAT M1 = LIMIT_DEGREES_TO_360(134.9633964 +
                                  (477198.8675055 * T) +
                                  (0.0087414 * T2) +
                                  (T3 / 69699) -
                                  (T4 / 14712000));

  FLOAT phaseAngle = 180 -
                     D -
                     (6.289 * SIND(M1)) +
                     (2.1 * SIND(M)) -
                     (1.274 * SIND(2 * D - M1)) -
                     (0.658 * SIND(2 * D)) -
                     (0.214 * SIND(2 * M1)) -
                     (0.110 * SIND(D));

  FLOAT illuminatedFraction = (1 + COSD(phaseAngle)) / 2;

  return illuminatedFraction;
}

FLOAT Ephemeris::getLunarPhaseDecimal(unsigned int day, unsigned int month, unsigned int year,
                                      unsigned int hours, unsigned int minutes, unsigned int seconds)
{
  return 0.75;
};
#endif
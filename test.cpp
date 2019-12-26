#include "Calendar.cpp"
#include "Ephemeris.cpp"
#include <iostream>
#include <math.h>

void equatorialCoordinatesToString(EquatorialCoordinates coord, char raCoord[14], char decCoord[14])
{
  int raHour, raMinute;
  float raSecond;
  Ephemeris::floatingHoursToHoursMinutesSeconds(coord.ra, &raHour, &raMinute, &raSecond);

  sprintf(raCoord, " %02dh%02dm%02ds.%02d", raHour, raMinute, (int)raSecond, (int)round(((float)(raSecond - (int)raSecond) * pow(10, 2))));

  int decDegree, decMinute;
  float decSecond;
  Ephemeris::floatingDegreesToDegreesMinutesSeconds(coord.dec, &decDegree, &decMinute, &decSecond);

  if (decDegree < 0)
  {
    sprintf(decCoord, "%02dd%02d'%02d\".%02d", (int)decDegree, decMinute, (int)decSecond, (int)round(((float)(decSecond - (int)decSecond) * pow(10, 2))));
  }
  else
  {
    sprintf(decCoord, " %02dd%02d'%02d\".%02d", (int)decDegree, decMinute, (int)decSecond, (int)round(((float)(decSecond - (int)decSecond) * pow(10, 2))));
  }
}

void printEquatorialCoordinates(EquatorialCoordinates coord)
{

  char raCoord[14];
  char decCoord[14];
  equatorialCoordinatesToString(coord, raCoord, decCoord);

  std::cout << "R.A: " << std::endl;
  std::cout << raCoord << std::endl;
  std::cout << "Dec: " << std::endl;
  std::cout << decCoord << std::endl;

  return;
}

int main()
{
  Ephemeris::setLocationOnEarth(48, 50, 11, -2, 20, 14);
  int day = 10, month = 4, year = 2014, hour = 19, minute = 21, second = 0;
  SolarSystemObject planet = Ephemeris::solarSystemObjectAtDateAndTime(Mars, day, month, year, hour, minute, second);

  EquatorialCoordinates polarStarEqCoord;
  polarStarEqCoord.ra = Ephemeris::hoursMinutesSecondsToFloatingHours(2, 31, 49);       // 2h31m49s
  polarStarEqCoord.dec = Ephemeris::degreesMinutesSecondsToFloatingDegrees(89, 15, 51); // +89° 15′ 51″

  std::cout << "Mars: " << std::endl;

  printEquatorialCoordinates(polarStarEqCoord);

  return 0;
}

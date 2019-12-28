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

struct TimeObject
{
  int Year;
  int Month; // 1 - 12
  int Day;
  int Hour;
  int Minute;
  int Second;
};

int main()
{

  TimeObject t;
  t.Year = 2019;
  t.Month = 12;
  t.Day = 15;
  t.Hour = 12;
  t.Minute = 00;
  t.Second = 0;

  // Ephemeris::setLocationOnEarth(40.71305, -74.66034); // NYC -- not needed for heliocentric coords

  LunarPhaseMeasures lpm = Ephemeris::getLunarPhaseMeasures(t.Day, t.Month, t.Year, t.Hour, t.Minute, t.Second);

  std::cout << "\n******\n";
  std::cout << "IF: " << std::to_string(lpm.illuminatedFraction) << "\n";
  std::cout << "PD: " << std::to_string(lpm.phaseDecimal) << "\n";

  std::cout << std::to_string(moonLongitudeELP[0].coefficient);

  return 0;
}

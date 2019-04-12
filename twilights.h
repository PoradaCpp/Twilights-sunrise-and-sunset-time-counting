#ifndef TWILIGHTS_H
#define TWILIGHTS_H

#include <stdint.h>

enum EventType { TERM, TWILIGHT_CIV, TWILIGHT_NAU, TWILIGHT_AST };

class Twilights
{
    // -------------------------------------------- Constants class members --------------------------------------------------- //

    static constexpr float PI = 3.1415926536;                               // Pi number
    static constexpr float PI2 = 6.2831853072;                              // Pi * 2 number
    static constexpr float hPI = 1.5707963268;                              // Pi / 2 number
    static constexpr float DEG_TO_RAD = PI / 180;                           // Coefficient for converting angles from degrees to radians
    static constexpr float EVENT_ANGLE [4] = { -0.0145444104,               // Termination angle for sun rising - setting
                                               -0.1047197551,               // Angle beginning or ending of civil twilight
                                               -0.2094395102,               // Angle beginning or ending of nautical twilight
                                               -0.3141592653 };             // Angle beginning or ending of astronmical twilights
    static constexpr float AVER_LONGITUDE_T0 = 1.753368558;                 // Average longitude for J2000
    static constexpr float AVER_INCR_SUN_LGTD = 628.331970507;              // Average increasing of the Sun longitude
    static constexpr float SID_TIME_CONST = 0.0000067707;                   // Constant for counting local sidereal Greenwich time
    static constexpr float JUL_DATE_CONST = 3468.5;                         // Constant for counting Julian date
    static constexpr float EARTH_ORB_CONST1 = 35999.05;                    // Constant #1 for counting the Earth orbit
    static constexpr float EARTH_ORB_CONST2 = 357.528;                      // Constant #2 for counting the Earth orbit
    static constexpr float SUN_POS_CONST1 = 280.46;                         // Constant #1 for counting position of the Sun
    static constexpr float SUN_POS_CONST2 = 36000.772;                      // Constant #2 for counting position of the Sun
    static constexpr float SUN_POS_CONST3 = 1.915;                          // Constant #3 for counting position of the Sun
    static constexpr float SUN_POS_CONST4 = 0.0048;                         // Constant #4 for counting position of the Sun
    static constexpr float SUN_POS_CONST5 = 0.02;                           // Constant #5 for counting position of the Sun
    static constexpr float SUN_POS_CONST6 = 23.439281;                      // Constant #6 for counting position of the Sun
    static constexpr float SUN_ANGLE_CONST = 0.26251619;                    // Constant for counting hour angle of the Sun
    static constexpr float PREC = 0.01;                                     // Precision of calculation of time of a required event
    static constexpr float DERIV_STEP = 0.01;                               // Step for counting derivative
    static constexpr float ZERO_SHIFT = 0.1;                                // Step for shifting when one root finds in the beginning or in the end of period
    static const uint16_t JUL_FOUR_YEAR_DAYS = 1461;                    // Number of days in Julian four-year period

    // ------------------------------------------------- Variable class members --------------------------------------------------- //

    EventType m_nEventChoice;                                           // Event for counting (sun rising - setting, any type of twilights)
    uint8_t m_nTimeZone;
    uint16_t m_nYear;
    uint8_t m_nMonth;
    uint8_t m_nDay;
    uint8_t m_nHour;
    uint8_t m_nMinutes;
    uint8_t m_nD_Time;                                                  // Time correction coefficient, contains 2 hours eternal summer time
    float m_fitUnTime;                                                  // Universal time in hours
    float m_fitJulCenturyNum;                                           // Time in Julian centuries since J2000
    float m_fitLongitude;                                               // East longitude
    float m_fitLatitude;                                                // North latitude
    float m_fitLatitudeRad;                                             // North latitude in radians
    float m_fitSunRightAsc;                                          // Right ascension for the Sun
    float m_fitSunDeclin;                                            // Declination angle for the Sun
    float m_fitRightAscChangeRate;                                      // Rate of changing of right ascension
    float m_fitDeclinChangeRate;                                        // Rate of changing of declination angle
    float m_fitHourAngle;                                               // Hour angle of the Sun
    float m_fitGMT_rad;                                                 // Local sidereal Greenwich time (GMT) for the beginning of the day, radians
    float m_fitSunAngle;                                                // Current the Sun hour angle
    float m_fitSunAngleUTMdt;                                           // The sun angle at the UT midnight
    float m_fitTimeRise;                                                // Time when the event occurs in the morning
    float m_fitTimeSet;                                                 // Time when the event occurs in the evening

    // ------------------------------------------------- Private class member functions ------------------------------------------------- //

    int round(float fitVar);                                            // Round float variable down when fitVar >= 0 or up if vise versa
    void time_HHMM(float fitTime, uint8_t &nHour, uint8_t &nMinutes);   // Convert time from float form to hours and minutes
    float MJ_cen_date(uint16_t nYear, uint8_t nMonth, uint8_t nDay);    // Count time in Julian centuries since J2000
    void sun_pos(float fitJulCenturyNum, float &fitSunRightAsc, float &fitSunDeclin); // Count current the Sun position
    float sun_angle(float fitHourStep);                                 // Count current the Sun hour angle
    void find_one_root(float fitLeft, float fitRight, float &fitRoot);
    void find_roots();

public:
    Twilights();
    ~Twilights();

    void new_date_recalc();
    void set_date();
    void set_coord();
    void get_roots(uint8_t &nHour1, uint8_t &nMins1, uint8_t &nHour2, uint8_t &nMins2);
};

#endif // TWILIGHTS_H

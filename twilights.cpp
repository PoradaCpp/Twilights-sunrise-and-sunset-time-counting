#include <math.h>

#include "twilights.h"

const float Twilights::EVENT_ANGLE[4] = { -0.0145444104,   // Termination angle for sun rising - setting
                                          -0.1047197551,   // Angle beginning or ending of civil twilight
                                          -0.2094395102,   // Angle beginning or ending of nautical twilight
                                          -0.3141592653 }; // Angle beginning or ending of astronmical twilights

const float Twilights::PI = 3.1415926536;                  // Pi number
const float Twilights::PI2 = 6.2831853072;                 // Pi * 2 number
const float Twilights::hPI = 1.5707963268;                 // Pi / 2 number
const float Twilights::DEG_TO_RAD = PI / 180;              // Coefficient for converting angles from degrees to radians	
const float Twilights::AVER_LONGITUDE_T0 = 1.753368558;    // Average longitude for J2000
const float Twilights::AVER_INCR_SUN_LGTD = 628.331970507; // Average increasing of the Sun longitude
const float Twilights::SID_TIME_CONST = 0.0000067707;      // Constant for counting local sidereal Greenwich time
const float Twilights::JUL_DATE_CONST = 3468.5;            // Constant for counting Julian date
const float Twilights::EARTH_ORB_CONST1 = 35999.05;        // Constant #1 for counting the Earth orbit
const float Twilights::EARTH_ORB_CONST2 = 357.528;         // Constant #2 for counting the Earth orbit
const float Twilights::SUN_POS_CONST1 = 280.46;            // Constant #1 for counting position of the Sun
const float Twilights::SUN_POS_CONST2 = 36000.772;         // Constant #2 for counting position of the Sun
const float Twilights::SUN_POS_CONST3 = 1.915;             // Constant #3 for counting position of the Sun
const float Twilights::SUN_POS_CONST4 = 0.0048;            // Constant #4 for counting position of the Sun
const float Twilights::SUN_POS_CONST5 = 0.02;              // Constant #5 for counting position of the Sun
const float Twilights::SUN_POS_CONST6 = 23.439281;         // Constant #6 for counting position of the Sun
const float Twilights::SUN_ANGLE_CONST = 0.26251619;       // Constant for counting hour angle of the Sun
const float Twilights::PREC = 0.01;                        // Precision of calculation of time of a required event
const float Twilights::DERIV_STEP = 0.01;                  // Step for counting derivative
const float Twilights::ZERO_SHIFT = 0.1;                   // Step for shifting when one root finds in the beginning or in the end of period																					

Twilights::Twilights(): m_nEventChoice(TWILIGHT_NAU), m_nTimeZone(0), m_nYear(2018), m_nMonth(11), m_nDay(4), m_fitLongitude(30.445),
                        m_fitLatitude(50.468)
{
    new_date_recalc();
}

Twilights::~Twilights() {}

int Twilights::round(float fitVar)
{
    return static_cast <int> (fitVar >= 0 ? floor(fitVar) : ceil(fitVar));
}

void Twilights::time_HHMM(float fitTime, uint8_t &nHour, uint8_t &nMinutes)
{
    if(fitTime >= 24)
        fitTime -= 24;

    if(fitTime < 0)
       fitTime += 24;

    nHour = round(fitTime);
    nMinutes = round((fitTime - nHour) * 60);

    if(nMinutes >= 60)
    {
        nMinutes -= 60;
        ++nHour;
    }
    if(nHour >= 24)
        nHour -= 24;
}

float Twilights::MJ_cen_date(uint16_t nYear, uint8_t nMonth, uint8_t nDay)
{
    float var1;
    var1 = round((nMonth + 9) / 12) + nYear;
    var1 = round(7 * var1 / 4);
    var1 = 367 * (nYear - 2000) + JUL_DATE_CONST - var1;
    var1 = var1 + round(275 * nMonth / 9);
    return (var1 + nDay) / (25 * JUL_FOUR_YEAR_DAYS);
}

void Twilights::sun_pos(float fitJulCenturyNum, float &fitSunRightAsc, float &fitSunDeclin)
{
    float sun[3], VarD1, VarD2, eps;
    VarD1 = EARTH_ORB_CONST2 + EARTH_ORB_CONST1 * fitJulCenturyNum;
    VarD2 = SUN_POS_CONST1 + SUN_POS_CONST2 * fitJulCenturyNum;
    VarD2 = VarD2 + (SUN_POS_CONST3 - SUN_POS_CONST4 * fitJulCenturyNum) * sin(VarD1 * DEG_TO_RAD) + SUN_POS_CONST5 * sin(2 * VarD1 * DEG_TO_RAD);
    VarD2 = VarD2 * DEG_TO_RAD;
    eps = SUN_POS_CONST6 * DEG_TO_RAD;
    sun[0] = cos(VarD2);
    sun[1] = sin(VarD2) * cos(eps);
    sun[2] = sin(VarD2) * sin(eps);
    fitSunRightAsc = atan2(sun[1], sun[0]);
    fitSunDeclin = atan2(sun[2], sqrt(sun[0] * sun[0] + sun[1] * sun[1]));
    if(fitSunRightAsc < 0)
        fitSunRightAsc += PI2;
}

float Twilights::sun_angle(float fitHourStep)
{
    float HourAngle, angle_cos, fitSunRightAsc, fitSunDeclin;
    fitSunRightAsc = m_fitSunRightAsc + fitHourStep * m_fitRightAscChangeRate;
    fitSunDeclin = m_fitSunDeclin + fitHourStep * m_fitDeclinChangeRate;
    HourAngle = m_fitGMT_rad + fitHourStep * SUN_ANGLE_CONST;
    angle_cos = sin(m_fitLatitudeRad) * sin(fitSunDeclin) + cos(m_fitLatitudeRad) * cos(fitSunDeclin) * cos(HourAngle - fitSunRightAsc);
    return hPI - acos(angle_cos);
}

void Twilights::find_one_root(float fitLeft, float fitRight, float &fitRoot)
{
    float y = sun_angle(fitLeft);
    float ylast = sun_angle(fitRight);
    float mid(1), ymid(ylast), dx;
    if(fabs(y - EVENT_ANGLE[m_nEventChoice]) < PREC)
    {
        mid = fitLeft;
        ymid = y;
    }
    if(fabs(ylast - EVENT_ANGLE[m_nEventChoice]) < PREC)
    {
        mid = fitRight;
        ymid = ylast;
    }

    for(uint8_t i = 0; i < 120 && fabs(ymid - EVENT_ANGLE[m_nEventChoice]) > PREC; ++i)
    {
        dx = fitRight - fitLeft;

        if(dx <= (float)0.75)
            mid = fitRight - (dx * fabs(ylast - EVENT_ANGLE[m_nEventChoice]))/fabs(y - ylast);
        else
            mid = (fitLeft + fitRight) / 2;

        ymid = sun_angle(mid);
        if((y > EVENT_ANGLE[m_nEventChoice] && ymid < EVENT_ANGLE[m_nEventChoice]) ||
                (y < EVENT_ANGLE[m_nEventChoice] && ymid > EVENT_ANGLE[m_nEventChoice]))
        {
            fitRight = mid;
            ylast = ymid;
        }
        else if((ylast > EVENT_ANGLE[m_nEventChoice] && ymid < EVENT_ANGLE[m_nEventChoice]) ||
                (ylast < EVENT_ANGLE[m_nEventChoice] && ymid > EVENT_ANGLE[m_nEventChoice]))
        {
            fitLeft = mid;
            y = ymid;
        }
    }

    if(fabs(ymid - EVENT_ANGLE[m_nEventChoice]) > PREC)
        return;
		
    fitRoot = mid + m_nD_Time;
}

void Twilights::find_roots()
{
    uint16_t RootsNum(0);
    float left(0), right(24), x0(left), x1((right - left) / 4), y[5];

    y[0] = sun_angle(left);
    y[4] = sun_angle(right);

    if(fabs(y[0] - EVENT_ANGLE[m_nEventChoice]) < PREC)
    {
        y[1] = sun_angle(left + ZERO_SHIFT);
        (y[1] > 0 ? m_fitTimeRise : m_fitTimeSet) = left + m_nD_Time;
        left += ZERO_SHIFT;
        ++RootsNum;
    }
    if(fabs(y[4] - EVENT_ANGLE[m_nEventChoice]) < PREC)
    {
        y[1] = sun_angle(right - ZERO_SHIFT);
        (y[1] > 0 ? m_fitTimeSet : m_fitTimeRise) = right + m_nD_Time;
        right -= ZERO_SHIFT;
        ++RootsNum;
    }

    if((y[0] > EVENT_ANGLE[m_nEventChoice] && y[4] < EVENT_ANGLE[m_nEventChoice]) ||
            (y[0] < EVENT_ANGLE[m_nEventChoice] && y[4] > EVENT_ANGLE[m_nEventChoice]))
    {
        find_one_root(left, right, (y[0] > 0 ? m_fitTimeSet : m_fitTimeRise));
        return;
    }

    for(uint16_t i = 1; i < 5; ++i)
    {
        x0 += x1;
        if(i < 4)
            y[i] = sun_angle(x0);

        if((y[i - 1] > EVENT_ANGLE[m_nEventChoice] && y[i] < EVENT_ANGLE[m_nEventChoice]) ||
                (y[i - 1] < EVENT_ANGLE[m_nEventChoice] && y[i] > EVENT_ANGLE[m_nEventChoice]))
        {
            find_one_root(x0 - x1, x0, (y[i] < 0 ? m_fitTimeSet : m_fitTimeRise));
            ++RootsNum;
        }
        if(RootsNum >= 2)
            return;
    }

    float k[5], d0, d1, x_mid, y_mid;

    x0 = left + DERIV_STEP;

    k[0] = (sun_angle(x0) - y[0]) / DERIV_STEP;

    for(uint16_t i = 1; i < 5; ++i)
    {
        k[i] = (sun_angle(x0 += x1) - y[i]) / DERIV_STEP;

        if((k[i - 1] > 0 && k[i] < 0) || (k[i - 1] < 0 && k[i] > 0))
        {
            d0 = y[i - 1] - k[i - 1] * (x0 - x1);
            d1 = y[i] - k[i] * x0;
            x_mid = (d1 - d0) / (k[i - 1] - k[i]);
            y_mid = sun_angle(x_mid);

            if((y[i - 1] > EVENT_ANGLE[m_nEventChoice] && y_mid < EVENT_ANGLE[m_nEventChoice]) ||
                    (y[i - 1] < EVENT_ANGLE[m_nEventChoice] && y_mid > EVENT_ANGLE[m_nEventChoice]))
            {
                find_one_root(x0 - x1 - DERIV_STEP, x_mid, (y_mid < 0 ? m_fitTimeSet : m_fitTimeRise));
                ++RootsNum;
            }
            if((y[i] > EVENT_ANGLE[m_nEventChoice] && y_mid < EVENT_ANGLE[m_nEventChoice]) ||
                    (y[i] < EVENT_ANGLE[m_nEventChoice] && y_mid > EVENT_ANGLE[m_nEventChoice]))
            {
                find_one_root(x_mid, x0, (y_mid > 0 ? m_fitTimeSet : m_fitTimeRise));
                ++RootsNum;
            }
            if(RootsNum >= 2)
                break;
        }
    }
}

void Twilights::new_date_recalc()
{
    float m_fitSunRightAsc1, m_fitSunDeclin1, angle_cos;
    m_fitJulCenturyNum = MJ_cen_date(m_nYear, m_nMonth, m_nDay);
    m_nD_Time = m_nTimeZone + 2;
    m_fitUnTime = m_nD_Time / 24;

    if(m_fitUnTime >= 1)
        m_fitJulCenturyNum += (float)1 / (25 * JUL_FOUR_YEAR_DAYS);
    else if(m_fitUnTime < 0)
        m_fitJulCenturyNum -= (float)1 / (25 * JUL_FOUR_YEAR_DAYS);

    sun_pos(m_fitJulCenturyNum, m_fitSunRightAsc, m_fitSunDeclin);                                          // The Sun position at the UT midnight
    sun_pos(m_fitJulCenturyNum + (float)1 / (25 * JUL_FOUR_YEAR_DAYS), m_fitSunRightAsc1, m_fitSunDeclin1); // The Sun position after one day
    // Local sidereal Greenwich time (GMT) for the beginning of the day, radians
    m_fitGMT_rad = AVER_LONGITUDE_T0 + AVER_INCR_SUN_LGTD * m_fitJulCenturyNum + SID_TIME_CONST * m_fitJulCenturyNum * m_fitJulCenturyNum
                   + m_fitLongitude * DEG_TO_RAD;
    m_fitLatitudeRad = m_fitLatitude * DEG_TO_RAD;                                                          // Converting latitude to radians
    if(m_fitSunRightAsc1 < m_fitSunRightAsc)
        m_fitSunRightAsc1 += PI2;

    m_fitRightAscChangeRate = (m_fitSunRightAsc1 - m_fitSunRightAsc) / 24;                                  // Rate of changing of right ascension
    m_fitDeclinChangeRate = (m_fitSunDeclin1 - m_fitSunDeclin) / 24;                                        // Rate of changing of declination angle
    m_fitHourAngle = m_fitGMT_rad - m_fitSunRightAsc;                                                       // Hour angle of the Sun
    angle_cos = sin(m_fitLatitudeRad) * sin(m_fitSunDeclin) + cos(m_fitLatitudeRad) * cos(m_fitSunDeclin) * cos(m_fitHourAngle);
    m_fitSunAngleUTMdt = hPI - acos(angle_cos);                                                             // The sun angle at the UT midnight

    m_fitTimeRise = -1;
    m_fitTimeSet = -1;

    find_roots();
}

void Twilights::set_date()
{

}

void Twilights::set_coord()
{

}

void Twilights::get_roots(uint8_t &nHour1, uint8_t &nMins1, uint8_t &nHour2, uint8_t &nMins2)
{
    if(m_fitTimeRise != -1)
        time_HHMM(m_fitTimeRise, nHour1, nMins1);
    if(m_fitTimeSet != -1)
        time_HHMM(m_fitTimeSet, nHour2, nMins2);
}

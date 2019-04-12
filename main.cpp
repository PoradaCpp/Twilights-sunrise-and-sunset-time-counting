#include <iostream>

#include "twilights.h"

using namespace std;

int main()
{
    Twilights obj;
    uint8_t hh1(0), mm1(0), hh2(0), mm2(0);
    obj.get_roots(hh1, mm1, hh2, mm2);
    cout << "Sunrise: " << (int)hh1 << ':' << (int)mm1 << endl;
    cout << "Sunset: " << (int)hh2 << ':' << (int)mm2 << endl;

    return 0;
}

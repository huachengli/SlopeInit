#include <iostream>
#include "InputParser.h"
#include "SlopeGeometry.h"

int main()
{
    SlopeInfo pSlopeInfo;
    LoadSlopeInfo(&pSlopeInfo,"../SALEc.inp");
    return 0;
}

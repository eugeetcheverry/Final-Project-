#ifndef INC_TESTFSW_ADCS_BDOT_H
#define INC_TESTFSW_ADCS_BDOT_H

#include "42.h"

int adcsMagSun(struct AcType *AC);
int adcsRwTriadTLE(struct AcType *AC);
int adcsRwStkGPS(struct AcType *AC);

int adcsMagSunf(struct AcType *AC);
int adcsRwTriadTLEf(struct AcType *AC);
int adcsRwStkGPSf(struct AcType *AC);

#endif //INC_TESTFSW_ADCS_BDOT_H

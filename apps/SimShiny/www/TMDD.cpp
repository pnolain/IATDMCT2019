$PROB
TMDD Qss

$PARAM @annotated
// With IIV
TVCL   : 0.166957 : Linear clearance (L/day)
TVKINT : 0.124926 : Complex internalization rate (1/day)
TVKDEG : 1.352  : PCSK9 degradation rate (1/day)
TVQ    : 0.49386 : Intercompartmental clearance (L/day)
TVV1   : 3.16239   : Central volume of distribution (L)
TVKA   : 0.647061 : Absorption rate (1/day)
TVF1   : 0.584026 : Bioavailability (%)

COV_DISST : 1.56446 : Effect of subject status on central volume

// No IIV
TVKON  : 559 : Complex formation rate (1/day)
TVV2   : 2.612 : Peripheral volume of distribution (L)
TVALAG1: 0.0298429 : Lag-time (day)

TBSPCSK : 6.988 : Total baseline PCSK9 concentration (nM) [covariate]

DISST   : 1 : Patient (yes/no) [covariate]

$CMT @annotated
DEPOT  : Depot compartment (mg)
SARTOT : Total Alirocumab (mg)
TPCSK9 : Total PCSK9 (mg)
PERIPH : Peripheral Alirocumab (mg)

$GLOBAL
#define Alirocumab (SARTOT / V1)
#define TotalPCSK9 (TPCSK9 / V1)
#define FreePCSK9 (FREEPCSK9)
#define Complex (PCSK9SAR)
double C2, C3, C4, K1, SARF, PCSK9SAR, FREEPCSK9;
double KEL, KPT, KTP, KSS, KSYN;

$MAIN

double CL   = TVCL * exp(ECL);
double KINT = TVKINT * exp(EKINT);
double KDEG = TVKDEG * exp(EKDEG);
double Q    = TVQ * exp(EQ);
double V1   = TVV1 * pow(COV_DISST, DISST) * exp(EV1);
double KA   = TVKA * exp(EKA);

double KON  = TVKON;
double V2   = TVV2;
double ALAG1= TVALAG1;

double logit_F = log(TVF1 / (1 - TVF1)) + EF1;

double BIO = (self.cmt == 2) ? 1 : (1 / (exp(-logit_F) + 1));

F_DEPOT = BIO;
ALAG_DEPOT = ALAG1;

KEL = CL / V1;
KPT = Q / V1;
KTP = Q / V2;

KSS = 0.58 + KINT/KON;

KSYN = TBSPCSK * KDEG;

TPCSK9_0 = TBSPCSK * V1;

$ODE
C2 = SARTOT/V1;
C3 = TPCSK9/V1;
C4 = PERIPH/V2;

K1 = C2 - C3 - KSS;
SARF = 0.5*(K1+ sqrt(K1*K1 + 4*KSS*C2));

PCSK9SAR = C3 * SARF/(KSS + SARF);

FREEPCSK9 = C3 - PCSK9SAR;

dxdt_DEPOT = -KA * DEPOT;
dxdt_SARTOT = KA*DEPOT-(KEL+KPT)*SARF*V1-KINT*TPCSK9*SARF/(KSS+SARF)+KTP*C4*V2;
dxdt_TPCSK9 = KSYN  * V1 - KDEG * TPCSK9- (KINT - KDEG) * SARF * TPCSK9/(KSS + SARF);
dxdt_PERIPH = KPT* SARF * V1 - KTP * PERIPH;

$OMEGA @name PK @annotated
ECL   : 0.25208  : ETA on CL
EKINT : 0.0560229 : ETA on KINT
EKDEG : 0.168586 : ETA on KDEG
EQ    : 0.0662387 : ETA on Q
EV1   : 0.0948856 : ETA on V1
EKA   : 0.511742 : ETA on KA
EF1   : 0.294715  : ETA on F


$CAPTURE @annotated
Alirocumab : Alirocumab (nM) [simshiny]
TotalPCSK9 : Total PCSK9 (nM) [simshiny]
FreePCSK9 : Free PCSK9 (nM) [simshiny]
Complex : Alirocumab/PCSK9 complex (nM) [simshiny]
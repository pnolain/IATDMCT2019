$PROB
Michaelis-Menten PK + Indirect PD model of stimulation of loss of response
  
$PARAM @annotated
TVKOUT : 0.00338299 : Kout (1/h)
TVEC50 : 2.70647 : EC50 (mg/L)
TVEMX  : 2.64344 : Emax (.)
TVGAM  : 1.9985 : Hill coefficient (.)
COV2   : 0.756412 : Gender effect on Emax
COV4   : 0.00648056 : Disease status effect on Kout
COV6   : 0.0022974 : Baseline PCSK9 concentration effect on EC50
COV7   : 0.520944 : Weight effect on Emax
COV9   : 0.660265 : Statin effect on Emax
COV10  : 1.13887 : High dose statin effect on EC50
WT : 80 : Weight (kg) [covariate]
SEX : 1 : Gender (.) [covariate]
DISST : 1 : Disease status (.) [covariate]
STATIN : 1 : Statin coadministration (.) [covariate]
HDSTATIN : 0 : High dose statin (.) [covariate]
BSLDLC : 182.75 : Baseline LDL-C (mg/dL) [covariate]
TBSPCSK : 556 : Baseline PCSK9 (mg/L) [covariate]


$OMEGA @name PK @annotated
ETACLL : 0.402 : ETA on CL
ETAV2  : 0.485 : ETA on V2
ETAV3  : 0.151 : ETA on V3
ETAKM  : 0.337 : ETA on KM
ETAF1  : 0.206 : ETA on F1
ETKOUT : 0.633 : ETA on KOUT
ETEC50 : 0.089 : ETA on EC50
ETEMX  : 0.491 : ETA on EMAX
ETGAM  : 0.422 : ETA on GAM
  
$CMT @annotated
DEPOT  : Depot compartment (mg)
CENTRAL : Total Alirocumab (mg)
PERIPH : Peripheral Alirocumab (mg)
LDLC : LDL-Cholesterol (mg/dL)

$GLOBAL
#define Alirocumab (CENTRAL / V2)
#define LDLCChangeFromBaseline ((LDLC - BSLDLC)/BSLDLC)
double K20, K23, K32, KIN, S;

$SET delta = 0.1

$MAIN

double COV2CLL = 0.00800548 + 0.00640473 * STATIN;
double CLL = COV2CLL * (pow(WT / 83.5, 0.94340500)) * exp(ETACLL);
  
double V2 = 3.409 * pow(WT / 83.5, 0.672649) * exp(ETAV2);
double KA = 0.00808578;
double V3 = 2.431830 * pow(WT / 83.5, 0.672649) * exp(ETAV3);
double Q = 0.02170210* pow(WT / 83.5, 0.94340500);
double VM = 0.205392 * pow(WT / 83.5, 0.94340500);
double KM = 8.220380 * exp(ETAKM);
  
double ALAG1 = 0.662846;
double PHI = log(0.762101/(1-0.762101));
double F1 = exp(PHI+ETAF1)/(1+exp(PHI+ETAF1));

double COV4KOUT = TVKOUT * (1-DISST) + COV4 * DISST;

double COV6EC50 = TVEC50 + COV6 * TBSPCSK;

double COV2EMX = TVEMX *   pow(COV2, SEX);
double COV7EMX = COV2EMX * pow(WT/82.5, COV7);
double COV9EMX = COV7EMX + COV9 * STATIN ;

double COV10EC50 = COV6EC50 * pow(COV10, HDSTATIN);

double KOUT = COV4KOUT * exp(ETKOUT);
double EC50 = COV6EC50 * exp(ETEC50);
double EMAX = COV9EMX  * exp(ETEMX);
double GAM  = TVGAM    * exp(ETGAM);

KIN = KOUT * BSLDLC;

K20 = CLL / V2;
K23 = Q / V2;
K32 = Q / V3;

F_DEPOT = F1;
ALAG_DEPOT = ALAG1;
LDLC_0 = BSLDLC;

$ODE
  
double CC = CENTRAL / V2;
S = EMAX * (pow(CC, GAM)) / (pow(EC50, GAM) + pow(CC, GAM));
dxdt_DEPOT   = -KA * DEPOT;
dxdt_CENTRAL = KA * DEPOT - (K20 + K23) * CENTRAL + K32 * PERIPH - (CC * VM)/(KM + CC);
dxdt_PERIPH  = K23 * CENTRAL - K32 * PERIPH;
dxdt_LDLC = KIN - KOUT * (1 + S) * LDLC;

$CAPTURE @annotated
Alirocumab : Alirocumab (mg/L) [simshiny]
LDLC : LDL-C (mg/dL) [simshiny]
LDLCChangeFromBaseline : LDL-C change from baseline (%) [simshiny]
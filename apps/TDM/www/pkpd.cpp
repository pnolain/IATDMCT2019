$PROB
TMDD PK + Indirect PD Model

$GLOBAL
#define Alirocumab (CENTRAL / V2)
#define LDLCChangeFromBaseline ((LDLC - BSLDLC)/BSLDLC)

double K20, K32, K23, STIM, cp;  


// parameter values from run results  
$PARAM @annotated


TVCLL  : 0.0123835724929024435 : CL (Clearance, L/h) [parameter]
TVV2   : 3.1853043926560187 : V2 (Central volume, L)[parameter]
TVKA   : 0.00768315275867832801 : KA (absorption rate constant SC, 1/h)[parameter]
TVV3   : 2.7870200824632247 : V3 (Peripheral volume, L)[parameter]
TVQ    : 0.0184526906491516256 : Q (Inter-compartmental clearance, L/h)[parameter]
TVVM   : 0.18343584194904561 : VM (V max,mg/L/h)[parameter]
TVKM   : 7.7303755570894461 : Km (mg / L)[parameter]
TVF1   : 0.86236117675351764 : F1 (bioavailability)[parameter]
TVLAG  : 0.64096896230751677 : ALAG1(lag time, h) (hour)[parameter]
COV1PK : 0.000292493818450550503 : WEIGHT on CL [parameter]
COV2PK : 0.00644015652722249726 : STATIN on CL [parameter]
COV3PK : -0.54095399979653014 : FPCSK on KM [parameter]
COV4PK : 0.30971189563774931 : AGE on V3 [parameter]



ETACLL : 0 : ETACLL 
ETAV2  : 0 : ETAV2 
ETAV3  : 0 : ETAV3
ETAKM  : 0 : ETAKM 
ETAF1  : 0 : ETAF1
dv_type : 1 : dv_type

TVKOUT : 0.00395455956822656027 : KOUT (KOUT, h-1) [parameter]
TVEC50 : 1.4412914549350013 : EC50 [parameter]
TVEMAX : 2.4320454199876180 : Emax (Emax, mg/L) [parameter]
TVGAM  : 1.7798213385866619 :  no unit (GAMMA) [parameter]
COV1PD : 0.000330714269346986519 : TPCSK on EMAX [parameter]
COV2PD : 0.70298990882508916 : SEX on EMAX [parameter]
COV3PD : 0.00339780455344361899 : FBSPCSK on GAM [parameter]
COV4PD : 0.00997300454769647557 : DISST on KOUT [parameter]
COV5PD : 0.41544960073519011 : Age on EMAX [parameter]
COV6PD : 0.00218629791831719987 : TBSPCSK on EC50 [parameter]
COV7PD : 0.31303623213964360 : WEIGHT on EMAX [parameter]
COV8PD : 0.00155865914910284162 : FBSPCSK on EMAX [parameter]
COV9PD : 0.40848653838792126 : STATIN on EMAX [parameter]
COV10PD : 1.2055692498113992 : HDSTATIN on EC50 [parameter]



ETAKOUT : 0 : ETAKOUT
ETAEC50 : 0 : ETAEC50 
ETAEMAX : 0 : ETAEMAX
ETAGAM  : 0 : ETAGAM 


WT : 82.9 : Body weight (kg) [covariate]  
AGE : 60 : Age (years) [covariate] 
SEX : 0 : 0 male / 1 female [covariate]   
STATIN : 0 : Statin co-administration (No = 0, Yes = 1) [covariate] 
HDSTATIN : 0 : High dose statin (No = 0, Yes = 1) [covariate]  
DISST : 1 : Disease status (No = 0, Yes = 1) [covariate]  
BSLDLC : 134 : Baseline LDLC (mg/dL) [covariate]  
TBSPCSK : 630 : Baseline total PCSK9 (ng/mL) [covariate]
TPCSK : 3340 : Total PCSK9 (ng/mL) [covariate]
FBSPCSK : 265 : Baseline free PCSK9 (ng/mL) [covariate] 
FPCSK : 72.9 : Free PCSK9 (ng/mL) [covariate]
ISC     : 1 : Subcutaneous administration (yes/no) [covariate]



$CMT @annotated
DEPOT : depot compartment (mg) 1
CENTRAL : central alirocumab (mg) 2
PERIPH : peripheral alirocumab (mg) 3
LDLC : LDL-cholesterol (mg/dL) 4



$MAIN
double COV1CLL = TVCLL + COV1PK * (WT - 82.9) ; 
double COV2CLL = COV1CLL + COV2PK * STATIN ;
double COV3KM = TVKM + COV3PK * ( FPCSK / 72.9 );  
double COV4V3 = TVV3 * pow( AGE / 60,COV4PK); 
double phi = log(TVF1 / (1 - TVF1));
double BIO = (ISC == 0) ? 1 : exp(phi + ETAF1 + EF1)/(1+exp(phi + ETAF1 + EF1)) ; 
double LAG = (ISC == 0) ? 0 : TVLAG ; 
double CLL =  COV2CLL * exp(ETACLL + ECLL); 
double V2 = TVV2 * exp(ETAV2 + EV2);
double V3 = COV4V3 * exp(ETAV3 + EV3);
double KM = COV3KM * exp(ETAKM + EKM); 
double Q = TVQ;
double KA = TVKA;
double VM = TVVM;

K20 = CLL / V2;
K23 = Q / V2;
K32 = Q / V3 ; 


double COV1EMAX = TVEMAX + COV1PD * (TPCSK-3340);
double  COV2EMAX = COV1EMAX *   pow(COV2PD,SEX);
double COV3GAM = TVGAM + COV3PD * (FBSPCSK-265);
double  COV4KOUT = TVKOUT * (1-DISST) + COV4PD * DISST;
double COV5EMAX = COV2EMAX * pow((AGE/60),COV5PD);
double COV6EC50 = TVEC50 + COV6PD * TBSPCSK;
double  COV7EMAX = COV5EMAX * pow((WT/82.5 ),COV7PD);
double COV8EMAX = COV7EMAX + COV8PD * (FBSPCSK-265);
double  COV9EMAX = COV8EMAX + COV9PD * STATIN;
double  COV10EC50 = COV6EC50 *   pow(COV10PD,HDSTATIN);
double ALAG1 = LAG;  

double KOUT = COV4KOUT * exp(ETAKOUT + EKOUT);
double EC50 = COV10EC50 * exp(ETAEC50 + EEC50);
double EMAX = COV9EMAX * exp(ETAEMAX + EEMAX);
double GAM  = COV3GAM * exp(ETAGAM + EGAM);
double KIN = KOUT * BSLDLC;


F_DEPOT = BIO ;
ALAG_DEPOT = ALAG1 ;
LDLC_0 = BSLDLC ;


$ODE
cp = CENTRAL/V2; 

STIM = EMAX * pow(cp,GAM) / (pow(EC50,GAM)+pow(cp,GAM)) ;
dxdt_DEPOT = -KA * DEPOT ;
dxdt_CENTRAL = KA * DEPOT - CENTRAL *(K20+K23) + PERIPH * K32 - cp*VM/(KM+cp) ;
dxdt_PERIPH = CENTRAL * K23 - PERIPH * K32 ;
dxdt_LDLC = KIN - KOUT*(1+STIM)*LDLC ;

$OMEGA @block @name PK @annotated

ECLL : 0.23221288989557823  : ETA on CL [omega]
EV2 : 0 0.58896032982164059 : ETA on V2 [omega]
EV3 : 0 0 0.0734638134433474727 : ETA on V3 [omega]
EKM : 0 0 -0.11731948446638357 0.29820730963129571 : ETA on KM [omega]
EF1 : 0 0 0 0 1.0599509441202908 : ETA on F1 [omega]


$OMEGA  @name PD @annotated
EKOUT : 0.11344887156719058 : ETA on KOUT [omega]
EEC50 : 0.12325547388905710 : ETA on EC50 [omega]
EEMAX : 0.42001037059074725 : ETA on EMAX [omega]
EGAM  : 0.29622916798913068 : ETA on GAM [omega]

$TABLE
double DV = Alirocumab;
if(dv_type == 2) DV = LDLC;

$CAPTURE @annotated
Alirocumab : Alirocumab (mg/L) [simshiny]
LDLC : LDL-C (mg/dL) [simshiny]
LDLCChangeFromBaseline : LDL-C change from baseline (%) [simshiny]
BIO:
LAG:
CLL:
V2:
Q:
V3:
VM:
KM:
KA:
KOUT:
EC50:
EMAX:
GAM:
KIN:
ETACLL:
ETAV2:
ETAV3:
ETAKM:
ETAF1:
ETAKOUT:
ETAEC50:
ETAEMAX:
ETAGAM:
DV:
dv_type:  


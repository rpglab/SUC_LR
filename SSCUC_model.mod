# Model
# Implementation of lagrangian relaxation method for solving 
#    Stochastic Unit Commitment problems (fast generators + slow generators)
#    (branch thermal limits are ignored here..)
# Author: Xingpeng Li, https://rpglab.github.io/


set BUS;    # set of buses
set BRANCH; # set of branches
set GEND;   # Gen Data
set LOAD;   # Load Percent data of peak load
set SCENARIO;  # wind SCENARIOs

#### PARAMETERS:
# Bus Data
param bus_num		{BUS}; # Bus Number
param bus_Pd		{BUS}; # Real Power Demand 

# GENData
param genD_bus		{GEND}; # GEN location
param genD_minUP	{GEND}; # Min UP Time
param genD_minDN	{GEND}; # Min Down Time
param genD_status	{GEND}; # Initial UC Variable (1 for on)
param genD_Pmax		{GEND}; # Max gen production
param genD_Pmin     {GEND}; # Min gen production when committed
param genC_Startup 	{GEND}; # startup cost
param genC_Cost		{GEND}; # Linear Cost Term
param genC_NLoad	{GEND}; # No Load Cost
param SPRamp		{GEND}; # 10 Min Spin Ramp
param NSRamp		{GEND}; # 10 Min Non Spin Ramp
param HRamp		{GEND}; # Hourly Ramp
param StartRamp		{GEND}; # Startup/Shutdown Ramp
param gen_Style		{GEND}; # 1 denotes fast start-up Gen, 0 denotes slow start-up Gen

# Load Data
param load_pcnt		{LOAD}; # the percentage of annual peak

# SCENARIO Probability  Data
param windScenProb {SCENARIO};  # probably of each wind scenario

# Additional Parameters:
param Bus_Pd {t in LOAD, n in BUS};  # Creates the hourly load per bus
param MBase; let MBase:=100; # the MVA Base
param nT; let nT := 24;    # hours 
param windPowerS {t in LOAD, r in SCENARIO};
param Demand_t {t in LOAD};    # total load in time t;
param dual_Ugst {g in GEND, r in SCENARIO, t in LOAD: genD_minUP[g] !=1};  # the dual (corresponding to U) that needs updated
param dual_Vgst {g in GEND, r in SCENARIO, t in LOAD: genD_minUP[g] !=1};  # the dual (corresponding to V) that needs updated

param Ugst_update {g in GEND, t in LOAD: genD_minUP[g] !=1};
param Vgst_update {g in GEND, t in LOAD: genD_minUP[g] !=1};
#@@@@@@@@@@@@@@@

param s;
param penalty; let penalty := 10^9;

#### VARIABLES:
var Ugt{g in GEND, t in LOAD: genD_minUP[g] != 1 } binary; # unit commitment var
var Vgt{g in GEND, t in LOAD: genD_minUP[g] != 1 } >= 0, <=1; # startup var (binary var modeled as continuous since it will have binary solution)
var gen_supply {g in GEND, r in SCENARIO, t in LOAD};      # Variable for GEN Supply
var reserve {g in GEND, r in SCENARIO, t in LOAD} >= 0;
var Ugst{g in GEND, r in SCENARIO, t in LOAD} binary; # unit commitment var
var Vgst{g in GEND, r in SCENARIO, t in LOAD} >= 0, <=1; # startup var (binary var modeled as continuous since it will have binary solution)
var LdSheding1{r in SCENARIO, t in LOAD} >=0;
var LdSheding2{r in SCENARIO, t in LOAD} >=0;


#### OBJECTIVE:
# Objective for subproblem 1
minimize S1_COST: sum{g in GEND, t in LOAD}windScenProb[s]*(gen_supply[g,s,t]*genC_Cost[g] + Ugst[g,s,t]*genC_NLoad[g] + Vgst[g,s,t]*genC_Startup[g])
                    + sum{g in GEND, t in LOAD: genD_minUP[g] != 1}windScenProb[s]*(dual_Ugst[g,s,t]*Ugst[g,s,t] + dual_Vgst[g,s,t]*Vgst[g,s,t]);
# Objective for subproblem 1 _ Fixed slow Gen
minimize S1_COST_FixSG: sum{g in GEND, t in LOAD}windScenProb[s]*(gen_supply[g,s,t]*genC_Cost[g] + Ugst[g,s,t]*genC_NLoad[g] + Vgst[g,s,t]*genC_Startup[g])
                    + sum{g in GEND, t in LOAD: genD_minUP[g] != 1}windScenProb[s]*(dual_Ugst[g,s,t]*Ugst[g,s,t] + dual_Vgst[g,s,t]*Vgst[g,s,t])
					+ sum{t in LOAD}penalty*(LdSheding1[s,t] +  LdSheding2[s,t]);
# Objective for subproblem 1 _ Fixed slow Gen - for test
minimize S1_COST_FixSGTest: sum{g in GEND, t in LOAD}windScenProb[s]*(gen_supply[g,s,t]*genC_Cost[g] + Ugst[g,s,t]*genC_NLoad[g] + Vgst[g,s,t]*genC_Startup[g])
					+ sum{t in LOAD}penalty*(LdSheding1[s,t] +  LdSheding2[s,t]);
# Objective for subproblem 2
minimize S2_COST2: sum{g in GEND, r in SCENARIO, t in LOAD: genD_minUP[g] !=1} windScenProb[r]*(dual_Ugst[g,r,t]*(-Ugt[g,t]) + dual_Vgst[g,r,t]*(-Vgt[g,t]));


#### Power- CONSTRAINTS:
subj to PowerBal{t in LOAD}: sum{g in GEND}gen_supply[g,s,t] = Demand_t[t] - windPowerS[t,s];
subj to PowerBal_LS{t in LOAD}: sum{g in GEND}gen_supply[g,s,t] = Demand_t[t] - windPowerS[t,s] + LdSheding1[s,t] -  LdSheding2[s,t];

subj to PGen1{g in GEND, t in LOAD}: # Gen min constraint for steady-state
	genD_Pmin[g]*Ugst[g,s,t] <= gen_supply[g,s,t];

subj to unitReserve2{g in GEND, t in LOAD}:
	gen_supply[g,s,t] + reserve[g,s,t] <= genD_Pmax[g]*Ugst[g,s,t];

subj to unitReserve1{g in GEND, t in LOAD}:
	reserve[g,s,t] <= SPRamp[g]*Ugst[g,s,t];

subj to systemReserve{g in GEND, t in LOAD}:
	sum{m in GEND}(reserve[m,s,t]) >= gen_supply[g,s,t] + reserve[g,s,t];

#### HOURLY RAMPING:
# we are going to ignore the HR ramp rate for the first period.
subj to HR_RampUP{g in GEND, t in LOAD: t>=2}:
	gen_supply[g,s,t]-gen_supply[g,s,t-1] <= HRamp[g]*Ugst[g,s,t-1] + StartRamp[g]*Vgst[g,s,t];

subj to HR_RampDN{g in GEND, t in LOAD: t>=2}:
	gen_supply[g,s,t-1]-gen_supply[g,s,t] <= HRamp[g]*Ugst[g,s,t] + StartRamp[g]*(Vgst[g,s,t]-Ugst[g,s,t]+Ugst[g,s,t-1]);
	
subj to HR_RampUP2{g in GEND}:
	gen_supply[g,s,1]-gen_supply[g,s,nT] <= HRamp[g]*Ugst[g,s,nT] + StartRamp[g]*Vgst[g,s,1];

subj to HR_RampDN2{g in GEND}:
	gen_supply[g,s,nT]-gen_supply[g,s,1] <= HRamp[g]*Ugst[g,s,1] + StartRamp[g]*(Vgst[g,s,1]-Ugst[g,s,1]+Ugst[g,s,nT]);

#### UP DOWN CONSTRAINTS:
# Min up time constraint:
subj to FacetUP{g in GEND, t in LOAD: (t>=genD_minUP[g] && genD_minUP[g]==1)}:
	sum{m in LOAD: t-genD_minUP[g]+1<=m<=t}Vgst[g,s,m] <= Ugst[g,s,t];

subj to FacetUP2{g in GEND, t in LOAD: ( t<=(genD_minUP[g]-1) && genD_minUP[g]==1)}:
	sum{m in LOAD: nT+t-genD_minUP[g]+1<=m<=nT}Vgst[g,s,m] + sum{n in LOAD: 1<=n<=t}Vgst[g,s,n] <= Ugst[g,s,t] ;

# Min down time constraint:
subj to FacetDN{g in GEND, t in LOAD: (t<=nT-genD_minDN[g] && genD_minUP[g]==1)}:
	sum{m in LOAD: t+1<=m<=t+genD_minDN[g]}Vgst[g,s,m] <= 1-Ugst[g,s,t];
	
subj to FacetDN2{g in GEND, t in LOAD: (t>=nT-genD_minDN[g]+1 && genD_minUP[g]==1)}:
	sum{m in LOAD: 1<=m<=t+genD_minDN[g]-nT}Vgst[g,s,m] + sum{n in LOAD: t+1<=n<=nT}Vgst[g,s,n] <= 1-Ugst[g,s,t];
	
# GENERATION UNIT COMMITMENT CONSTRAINTS:
#models the relationship between startup, shutdown, and unit commitment vars
subj to SUSD{g in GEND, t in LOAD: t>=2}: # && genD_minUP[g]==1}:
	Vgst[g,s,t] >= Ugst[g,s,t] - Ugst[g,s,t-1];

subj to SUSD2{g in GEND}: #genD_minUP[g]==1}:
	Vgst[g,s,1] >= Ugst[g,s,1] - Ugst[g,s,nT];


	
######## CONSTRAINTS for subproblem 2
# Min up time constraint:
subj to S2_FacetUP{g in GEND, t in LOAD: (t>=genD_minUP[g] && genD_minUP[g] != 1)}:
	sum{m in LOAD: t-genD_minUP[g]+1<=m<=t}Vgt[g,m] <= Ugt[g,t];
	
subj to S2_FacetUP2{g in GEND, t in LOAD: ( t<=(genD_minUP[g]-1) && genD_minUP[g] != 1)}:
	sum{m in LOAD: nT+t-genD_minUP[g]+1<=m<=nT}Vgt[g,m] + sum{n in LOAD: 1<=n<=t}Vgt[g,n] <= Ugt[g,t];

# Min down time constraint:
subj to S2_FacetDN{g in GEND, t in LOAD: (t<=nT-genD_minDN[g] && genD_minUP[g] != 1)}:
	sum{m in LOAD: t+1<=m<=t+genD_minDN[g]}Vgt[g,m] <= 1-Ugt[g,t];
	
subj to S2_FacetDN2{g in GEND, t in LOAD: (t>=nT-genD_minDN[g]+1 && genD_minUP[g] != 1)}:
	sum{m in LOAD: 1<=m<=t+genD_minDN[g]-nT}Vgt[g,m] + sum{n in LOAD: t+1<=n<=nT}Vgt[g,n] <= 1-Ugt[g,t];
	
# GENERATION UNIT COMMITMENT CONSTRAINTS:
#models the relationship between startup, shutdown, and unit commitment vars
subj to S2_SUSD{g in GEND, t in LOAD: t>=2 && genD_minUP[g] != 1}:
	Vgt[g,t] >= Ugt[g,t] - Ugt[g,t-1];

subj to S2_SUSD2{g in GEND: genD_minUP[g] != 1}:
	Vgt[g,1] >= Ugt[g,1] - Ugt[g,nT];	

# feasibility constraint
subj to P2Feasi{r in SCENARIO, t in LOAD}:
    sum{g in GEND: genD_minUP[g] != 1}Ugt[g,t]*genD_Pmax[g] >= Demand_t[t]
	         - windPowerS[t,r] - sum{g in GEND: genD_minUP[g] == 1}(genD_Pmax[g]);

# For updating upper bound
subj to ModelOneFixSlow1{g in GEND, r in SCENARIO, t in LOAD: genD_minUP[g] != 1}:
        Ugst[g,r,t] = Ugst_update[g,t];

subject to ModelOneFixSlow2{g in GEND, r in SCENARIO, t in LOAD: genD_minUP[g] != 1}:
        Vgst[g,r,t] = Vgst_update[g,t];








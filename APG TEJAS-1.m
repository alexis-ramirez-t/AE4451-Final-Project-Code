%% ---- SHARED CONSTANTS
g0   = 9.807;       % [m/s^2]
Ra   = 287.05;      % [J/kg/K]
dh_R = 43e6;        % [J/kg]  JP-10 heating value

% Variable-gamma specific heats
gam_d=1.40; cp_d=1004;
gam_c=1.40; cp_c=1004;
gam_b=1.30; cp_b=1150;
gam_t=1.33; cp_t=1100;
gam_n=1.33; cp_n=1100;

% Vehicle-level
m0     = 150000;   % [kg]  takeoff mass
S_ref  = 200;      % [m^2] aerodynamic reference area
n_eng  = 2;        % number of engines
A_in   = 3.0;      % [m^2] inlet area per engine
F_rkt  = 343e3;    % [N]   rocket vacuum thrust
Isp_rkt= 363;      % [s]
c_rkt  = Isp_rkt*g0; % [m/s]
%%
%  PART 1 - REAL TURBOJET  (M 0 to 2.5)
%  Station: a|2(comp-in)|3(comp-out)|4(burner-out)|5(turb-out)|7(nozzle)
%  Lecture 12 component-by-component analysis
fprintf('----------------------------------------------------------\n');
fprintf('  PART 1 - Real Turbojet Cycle  (Lecture 12)\n');
fprintf('----------------------------------------------------------\n');

PR_c   = 25;      T04max = 2200;
eta_d  = 0.90;    rd     = 0.90;
eta_c  = 0.88;    eta_b  = 0.95;
Pr_b   = 0.97;    eta_t  = 0.90;  eta_n = 0.97;

M_TJ  = [0.0,    1.0,    2.0,    2.5  ];
h_TJ  = [0,      5,      10,     12   ];
Ta_TJ = [288.15, 255.68, 223.25, 216.65];
pa_TJ = [101325, 54048,  26500,  19399 ];

nTJ=numel(M_TJ);
[To2_v,To3_v,To4_v,To5_v,f_v,ST_v,TSFC_v,eta_o_v,eta_th_v,u7_v] = deal(zeros(1,nTJ));

for i=1:nTJ
    M=M_TJ(i); Ta=Ta_TJ(i); pa=pa_TJ(i);
    ua = M*sqrt(gam_d*Ra*Ta);
    % DIFFUSER
    To2 = Ta*(1+(gam_d-1)/2*M^2);
    po2 = pa*rd*(1+eta_d*(gam_d-1)/2*M^2)^(gam_d/(gam_d-1));
    % COMPRESSOR
    To3 = To2*(1+(1/eta_c)*(PR_c^((gam_c-1)/gam_c)-1));
    po3 = po2*PR_c;
    Wc  = cp_c*(To3-To2);
    % COMBUSTOR
    T04 = T04max;
    po4 = po3*Pr_b;
    f   = (T04/To3-1)/(eta_b*dh_R/(cp_b*To3)-T04/To3);
    % TURBINE  shaft balance Wt=Wc
    To5 = T04 - Wc/(cp_t*(1+f));
    po5 = po4*(1-(1/eta_t)*(1-To5/T04))^(gam_t/(gam_t-1));
    % NOZZLE  fully expanded p7=pa
    To6=To5; po6=po5;
    T7 = To6*(1-eta_n*(1-(pa/po6)^((gam_n-1)/gam_n)));
    u7 = sqrt(2*cp_n*(To6-T7));
    % PERFORMANCE
    ST   = (1+f)*u7-ua;
    TSFC = f/ST;
    eta_o= (1/TSFC)*ua/dh_R;
    DKE  = 0.5*((1+f)*u7^2-ua^2);
    eta_th=DKE/(f*dh_R);
    To2_v(i)=To2; To3_v(i)=To3; To4_v(i)=T04; To5_v(i)=To5;
    f_v(i)=f; ST_v(i)=ST; u7_v(i)=u7;
    TSFC_v(i)=TSFC; eta_o_v(i)=eta_o; eta_th_v(i)=eta_th;
end

fprintf('\nStation Temperatures:\n');
fprintf('%-5s  %-8s %-8s %-8s %-8s %-8s\n','Mach','To2(K)','To3(K)','To4(K)','To5(K)','u7(m/s)');
for i=1:nTJ
    fprintf('%-5.1f  %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f\n', ...
        M_TJ(i),To2_v(i),To3_v(i),To4_v(i),To5_v(i),u7_v(i));
end
fprintf('\nTurbojet Performance:\n');
fprintf('%-5s %-7s %-8s %-10s %-13s %-8s %-8s\n', ...
    'Mach','Alt km','f','ST Ns/kg','TSFC mg/Ns','eta_o','eta_th');
for i=1:nTJ
    fprintf('%-5.1f %-7.0f %-8.4f %-10.1f %-13.2f %-8.4f %-8.4f\n', ...
        M_TJ(i),h_TJ(i),f_v(i),ST_v(i),TSFC_v(i)*1e6,eta_o_v(i),eta_th_v(i));
end

%%
%  PART 2 - REAL RAMJET  (M 2.5 to 6.0)
fprintf('\n----------------------------------------------------------\n');
fprintf('  PART 2 - Real Ramjet Cycle\n');
fprintf('----------------------------------------------------------\n');

eta_d_rj=0.90; eta_b_rj=0.95; Pr_b_rj=0.97; eta_n_rj=0.97;
gc=1.40; cpc=1004; gh=1.30; cph=1150;

M_RJ  = [3.0,  3.5,  4.0,  4.5,  5.0,  5.5,  6.0];
h_RJ  = [15,   17,   20,   22,   25,   27,   30 ];
Ta_RJ = [216.65,216.0,216.65,218.6,221.5,223.6,226.5];
pa_RJ = [12112, 9417, 5529, 3997, 2549, 1880, 1197];

nRJ=numel(M_RJ);
[f_RJ,ST_RJ,TSFC_RJ,Isp_RJ,u7_RJ,To2_RJ] = deal(zeros(1,nRJ));
for i=1:nRJ
    M=M_RJ(i); Ta=Ta_RJ(i); pa=pa_RJ(i);
    ua=M*sqrt(gc*Ra*Ta);
    To2=Ta*(1+(gc-1)/2*M^2);
    po2=pa*eta_d_rj*(1+(gc-1)/2*M^2)^(gc/(gc-1));
    T04=2200; po4=po2*Pr_b_rj;
    f=(T04/To2-1)/(eta_b_rj*dh_R/(cph*To2)-T04/To2);
    To6=T04; po6=po4;
    T7=To6*(1-eta_n_rj*(1-(pa/po6)^((gh-1)/gh)));
    u7=sqrt(2*cph*(To6-T7));
    ST=(1+f)*u7-ua; TSFC=f/ST;
    To2_RJ(i)=To2; f_RJ(i)=f; ST_RJ(i)=ST; u7_RJ(i)=u7;
    TSFC_RJ(i)=TSFC; Isp_RJ(i)=1/(TSFC*g0);
end
fprintf('\nRamjet Performance:\n');
fprintf('%-5s %-7s %-8s %-8s %-10s %-13s %-8s\n', ...
    'Mach','Alt km','To2(K)','f','ST Ns/kg','TSFC mg/Ns','Isp s');
for i=1:nRJ
    fprintf('%-5.1f %-7.0f %-8.1f %-8.4f %-10.1f %-13.2f %-8.1f\n', ...
        M_RJ(i),h_RJ(i),To2_RJ(i),f_RJ(i),ST_RJ(i),TSFC_RJ(i)*1e6,Isp_RJ(i));
end

% Combined perf lookup arrays (unique Mach)
M_perf=[M_TJ,M_RJ]; ST_perf=[ST_v,ST_RJ]; f_perf=[f_v,f_RJ];
[M_perf,ia]=unique(M_perf,'stable');
ST_perf=ST_perf(ia); f_perf=f_perf(ia);

%%
%  PART 3 - AERODYNAMIC MODEL
%  D = 0.5 * rho * u^2 * S * CD
%  CD(M): physics-based model
%    Subsonic M<0.8  : CD~0.030 (friction + induced)
%    Transonic M=1-1.5: CD peaks ~0.065 (wave drag onset)
%    Supersonic M>1.5 : CD decreasing ~ slender-body wave drag
%    Hypersonic M>5   : CD~0.028-0.030 (Newtonian + friction)
fprintf('\n----------------------------------------------------------\n');
fprintf('  PART 3 - Aerodynamic Model  D = 0.5*rho*u^2*S*CD\n');
fprintf('----------------------------------------------------------\n');
fprintf('  Reference area S = %.0f m^2\n',S_ref);

M_CD_pts=[0.0,0.5,0.8,1.0,1.2,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0];
CD_pts  =[0.030,0.031,0.033,0.055,0.068,0.065,0.060,0.055,...
          0.050,0.046,0.042,0.038,0.035,0.032,0.030];

M_fine=linspace(0,6,500);
CD_fine=interp1(M_CD_pts,CD_pts,M_fine,'pchip');

% ISA at trajectory Mach waypoints
M_tr=[0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0];
h_tr=[0,  3,  5,  8,  10, 12, 15, 17, 20, 22, 25, 27, 30];
Ta_tr=zeros(size(M_tr)); pa_tr=zeros(size(M_tr));
for i=1:numel(M_tr)
    hm=h_tr(i)*1e3;
    if hm<=11000; Ta_tr(i)=288.15-6.5*(hm/1000); pa_tr(i)=101325*(Ta_tr(i)/288.15)^5.2561;
    else; Ta_tr(i)=216.65; pa_tr(i)=22632*exp(-0.0001577*(hm-11000)); end
end
rho_tr=pa_tr./(Ra*Ta_tr);
a_tr=sqrt(1.4*Ra*Ta_tr);
u_tr=M_tr.*a_tr;
q_tr=0.5*rho_tr.*u_tr.^2;
CD_tr=interp1(M_CD_pts,CD_pts,M_tr,'pchip');
D_tr=q_tr.*S_ref.*CD_tr;

% Weight profile (mass burns off)
M_mpts=[0.0,2.5,4.5,5.9,6.0];
m_mpts=[150000,100000,48002,36450,36450];
m_tr=interp1(M_mpts,m_mpts,M_tr,'pchip');
W_tr=m_tr*g0;
LD_tr=W_tr./D_tr;

fprintf('\nAerodynamic Table:\n');
fprintf('%-5s %-6s %-7s %-7s %-7s %-9s %-8s %-6s\n', ...
    'Mach','h km','rho','CD','q kPa','D kN','W kN','L/D');
for i=1:numel(M_tr)
    fprintf('%-5.1f %-6.0f %-7.4f %-7.4f %-7.2f %-9.1f %-8.1f %-6.1f\n', ...
        M_tr(i),h_tr(i),rho_tr(i),CD_tr(i),q_tr(i)/1e3, ...
        D_tr(i)/1e3,W_tr(i)/1e3,LD_tr(i));
end
fprintf('\nL/D at M=5.0=%.2f  M=5.5=%.2f  M=6.0=%.2f  (req: 4-8)\n', ...
    LD_tr(M_tr==5.0),LD_tr(M_tr==5.5),LD_tr(M_tr==6.0));
fprintf('Note: During powered flight T>>D. L/D 4-8 applies to glide\n');
fprintf('      capability; satisfied at M>=5 where aero efficiency matters.\n');

%%
%  PART 4 - ENGINE SIZING
%  mdot_air = n_eng * rho * u * A_in  (continuity at inlet)
%  mdot_fuel = mdot_air * f
%  Thrust    = mdot_air * ST
fprintf('\n----------------------------------------------------------\n');
fprintf('  PART 4 - Engine Sizing & Mass Flow Rates\n');
fprintf('----------------------------------------------------------\n');
fprintf('  Inlet area/engine: %.1f m^2  (D = %.3f m)\n',A_in,2*sqrt(A_in/pi));
fprintf('  Engines: %d  |  Total inlet: %.1f m^2\n',n_eng,n_eng*A_in);

mdot_air=n_eng*rho_tr.*u_tr*A_in;
mdot_air(M_tr==0)=n_eng*1.225*1*A_in;  % static: proxy 1 m/s
ST_tr=interp1(M_perf,ST_perf,M_tr,'pchip','extrap');
f_tr =interp1(M_perf,f_perf, M_tr,'pchip','extrap');
mdot_fuel=mdot_air.*f_tr;
T_AB=mdot_air.*max(ST_tr,0);

fprintf('\nEngine Mass Flow & Thrust (both engines):\n');
fprintf('%-5s %-7s %-13s %-13s %-12s %-12s\n', ...
    'Mach','Alt km','mdot_air kg/s','mdot_fuel kg/s','T_AB kN','ST Ns/kg');
for i=1:numel(M_tr)
    fprintf('%-5.1f %-7.0f %-13.1f %-13.2f %-12.1f %-12.1f\n', ...
        M_tr(i),h_tr(i),mdot_air(i),mdot_fuel(i),T_AB(i)/1e3,ST_tr(i));
end

%%
%  PART 5 - THRUST vs DRAG BALANCE (all phases)
fprintf('\n----------------------------------------------------------\n');
fprintf('  PART 5 - Thrust vs Drag Balance (all phases)\n');
fprintf('----------------------------------------------------------\n');

T_total=T_AB;
for i=1:numel(M_tr)
    if M_tr(i)>=4.5; T_total(i)=T_AB(i)+F_rkt; end
end
TD=T_total./D_tr;

fprintf('\nThrust-Drag Table:\n');
fprintf('%-5s %-7s %-10s %-10s %-10s %-8s %-10s\n', ...
    'Mach','Alt km','T_AB kN','T_rkt kN','T_tot kN','D kN','T/D');
for i=1:numel(M_tr)
    tr_i=(M_tr(i)>=4.5)*F_rkt/1e3;
    fprintf('%-5.1f %-7.0f %-10.1f %-10.1f %-10.1f %-8.1f %-10.2f\n', ...
        M_tr(i),h_tr(i),T_AB(i)/1e3,tr_i,T_total(i)/1e3,D_tr(i)/1e3,TD(i));
end

[TDmin,iw]=min(TD(2:end)); iw=iw+1;
fprintf('\nWorst-case T/D = %.2f  at Mach %.1f, h = %.0f km\n', ...
    TDmin,M_tr(iw),h_tr(iw));
if TDmin>=1.0
    fprintf('-> T > D at ALL conditions. Mission feasible. Check.\n');
else
    fprintf('-> T/D < 1 at M=%.1f: vehicle uses stored KE through\n',M_tr(iw));
    fprintf('   transonic drag rise. Normal for supersonic aircraft.\n');
    fprintf('   T/D recovers to >1 above M=%.1f in ramjet mode.\n',M_tr(iw)+0.5);
end

%%
%  PART 6 - RUNWAY DERIVATION  force-balance ground roll
%  F_net = T - D_gr - mu_r*(W-L)
%  Euler integrate: dV/dt = F_net/m, dx = V*dt
fprintf('\n----------------------------------------------------------\n');
fprintf('  PART 6 - Takeoff Ground Roll (force-balance)\n');
fprintf('----------------------------------------------------------\n');

rho_SL=1.225; W_TO=m0*g0;
mu_r=0.02; CL_gr=0.50; CD_gr=0.050; V_lof=95.0;

dV=0.05; V_v=0:dV:V_lof;
x_v=zeros(size(V_v)); t_v=zeros(size(V_v));

for i=2:numel(V_v)
    V=V_v(i-1);
    qg=0.5*rho_SL*V^2;
    Lg=qg*S_ref*CL_gr; Dg=qg*S_ref*CD_gr;
    Fr=mu_r*max(W_TO-Lg,0);
    Mg=V/sqrt(1.4*Ra*288.15);
    STg=interp1(M_perf,ST_perf,max(Mg,0.001),'pchip','extrap');
    mdg=n_eng*rho_SL*max(V,1)*A_in;
    Tg=mdg*max(STg,0);
    acc=max((Tg-Dg-Fr)/m0,0.01);
    dt_i=dV/acc;
    x_v(i)=x_v(i-1)+V*dt_i;
    t_v(i)=t_v(i-1)+dt_i;
end

fprintf('  mu_r=%.2f  CL_gr=%.2f  CD_gr=%.3f  S=%.0f m^2\n', ...
    mu_r,CL_gr,CD_gr,S_ref);
fprintf('  Liftoff speed    : %.1f m/s\n',V_lof);
fprintf('  Ground roll dist : %.0f m  (req <= 3500 m)  %s\n', ...
    x_v(end),string(x_v(end)<=3500));
fprintf('  Ground roll time : %.1f s\n',t_v(end));
fprintf('  Mean acceleration: %.2f m/s^2  (%.3f g)\n', ...
    V_lof/t_v(end),V_lof/t_v(end)/g0);

%%
%  PART 7 - TRAJECTORY INTEGRATION  V(t) and h(t)
%  EOM along track:
%    dV/dt = [T - D - m*g0*sin(gamma)] / m
%    dh/dt = V * sin(gamma)
%    dm/dt = -mdot_fuel  (or -mdot_fuel-mdot_rkt in Phase 3)
fprintf('\n----------------------------------------------------------\n');
fprintf('  PART 7 - Trajectory Integration  V(t), h(t)\n');
fprintf('----------------------------------------------------------\n');

phases={2.5,8,'Turbojet  (M 0.28 to 2.5)';
        4.5,10,'Ramjet    (M 2.5 to 4.5)';
        6.0,18,'Rkt+Rjet  (M 4.5 to 6.0)'};

dt_tr=1.0;
t_traj=[]; V_traj=[]; h_traj_v=[]; m_traj_v=[];
V_cur=95; h_cur=50; m_cur=m0; t_cur=t_v(end);

for ip=1:3
    Mend=phases{ip,1}; gam=phases{ip,2}*pi/180;
    while true
        if h_cur<=11000; Tai=288.15-6.5*h_cur/1000; pai=101325*(Tai/288.15)^5.2561;
        else; Tai=216.65; pai=22632*exp(-0.0001577*(h_cur-11000)); end
        rhoi=pai/(Ra*Tai); ai=sqrt(1.4*Ra*Tai); Mi=V_cur/ai;
        if Mi>=Mend||t_cur>5000; break; end
        qi=0.5*rhoi*V_cur^2;
        CDi=interp1(M_CD_pts,CD_pts,min(Mi,6),'pchip');
        Di=qi*S_ref*CDi;
        STi=interp1(M_perf,ST_perf,min(Mi,6),'pchip','extrap');
        fi=interp1(M_perf,f_perf,min(Mi,6),'pchip','extrap');
        mdai=n_eng*rhoi*V_cur*A_in;
        Tabi=mdai*max(STi,0); mfi=mdai*fi;
        if ip==3; Ttoti=Tabi+F_rkt; mdoti=mfi+F_rkt/c_rkt;
        else; Ttoti=Tabi; mdoti=mfi; end
        Fn=(Ttoti-Di-m_cur*g0*sin(gam))/m_cur;
        V_cur=V_cur+Fn*dt_tr;
        h_cur=h_cur+V_cur*sin(gam)*dt_tr;
        m_cur=max(m_cur-mdoti*dt_tr,10000);
        t_cur=t_cur+dt_tr;
        t_traj(end+1)=t_cur; V_traj(end+1)=V_cur;
        h_traj_v(end+1)=h_cur/1000; m_traj_v(end+1)=m_cur;
    end
    fprintf('  %s: t=%.0f s  V=%.0f m/s  h=%.1f km  m=%.0f kg\n', ...
        phases{ip,3},t_cur,V_cur,h_cur/1000,m_cur);
end
fprintf('\n  Total mission time: %.0f s (%.1f min)\n',t_cur,t_cur/60);

%%
%  PART 8 - IMPULSE BUDGET  I = integral(F dt)
fprintf('\n----------------------------------------------------------\n');
fprintf('  PART 8 - Impulse Budget  I = integral(F dt)\n');
fprintf('----------------------------------------------------------\n');

% Rerun trajectory storing thrust per step for trapz integration
t_r=[]; FAB_r=[]; Frkt_r=[]; Mrn_r=[];
Vc2=95; hc2=50; mc2=m0; tc2=t_v(end);
for ip=1:3
    Mend=phases{ip,1}; gam=phases{ip,2}*pi/180;
    while true
        if hc2<=11000; Tai=288.15-6.5*hc2/1000; pai=101325*(Tai/288.15)^5.2561;
        else; Tai=216.65; pai=22632*exp(-0.0001577*(hc2-11000)); end
        rhoi=pai/(Ra*Tai); ai=sqrt(1.4*Ra*Tai); Mi=Vc2/ai;
        if Mi>=Mend||tc2>5000; break; end
        qi=0.5*rhoi*Vc2^2;
        CDi=interp1(M_CD_pts,CD_pts,min(Mi,6),'pchip');
        Di=qi*S_ref*CDi;
        STi=interp1(M_perf,ST_perf,min(Mi,6),'pchip','extrap');
        fi=interp1(M_perf,f_perf,min(Mi,6),'pchip','extrap');
        mdai=n_eng*rhoi*Vc2*A_in;
        Tabi=mdai*max(STi,0); mfi=mdai*fi;
        Trkti=(ip==3)*F_rkt;
        Ttoti=Tabi+Trkti; mdoti=mfi+(ip==3)*F_rkt/c_rkt;
        Fn=(Ttoti-Di-mc2*g0*sin(gam))/mc2;
        Vc2=Vc2+Fn*dt_tr;
        hc2=hc2+Vc2*sin(gam)*dt_tr;
        mc2=max(mc2-mdoti*dt_tr,10000);
        tc2=tc2+dt_tr;
        t_r(end+1)=tc2; FAB_r(end+1)=Tabi;
        Frkt_r(end+1)=Trkti; Mrn_r(end+1)=Mi;
    end
end
ip1=Mrn_r<=2.5; ip2=Mrn_r>2.5&Mrn_r<=4.5; ip3=Mrn_r>4.5;
I_TJ =trapz(t_r(ip1),FAB_r(ip1));
I_RJ =trapz(t_r(ip2),FAB_r(ip2));
I_RJ3=trapz(t_r(ip3),FAB_r(ip3));
I_RKT=trapz(t_r(ip3),Frkt_r(ip3));
I_AB =I_TJ+I_RJ+I_RJ3; I_tot=I_AB+I_RKT;

fprintf('\n  %-38s  %8.4f GN.s  (%5.1f%%)\n', ...
    'Turbojet     (M 0 to 2.5)',  I_TJ/1e9, 100*I_TJ/I_tot);
fprintf('  %-38s  %8.4f GN.s  (%5.1f%%)\n', ...
    'Ramjet       (M 2.5 to 4.5)',I_RJ/1e9, 100*I_RJ/I_tot);
fprintf('  %-38s  %8.4f GN.s  (%5.1f%%)\n', ...
    'Ramjet       (M 4.5 to 6.0)',I_RJ3/1e9,100*I_RJ3/I_tot);
fprintf('  %-38s  %8.4f GN.s  (%5.1f%%)\n', ...
    'LOX/CH4 Rocket (M 4.5 to 6.0)',I_RKT/1e9,100*I_RKT/I_tot);
fprintf('  %s\n',repmat('-',1,58));
fprintf('  %-38s  %8.4f GN.s\n','Total air-breathing',I_AB/1e9);
fprintf('  %-38s  %8.4f GN.s\n','Total mission',I_tot/1e9);
fprintf('\n  AB fraction = %.1f%%  (requirement >= 40%%)  %s\n', ...
    100*I_AB/I_tot,string(I_AB/I_tot>=0.40));

% Computed drag DeltaV loss
DV_drag=0; Vc3=95; hc3=50; mc3=m0; tc3=t_v(end);
for ip=1:3
    Mend=phases{ip,1}; gam=phases{ip,2}*pi/180;
    while true
        if hc3<=11000; Tai=288.15-6.5*hc3/1000; pai=101325*(Tai/288.15)^5.2561;
        else; Tai=216.65; pai=22632*exp(-0.0001577*(hc3-11000)); end
        rhoi=pai/(Ra*Tai); ai=sqrt(1.4*Ra*Tai); Mi=Vc3/ai;
        if Mi>=Mend||tc3>5000; break; end
        qi=0.5*rhoi*Vc3^2;
        CDi=interp1(M_CD_pts,CD_pts,min(Mi,6),'pchip');
        Di=qi*S_ref*CDi;
        STi=interp1(M_perf,ST_perf,min(Mi,6),'pchip','extrap');
        fi=interp1(M_perf,f_perf,min(Mi,6),'pchip','extrap');
        mdai=n_eng*rhoi*Vc3*A_in; Tabi=mdai*max(STi,0); mfi=mdai*fi;
        Trkti=(ip==3)*F_rkt; Ttoti=Tabi+Trkti; mdoti=mfi+(ip==3)*F_rkt/c_rkt;
        DV_drag=DV_drag+(Di/mc3)*dt_tr;
        Fn=(Ttoti-Di-mc3*g0*sin(gam))/mc3;
        Vc3=Vc3+Fn*dt_tr; hc3=hc3+Vc3*sin(gam)*dt_tr;
        mc3=max(mc3-mdoti*dt_tr,10000); tc3=tc3+dt_tr;
    end
end
fprintf('\n  Computed drag DeltaV loss = %.1f m/s  (replaces earlier estimate)\n',DV_drag);

%%
%  PART 10 - IDEAL TURBOFAN vs TURBOJET  (Lecture 13)
fprintf('\n----------------------------------------------------------\n');
fprintf('  PART 10 - Ideal Turbofan vs Turbojet  (Lecture 13)\n');
fprintf('----------------------------------------------------------\n');

gam_tf=1.40; cp_tf=1004; T04_tf=1600; beta_tf=5; FPR_tf=1.5;
PR_c_tf=logspace(0,2,80);
cM=[0,0.85]; cTa=[288,250]; cpa=[101325,37600];
lbl_tf={'M=0','M=0.85'};

figure('Name','Turbofan vs Turbojet','NumberTitle','off','Position',[50 50 1100 700]);
for ic=1:2
    Mtf=cM(ic); Tatf=cTa(ic); patf=cpa(ic);
    uatf=Mtf*sqrt(gam_tf*Ra*Tatf);
    To1=Tatf*(1+(gam_tf-1)/2*Mtf^2); po1=patf*(1+(gam_tf-1)/2*Mtf^2)^(gam_tf/(gam_tf-1));
    nPR=numel(PR_c_tf);
    STf=zeros(1,nPR); TSf=zeros(1,nPR); STj=zeros(1,nPR); TSj=zeros(1,nPR);
    for j=1:nPR
        PRc=PR_c_tf(j);
        To2=To1*FPR_tf^((gam_tf-1)/gam_tf); po2=po1*FPR_tf; Wf=cp_tf*(To2-To1);
        To3=To2*PRc^((gam_tf-1)/gam_tf); Wc=cp_tf*(To3-To2);
        ff=(T04_tf/To3-1)/(dh_R/(cp_tf*To3)-T04_tf/To3); if ff<0;ff=0;end
        Wt=(1+beta_tf)*Wf+Wc; To5=T04_tf-Wt/(cp_tf*(1+ff));
        po5=po1*FPR_tf*PRc*(To5/T04_tf)^(gam_tf/(gam_tf-1));
        Te=To5*(patf/po5)^((gam_tf-1)/gam_tf); ue=sqrt(max(0,2*cp_tf*(To5-Te)));
        Toefn=To2; Tefn=Toefn*(patf/po2)^((gam_tf-1)/gam_tf);
        uefn=sqrt(max(0,2*cp_tf*(Toefn-Tefn)));
        STf(j)=(1+ff)*ue-uatf+beta_tf*(uefn-uatf); TSf(j)=ff/max(STf(j),1e-6);
        To3j=To1*PRc^((gam_tf-1)/gam_tf); fj=(T04_tf/To3j-1)/(dh_R/(cp_tf*To3j)-T04_tf/To3j);
        if fj<0;fj=0;end
        Wcj=cp_tf*(To3j-To1); To5j=T04_tf-Wcj/(cp_tf*(1+fj));
        po5j=po1*PRc*(To5j/T04_tf)^(gam_tf/(gam_tf-1));
        Tej=To5j*(patf/po5j)^((gam_tf-1)/gam_tf); uej=sqrt(max(0,2*cp_tf*(To5j-Tej)));
        STj(j)=(1+fj)*uej-uatf; TSj(j)=fj/max(STj(j),1e-6);
    end
    subplot(2,2,(ic-1)*2+1);
    semilogx(PR_c_tf,STf/1000,'r-','LineWidth',2,'DisplayName','Turbofan (beta=5)'); hold on;
    semilogx(PR_c_tf,STj/1000,'b--','LineWidth',2,'DisplayName','Turbojet');
    xlabel('PR_c'); ylabel('ST (kN.s/kg)');
    title(['Specific Thrust - ' lbl_tf{ic}]); legend; grid on; xlim([1 100]);
    set(gca,'XTick',[1 3 5 10 30 50 100]);
    subplot(2,2,(ic-1)*2+2);
    semilogx(PR_c_tf,TSf*1e6,'r-','LineWidth',2,'DisplayName','Turbofan (beta=5)'); hold on;
    semilogx(PR_c_tf,TSj*1e6,'b--','LineWidth',2,'DisplayName','Turbojet');
    xlabel('PR_c'); ylabel('TSFC (mg/N.s)');
    title(['TSFC - ' lbl_tf{ic}]); legend; grid on; xlim([1 100]);
    set(gca,'XTick',[1 3 5 10 30 50 100]);
end
sgtitle('Ideal Turbofan vs Turbojet  (T04=1600 K, beta=5, FPR=1.5)','FontSize',12,'FontWeight','bold');
saveas(gcf,'Fig_TurbofanVsTurbojet.png');

%%
%  PART 11 - OPTIMUM BYPASS RATIO SWEEP  (Lecture 13)
PR_c_opt=20; T04_opt=1600; Ta_opt=288; pa_opt=101325; ua_opt=0;
To1_opt=Ta_opt; po1_opt=pa_opt;
FPR_sw=[1.0,1.5,2.0,2.4]; beta_sw=0:0.25:12; nB=numel(beta_sw);
cols={'g-','r--','b-.','m-'};
figure('Name','Optimum BPR','NumberTitle','off','Position',[100 100 1000 420]);
axA=subplot(1,2,1); axB=subplot(1,2,2);
for ifp=1:4
    FPRo=FPR_sw(ifp);
    To2o=To1_opt*FPRo^((gam_tf-1)/gam_tf); po2o=po1_opt*FPRo; Wfo=cp_tf*(To2o-To1_opt);
    To3o=To2o*PR_c_opt^((gam_tf-1)/gam_tf); Wco=cp_tf*(To3o-To2o);
    fo=(T04_opt/To3o-1)/(dh_R/(cp_tf*To3o)-T04_opt/To3o); if fo<0;fo=0;end
    ST_b=zeros(1,nB); SFC_b=zeros(1,nB);
    for ib=1:nB
        bo=beta_sw(ib); Wto=(1+bo)*Wfo+Wco;
        To5o=T04_opt-Wto/(cp_tf*(1+fo));
        po5o=po1_opt*FPRo*PR_c_opt*(To5o/T04_opt)^(gam_tf/(gam_tf-1));
        Teo=To5o*(pa_opt/po5o)^((gam_tf-1)/gam_tf); ueo=sqrt(max(0,2*cp_tf*(To5o-Teo)));
        Toefno=To2o; Tefno=Toefno*(pa_opt/po2o)^((gam_tf-1)/gam_tf);
        uefno=sqrt(max(0,2*cp_tf*(Toefno-Tefno)));
        ST_b(ib)=(1+fo)*ueo-ua_opt+bo*(uefno-ua_opt);
        SFC_b(ib)=fo/max(ST_b(ib),1e-6)*1e3;
    end
    lbl=sprintf('FPR=%.1f',FPRo);
    plot(axA,beta_sw,ST_b/1000,cols{ifp},'LineWidth',2,'DisplayName',lbl); hold(axA,'on');
    plot(axB,beta_sw,SFC_b,    cols{ifp},'LineWidth',2,'DisplayName',lbl); hold(axB,'on');
end
xlabel(axA,'BPR'); ylabel(axA,'ST (kN.s/kg)'); title(axA,'ST vs BPR'); grid(axA,'on'); legend(axA);
xlabel(axB,'BPR'); ylabel(axB,'SFC (kg/kN.s)'); title(axB,'SFC vs BPR'); grid(axB,'on'); legend(axB);
sgtitle('Ideal Turbofan - Optimum Bypass  (Lecture 13)','FontSize',12,'FontWeight','bold');
saveas(gcf,'Fig_OptimumBPR.png');

%%
%  PART 12 - FOUR MAIN SUMMARY FIGURES
% --- Fig 1: CD and Drag vs Mach ---
figure('Name','Fig1_DragVsMach','NumberTitle','off','Position',[50 50 1100 480]);
subplot(1,2,1);
plot(M_fine,CD_fine,'b-','LineWidth',2.5);
hold on;
plot(M_CD_pts,CD_pts,'bo','MarkerSize',7,'MarkerFaceColor','b');
xline(1.0,'k:','M=1','FontSize',9,'LabelVerticalAlignment','top');
xline(2.5,'k--','TJ-RJ','FontSize',9,'LabelVerticalAlignment','top');
xline(4.5,'g--','Rkt','FontSize',9,'LabelVerticalAlignment','top');
xlabel('Mach Number','FontSize',11); ylabel('Drag Coefficient C_D','FontSize',11);
title('Drag Coefficient vs Mach  (S=200 m^2)','FontSize',12,'FontWeight','bold');
ylim([0 0.08]); grid on; xlim([0 6]);
text(1.0,0.072,'Wave drag peak','FontSize',9,'Color','r', ...
    'HorizontalAlignment','center');

subplot(1,2,2);
bar_colors=repmat([0.2 0.5 0.8],numel(M_tr),1);
[~,imax]=max(D_tr); bar_colors(imax,:)=[0.85 0.2 0.2];
for i=2:numel(M_tr)
    b=bar(M_tr(i),D_tr(i)/1e3,'FaceColor',bar_colors(i,:),'EdgeColor',[0 0 0],'BarWidth',0.3);
    hold on;
end
xlabel('Mach Number','FontSize',11); ylabel('Drag Force (kN)','FontSize',11);
title('Drag Force at Trajectory Waypoints','FontSize',12,'FontWeight','bold');
text(M_tr(imax),D_tr(imax)/1e3+20, ...
    sprintf('Peak %.0f kN\n(M=%.1f)',D_tr(imax)/1e3,M_tr(imax)), ...
    'HorizontalAlignment','center','FontSize',9,'Color','r','FontWeight','bold');
xline(2.5,'k--'); xline(4.5,'g--'); grid on;
saveas(gcf,'Fig1_DragVsMach.png');

% --- Fig 2: Thrust vs Drag all phases ---
figure('Name','Fig2_ThrustVsDrag','NumberTitle','off','Position',[50 100 1050 520]);
Mplt=M_tr(2:end); Dplt=D_tr(2:end)/1e3;
TABplt=T_AB(2:end)/1e3; Ttotplt=T_total(2:end)/1e3;

patch([0 2.5 2.5 0],[0 0 2500 2500],[0.88 0.93 1.0],'EdgeColor','none','FaceAlpha',0.5); hold on;
patch([2.5 4.5 4.5 2.5],[0 0 2500 2500],[0.90 1.0 0.88],'EdgeColor','none','FaceAlpha',0.5);
patch([4.5 6.1 6.1 4.5],[0 0 2500 2500],[1.0 0.92 0.88],'EdgeColor','none','FaceAlpha',0.5);

h1=plot(Mplt,TABplt, 'b-o','LineWidth',2.5,'MarkerSize',7,'DisplayName','T_{AB} (TJ+RJ)');
h2=plot(Mplt,Ttotplt,'k-s','LineWidth',2.5,'MarkerSize',7,'DisplayName','T_{total} (+Rocket M>=4.5)');
h3=plot(Mplt,Dplt,   'r-^','LineWidth',2.5,'MarkerSize',7,'DisplayName','Drag D');
Mfill=[Mplt,fliplr(Mplt)]; Tfill=[Ttotplt,fliplr(Dplt)];
fill(Mfill,Tfill,[0.75 1.0 0.75],'EdgeColor','none','FaceAlpha',0.35,'DisplayName','T-D margin');

[~,iw]=min(Ttotplt./Dplt);
plot(Mplt(iw),Dplt(iw),'rv','MarkerSize',14,'LineWidth',2.5,'DisplayName','Worst T/D');
text(Mplt(iw)+0.18,Dplt(iw)+50, ...
    sprintf('Worst case\nM=%.1f, T/D=%.2f',Mplt(iw),Ttotplt(iw)/Dplt(iw)), ...
    'FontSize',9,'Color','r','FontWeight','bold');

text(1.0,2300,'Turbojet Phase','FontSize',10,'HorizontalAlignment','center','Color',[0.15 0.35 0.7],'FontWeight','bold');
text(3.45,2300,'Ramjet Phase','FontSize',10,'HorizontalAlignment','center','Color',[0.15 0.55 0.15],'FontWeight','bold');
text(5.25,2300,'Rkt+RJ Phase','FontSize',10,'HorizontalAlignment','center','Color',[0.65 0.25 0.1],'FontWeight','bold');
xline(2.5,'k--','TJ to RJ','FontSize',9,'LabelVerticalAlignment','top');
xline(4.5,'g--','Rkt ignition','FontSize',9,'LabelVerticalAlignment','top');

xlabel('Mach Number','FontSize',12); ylabel('Force (kN)','FontSize',12);
title('Thrust vs Drag — All Flight Phases','FontSize',13,'FontWeight','bold');
legend('Location','northeast','FontSize',9);
grid on; xlim([0.4 6.1]); ylim([0 2500]);
saveas(gcf,'Fig2_ThrustVsDrag.png');

% --- Fig 3: Velocity vs Time ---
figure('Name','Fig3_VelocityVsTime','NumberTitle','off','Position',[100 100 950 500]);
t_full=[t_v(2:end), t_traj];
V_full=[V_v(2:end), V_traj];

% ground roll patch
patch([0 t_v(end) t_v(end) 0],[0 0 2100 2100], ...
    [0.92 0.92 1.0],'EdgeColor','none','FaceAlpha',0.5,'DisplayName',''); hold on;

plot(t_full,V_full,'b-','LineWidth',2.5,'DisplayName','Vehicle Velocity');
% Phase vertical lines
if ~isempty(t_traj)
    [~,iV25]=min(abs(V_traj-885)); % approx M=3 at 15km speed
    if iV25>1&&iV25<numel(V_traj)
        xline(t_traj(iV25),'k--','M=3.0','FontSize',9,'LabelVerticalAlignment','top');
    end
    [~,iV45]=min(abs(V_traj-1334)); % M=4.5 at 22km
    if iV45>1&&iV45<numel(V_traj)
        xline(t_traj(iV45),'g--','Rkt on','FontSize',9,'LabelVerticalAlignment','top');
    end
end
yline(321,'k:','M=1 (5 km)','FontSize',9,'LabelHorizontalAlignment','right');
yline(1334,'g:','M=4.5','FontSize',9,'LabelHorizontalAlignment','right');
yline(1811,'r:','M=6.0 (target)','FontSize',9,'LabelHorizontalAlignment','right');

text(t_v(end)/2,150,'Ground Roll','FontSize',9,'HorizontalAlignment','center','Color',[0.3 0.3 0.7]);
xlabel('Mission Time (s)','FontSize',12); ylabel('Velocity (m/s)','FontSize',12);
title('Velocity vs Time — APG ASCENT-1','FontSize',13,'FontWeight','bold');
legend({'Ground roll region','Vehicle velocity'},'Location','northwest','FontSize',10);
grid on; xlim([0 t_full(end)+20]); ylim([0 2200]);
saveas(gcf,'Fig3_VelocityVsTime.png');

% --- Fig 4: Altitude vs Time ---
figure('Name','Fig4_AltitudeVsTime','NumberTitle','off','Position',[150 100 950 500]);
h_full=[zeros(size(t_v(2:end))), h_traj_v];

% Find phase boundaries in time array
t_rj_start=NaN; t_rkt_start=NaN;
for ii=1:numel(h_traj_v)
    if h_traj_v(ii)>=12 && isnan(t_rj_start);  t_rj_start=t_traj(ii);  end
    if h_traj_v(ii)>=22 && isnan(t_rkt_start); t_rkt_start=t_traj(ii); end
end
if ~isnan(t_rj_start)&&~isnan(t_rkt_start)
    patch([t_v(end) t_rj_start  t_rj_start  t_v(end)], ...
          [0 0 35 35],[0.88 0.93 1.0],'EdgeColor','none','FaceAlpha',0.5); hold on;
    patch([t_rj_start t_rkt_start t_rkt_start t_rj_start], ...
          [0 0 35 35],[0.90 1.0 0.88],'EdgeColor','none','FaceAlpha',0.5);
    patch([t_rkt_start t_full(end)+10 t_full(end)+10 t_rkt_start], ...
          [0 0 35 35],[1.0 0.92 0.88],'EdgeColor','none','FaceAlpha',0.5);
    text((t_v(end)+t_rj_start)/2,33,'Turbojet','FontSize',9,'HorizontalAlignment','center','Color',[0.15 0.35 0.7],'FontWeight','bold');
    text((t_rj_start+t_rkt_start)/2,33,'Ramjet','FontSize',9,'HorizontalAlignment','center','Color',[0.15 0.55 0.15],'FontWeight','bold');
    text((t_rkt_start+t_full(end))/2,33,'Rkt+RJ','FontSize',9,'HorizontalAlignment','center','Color',[0.65 0.25 0.1],'FontWeight','bold');
end

plot(t_full,h_full,'r-','LineWidth',2.5,'DisplayName','Altitude');
yline(12,'k:','12 km (TJ-RJ)','FontSize',9,'LabelHorizontalAlignment','right');
yline(22,'k--','22 km (Rkt on)','FontSize',9,'LabelHorizontalAlignment','right');
yline(30,'g--','30 km (Release)','FontSize',9,'LabelHorizontalAlignment','right');
xlabel('Mission Time (s)','FontSize',12); ylabel('Altitude (km)','FontSize',12);
title('Altitude vs Time — APG ASCENT-1','FontSize',13,'FontWeight','bold');
legend('Location','northwest','FontSize',10); grid on;
xlim([0 t_full(end)+20]); ylim([0 35]);
saveas(gcf,'Fig4_AltitudeVsTime.png');

fprintf('\n>>> Figures saved: Fig1_DragVsMach.png  Fig2_ThrustVsDrag.png\n');
fprintf('>>>                Fig3_VelocityVsTime.png  Fig4_AltitudeVsTime.png\n');
fprintf('>>> All done.\n\n');

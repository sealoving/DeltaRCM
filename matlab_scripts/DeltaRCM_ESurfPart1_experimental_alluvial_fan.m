
clc
clear all
close all

load cmapeta
load cmaph


mode_flow = 1;

VER_MATLAB = 0; % 0-old, 1-new
itermax = 1;
omega_flow_iter = 2*1/itermax;

plotnumber = 1400;
plotinterval = 5;
totaltimestep = 2000;

L = 70;
W = 130;
CTR = floor(W/2);
dx = 0.02; % (m) not sure the effect
N0 = 5; % number of cells accross inlet channel
u0 = 0.3; % (m/s) - velocity scale (set by inlet velocity)
h0 = 0.02; % (m) - depth/length scale (set by inlet depth)
hB = 1*h0; % (m)

V0 = h0*(dx*dx);
dVs = 0.1*N0^2*V0;

Qw0 = u0*h0*N0*dx;
C0 = 0.1*1/100;
Qs0 = Qw0*C0;
dt = dVs/Qs0; % N_crossdiff should be a function of dt
qw0 = u0*h0;
qs0 = qw0*C0;
Np_water = 1000;
Qp_water = Qw0/Np_water;
Np_sed = 1000;
Vp_sed = dVs/Np_sed;

GRAVITY = 9.81;
S0 = 0.02;
u_max = 2.0*u0;
H_SL = 0;
SLR = 0/1000*h0/dt; % sea-level rise per time step
dry_depth = min(0.1,0.1*h0); % (m) critial depth to switch to "dry" node

U_ref = u0;
S_ref = S0;
gamma = GRAVITY*S_ref*dx/U_ref/U_ref;

alpha = 0.1; %0.05*(0.25*1*0.125); % topo-diffusion coefficient
beta = 3; % non-linear exponent of sediment flux to flow velocity
f_bedload = 1.0;
lambda = 1; %f_bedload; % "sedimentation lag" - 1.0 means no lag

theta_wq = 1.0;
theta_sq = theta_wq;
theta_wh = 1.0;
theta_sh = 2.0*theta_wh;

U_dep = 0.75*u0;
U_ero = 1.05*u0;

Nsmooth = 10; % iteration of surface smoothing per time step
Csmooth = 0.9; % center-weighted surface smoothing

omega_sfc = 0.1;
omega_flow = 0.9;

N_crossdiff = round(dVs/V0);

% --- storage preparation
eta = zeros(L,W); % bed elevation
H = zeros(L,W); % free surface elevation
h = zeros(L,W); % depth of water
qx = zeros(L,W);
qy = zeros(L,W);
qw = zeros(L,W);
ux = zeros(L,W);
uy = zeros(L,W);
uw = zeros(L,W);

DELTA_INITIAL


% % =================== time steps =============================
wet_flag = zeros(L,W);
px_start = 1; % x-axis of inlet cells
py_start = [CTR-round(N0/2)+1:CTR-round(N0/2)+N0]; % y-axis of inlet cells
dxn_iwalk_inlet = dxn_iwalk(1); % x-component of inlet flow direction
dxn_jwalk_inlet = dxn_jwalk(1); % y-component of inlet flow direction
itmax = 2*(L+W);

Hnew = zeros(L,W);
qxn = zeros(L,W);
qyn = zeros(L,W);
qwn = zeros(L,W);
sfc_visit = zeros(L,W);
sfc_sum = zeros(L,W);
prepath_flag = zeros(L,W);

iseq = zeros(itmax,1);
jseq = zeros(itmax,1);

weight = zeros(8,1);
weight_int = zeros(8,1);
weight_sfc = zeros(8,1);
dxn = zeros(8,1);
weight_val = zeros(8,1);

qs_x = zeros(L,W);
qs_y = zeros(L,W);
qs = zeros(L,W);

timestep = 0;
fignum = 0;

% prepare for recording strata
z0 = H_SL-hB*1.1;
dz = 0.01*h0;
zmax = round((H_SL+SLR*totaltimestep*dt+S0*L/2*dx-z0)/dz);
strata0 = -totaltimestep;
strata = ones(L,W,zmax)*strata0;
topz = zeros(L,W);
for i = 1:L
    for j = 1:W
        zn = round((eta(i,j)-z0)/dz);
        zn = max(1,zn);
        topz(i,j) = zn;
        topz(i,j) = min(zmax,topz(i,j));
    end
end
strata_age = zeros(L,W);

% --- make plots ---
Emax = H_SL+h0;
Emin = -hB*1.2;
vE = [Emin:(Emax-Emin)/25:Emax];

tic
while timestep <= totaltimestep

    timestep = timestep+1;
    timestep

    for iter = 1:itermax
        % =========== water parcels ===========
        qxn = 0*qxn;
        qyn = 0*qyn;
        qwn = 0*qwn;

        wet_flag = 0*wet_flag;
        for i = 1:L
            for j = 1:W
                if h(i,j) >= dry_depth
                    wet_flag(i,j) = 1;
                end
            end
        end

        Hnew = 0*Hnew;
        sfc_visit = 0*sfc_visit;
        sfc_sum = 0*sfc_sum;

        for np = 1:Np_water
            prepath_flag = 0*prepath_flag; % "pre"-path refers to one parcel
            iseq = 0*iseq;
            jseq = 0*jseq;
            water_continue = TRUE;
            sfccredit = TRUE;

            if mod(np,100) == 0
                np
            end

            px = px_start;
            if VER_MATLAB == 0
                py = py_start(randint(1,1,[1,length(py_start)]));
            else
                py = py_start(randi([1,length(py_start)],1,1));
            end
            qxn(px,py) = qxn(px,py)+dxn_iwalk_inlet*Qp_water/dx/2;
            qyn(px,py) = qyn(px,py)+dxn_jwalk_inlet*Qp_water/dx/2;
            qwn(px,py) = qwn(px,py)+Qp_water/dx/2;

            it = 1;
            iseq(it) = px;
            jseq(it) = py;
            while water_continue == TRUE && it < itmax

                prepath_flag(px,py) = 1;
                it = it+1;

                % ========== calculate routing weights =========
                nk = Nnbr(px,py);
                weight = 0*weight;
                weight_int = 0*weight_int;
                weight_sfc = 0*weight_sfc;
                weight_val = 0*weight_val;
                for k = 1:nk
                    dxn(k) = nbr(px,py,k);
                end



                % calculate weight_int and weight_sfc
                for k = 1:nk
                    pxn = px+dxn_iwalk(dxn(k));
                    pyn = py+dxn_jwalk(dxn(k));
                    dist = dxn_dist(dxn(k));

                    if wet_flag(pxn,pyn) == 1 && wall_flag(pxn,pyn) == 0
                        weight_sfc(k) = max(0,H(px,py)-H(pxn,pyn))/dist;
                        weight_int(k) = max(0,qx(px,py)*dxn_ivec(dxn(k))+qy(px,py)*dxn_jvec(dxn(k)))/dist;
                    end
                end

                % normalize and calculate weight
                if sum(weight_sfc) ~= 0
                    % require that weightsfc >= 0
                    weight_sfc = weight_sfc/sum(weight_sfc);
                end
                if sum(weight_int) ~= 0
                    % require that weightint >= 0
                    weight_int = weight_int/sum(weight_int);
                end
                weight = gamma*weight_sfc + weight_int;
                for k = 1:nk
                    pxn = px+dxn_iwalk(dxn(k));
                    pyn = py+dxn_jwalk(dxn(k));
                    dist = dxn_dist(dxn(k));

                    if wet_flag(pxn,pyn) == 1 %&& wall_flag(pxn,pyn) == 0
                        weight(k) = h(pxn,pyn)^theta_wh*weight(k);
                    end
                end                
                % if weight is not all zeros
                if sum(weight) > 0
                    weight = weight/sum(weight);
                    % choose target cell by probability
                    for k = 1:nk
                        weight_val(k) = sum(weight(1:k));
                    end
                    step_rand = rand();
                    for k = 1:nk
                        if step_rand < weight_val(k)
                            istep = dxn_iwalk(dxn(k));
                            jstep = dxn_jwalk(dxn(k));
                            break
                        end
                    end
                end
                % if weight is all zero, do a random walk
                % if weight is all zero, do a random walk
                if sum(weight) == 0
                    if VER_MATLAB == 0
                        pxn = px+randint(1,1,[-1,1]);
                        pyn = py+randint(1,1,[-1,1]);
                        pxn = max(1,pxn);
                    else
                        pxn = px+randi([-1,1],1,1);
                        pyn = py+randi([-1,1],1,1);
                        pxn = max(1,pxn);
                    end
                    ntry = 0;
                    while wet_flag(pxn,pyn) == 0 && ntry < 5
                        ntry = ntry+1;
                        if VER_MATLAB == 0
                            pxn = px+randint(1,1,[-1,1]);
                            pyn = py+randint(1,1,[-1,1]);
                            pxn = max(1,pxn);
                        else
                            pxn = px+randi([-1,1],1,1);
                            pyn = py+randi([-1,1],1,1);
                            pxn = max(1,pxn);
                        end
                    end
                    istep = pxn-px;
                    jstep = pyn-py;
                end
            % got istep, jstep
                pxn = px+istep;
                pyn = py+jstep;
                dist = sqrt(istep^2+jstep^2);
                if dist > 0
                    qxn(px,py) = qxn(px,py)+istep*Qp_water/dx/2/dist;
                    qyn(px,py) = qyn(px,py)+jstep*Qp_water/dx/2/dist;
                    qwn(px,py) = qwn(px,py)+Qp_water/dx/2;
                    qxn(pxn,pyn) = qxn(pxn,pyn)+istep*Qp_water/dx/2/dist;
                    qyn(pxn,pyn) = qyn(pxn,pyn)+jstep*Qp_water/dx/2/dist;
                    qwn(pxn,pyn) = qwn(pxn,pyn)+Qp_water/dx/2;
                end
                px = pxn;
                py = pyn;
                iseq(it) = px;
                jseq(it) = py;
                % deal with loops
                % deal with loops
                if prepath_flag(px,py) == TRUE && it > L0
                    sfccredit = FALSE;
                    Fx = px-1;
                    Fy = py-CTR;
                    Fw = sqrt(Fx.^2+Fy.^2);
                    px = px+round(Fx/Fw*5);
                    py = py+round(Fy/Fw*5);
                    px = max(px,L0+1); px = min(L-1,px);
                    py = max(2,py); py = min(W-1,py);
                end
                if boundflag(px,py) == TRUE
                    water_continue = FALSE;
                    itend = it;
                end
            end
            if dist > 0
                qxn(px,py) = qxn(px,py)+istep*Qp_water/dx/2/dist;
                qyn(px,py) = qyn(px,py)+jstep*Qp_water/dx/2/dist;
                qwn(pxn,pyn) = qwn(pxn,pyn)+Qp_water/dx/2;
            end

            % =========== calculate free surface =============
            % calcuate free surface along one water parcel path
            % not update yet
            itback = itend; %size(iseq,2);
            if boundflag(iseq(itback),jseq(itback)) == TRUE && sfccredit == TRUE
                Hnew(iseq(itback),jseq(itback)) = H_SL;
                it0 = 0;
                Ldist = 0;
                for it = itback-1:-1:1
                    i = iseq(it); 
                    ip = iseq(it+1);
                    j = jseq(it);
                    jp = jseq(it+1);
                    dist = ((ip-i)^2+(jp-j)^2)^0.5;

                    if it0 == 0
                        if uw(i,j) > u0*0.5 || h(i,j) < 0.1*h0 %see if it is shoreline
                            it0 = it;
                        end
                        dH = 0;
                    else
                        if uw(i,j) == 0
                            dH = 0;
                        else
                            if uw(i,j)>0
                                if h(ip,jp) > 0
                                    Fr2_loc = uw(ip,jp)^2/GRAVITY/h(ip,jp);
                                else
                                    Fr2_loc = 1;
                                end
%                                 if Fr2_loc < 0.9^2
%                                     dH = max(0,Cf/GRAVITY/h(i,j)*uw(i,j)*(ux(i,j)*(ip-i)*dx+uy(i,j)*(jp-j)*dx));
%                                     dH = dH+(eta(ip,jp)-eta(i,j));
%                                     dH = dH/(1-Fr2_loc);
%                                     dH = dH+eta(i,j)-eta(ip,jp);
%                                     dH = min(dH,S0*dx*dist);
%                                     dH = max(0,dH);
%                                 else
                                    dH = S0*(ux(i,j)*(ip-i)*dx+uy(i,j)*(jp-j)*dx)/uw(i,j);
%                                 end
                            end
                        end
                    end
                    Hnew(i,j) = Hnew(ip,jp)+dH;
                    sfc_visit(i,j) = sfc_visit(i,j)+1;
                    sfc_sum(i,j) = sfc_sum(i,j)+Hnew(i,j);
                end
            end
        end % --- end of one individual water parcel, go to next parcel

        % ======= all water parcels are done, update surface
        % update free surface
        Hnew = eta+h;
        Hnew = max(Hnew,H_SL);
        for i = 1:L
            for j = 1:W
                if sfc_visit(i,j) > 0
                    Hnew(i,j) = sfc_sum(i,j)/sfc_visit(i,j);
                end
            end
        end
        Htemp = Hnew; % smoother is applied to newly calculated free surface Hnew
        for itsmooth = 1:Nsmooth
            Hsmth = Htemp;
            for i = 1:L
                for j = 1:W
                    if boundflag(i,j) ~= 1 %&& wet_flag(i,j) == 1
                        sumH = 0;
                        nbcount = 0;
                        for k = 1:Nnbr(i,j)
                            dxn(k) = nbr(i,j,k);
                            inbr = i+dxn_iwalk(dxn(k));
                            jnbr = j+dxn_jwalk(dxn(k));
                            if wall_flag(inbr,jnbr) == 0
                                sumH = sumH+Hsmth(inbr,jnbr);
                                nbcount = nbcount+1;
                            end
                        end
                        if nbcount == 0
        %                     sprintf('nbcount is zero @ (%d, %d)',i,j)
                        else
                            Htemp(i,j) = Csmooth*Hsmth(i,j)+(1-Csmooth)*sumH/nbcount;
                        end
                    end
                end
            end
        end
        Hsmth = Htemp;
        if timestep > 1
            H = (1-omega_sfc)*H+omega_sfc*Hsmth; 
        end

        %  flooding/dry-wet correction
        for i = 1:L
            for j = 1:W
                if wet_flag(i,j) == 0 % locate dry nodes
                    for k = 1:Nnbr(i,j)
                        dxn(k) = nbr(i,j,k);
                        inbr = i+dxn_iwalk(dxn(k));
                        jnbr = j+dxn_jwalk(dxn(k));
                        if wet_flag(inbr,jnbr) == 1 && H(inbr,jnbr)>eta(i,j)
                            H(i,j) = H(inbr,jnbr);
                        end
                    end
                end
            end
        end

        h = max(0,H-eta);

        % ======= update flow field and velocity field ======
        % update flow field
        if timestep > 1
            if iter == 1
                qx = qxn*omega_flow+qx*(1-omega_flow);
                qy = qyn*omega_flow+qy*(1-omega_flow);
            else
                qx = qxn*omega_flow_iter+qx*(1-omega_flow_iter);
                qy = qyn*omega_flow_iter+qy*(1-omega_flow_iter);
            end    
        else
            qx = qxn;
            qy = qyn;
        end
        qw = (qx.^2+qy.^2).^0.5; 
        if mode_flow == 1
            qx = qx.*(qwn./max(1e-10,qw));
            qy = qy.*(qwn./max(1e-10,qw));
            qw = qwn;
        end
        % apply upstream constant flux boundary condition
        qx(px_start,py_start) = qw0;
        qy(px_start,py_start) = 0;
        qw(px_start,py_start) = qw0;
        % update velocity field
        for i = 1:L
            for j = 1:W
                if h(i,j) > dry_depth && qw(i,j) > 0
                    uw(i,j) = min(u_max,qw(i,j)/h(i,j));
                    ux(i,j) = uw(i,j)*qx(i,j)/qw(i,j);
                    uy(i,j) = uw(i,j)*qy(i,j)/qw(i,j);
                else
                    ux(i,j) = 0;
                    uy(i,j) = 0;
                    uw(i,j) = 0;
                end
            end
        end

    end

    % ============== sediment transport ==============
    qs = 0*qs;

    for np_sed = 1:Np_sed
        Vp_res = Vp_sed;

        itmax = 2*(L+W);
        if mod(np_sed,100) == 0
            np_sed
        end

        px = px_start;
        if VER_MATLAB == 0
            py = py_start(randint(1,1,[1,length(py_start)]));
        else
            py = py_start(randi([1,length(py_start)],1,1));
        end
        qs(px,py) = qs(px,py)+Vp_res/2/dt/dx;

        it = 1;
        iseq(it) = px;
        jseq(it) = py;
        sed_continue = TRUE;
        while sed_continue == TRUE && it < itmax
            clear weight
            it = it+1;

            % ======== decide the next step
            % get local out dxns 1:k
            nk = Nnbr(px,py);
            weight = zeros(nk,1);
            for k = 1:nk
                dxn(k) = nbr(px,py,k);
            end

            for k = 1:nk
                pxn = px+dxn_iwalk(dxn(k));
                pyn = py+dxn_jwalk(dxn(k));
                dist = dxn_dist(dxn(k));
                % ------- both theta_sq and theta_sh participate in the
                % calculation
                weight(k) = (max(0,qx(px,py)*dxn_ivec(dxn(k))+qy(px,py)*dxn_jvec(dxn(k))))^theta_sq*...
                    h(pxn,pyn)^theta_sh/dist;
                if wet_flag(pxn,pyn) ~= 1
                    weight(k) = 0; % doesn't allow dry nodes
                end
                if wall_flag(pxn,pyn) ~= 0
                    weight(k) = 0; % doesn't allow wall nodes
                end
            end
            if sum(weight) == 0 
                for k = 1:nk
                    pxn = px+dxn_iwalk(dxn(k));
                    pyn = py+dxn_jwalk(dxn(k));
                    dist = dxn_dist(dxn(k));

                    weight(k) = 1/dist;
                    if wall_flag(pxn,pyn) == TRUE
                        weight(k) = 0;
                    end
                end
            end

            weight = weight/sum(weight);

            for k = 1:nk
                weight_val(k) = sum(weight(1:k));
            end
            step_rand = 1-rand();
            for k = 1:nk
                if step_rand < weight_val(k)
                    istep = dxn_iwalk(dxn(k));
                    jstep = dxn_jwalk(dxn(k));
                    break
                end
            end
            dist = sqrt(istep^2+jstep^2);

            if dist > 0
                qs(px,py) = qs(px,py)+Vp_res/2/dt/dx;
            end % exit accumulation

            px = px+istep;
            py = py+jstep;

            if dist > 0
                qs(px,py) = qs(px,py)+Vp_res/2/dt/dx;
            end % entry accumulation            
            % =========== deposition and erosion at one step
            U_loc = uw(px,py);
            qs_cap = qs0/u0^beta*U_loc^beta;
            qs_loc = qs(px,py);
            eta_change_loc = 0;
            if qs_loc > qs_cap
                Vp_dep = Vp_res*lambda;
                Vp_dep = min(Vp_dep,(H(px,py)-eta(px,py))/4*(dx*dx));

                eta_change_loc = Vp_dep/(dx*dx);
                eta(px,py) = eta(px,py)+eta_change_loc;
                h(px,py) = H(px,py) - eta(px,py);
                uw(px,py) = min(u_max,qw(px,py)/h(px,py));
                if qw(px,py) > 0
                    ux(px,py) = uw(px,py)*qx(px,py)/qw(px,py);
                    uy(px,py) = uw(px,py)*qy(px,py)/qw(px,py);
                else
                    ux(px,py) = 0;
                    uy(px,py) = 0;
                end
                Vp_res = Vp_res-Vp_dep;
            elseif U_loc > U_ero && qs_loc < qs_cap
                Vp_ero = Vp_sed*(U_loc^beta-U_ero^beta)/U_ero^beta;
                Vp_ero = min(Vp_ero,(H(px,py)-eta(px,py))/4*(dx*dx));

                eta_change_loc = -Vp_ero/(dx*dx);
                eta(px,py) = eta(px,py)+eta_change_loc;
                h(px,py) = H(px,py) - eta(px,py);
                uw(px,py) = min(u_max,qw(px,py)/h(px,py));
                if qw(px,py) > 0
                    ux(px,py) = uw(px,py)*qx(px,py)/qw(px,py);
                    uy(px,py) = uw(px,py)*qy(px,py)/qw(px,py);
                else
                    ux(px,py) = 0;
                    uy(px,py) = 0;
                end
                Vp_res = Vp_res+Vp_ero;
            end

            zn = round((eta(px,py)-z0)/dz);
            zn = max(1,zn);
            new_topz_loc = max(1,zn);
            if eta_change_loc > 0 % deposition happens
                strata_age(px,py) = timestep;
                for z = topz(px,py):new_topz_loc
                    strata(px,py,z) = strata_age(px,py);
                end
            end
            if eta_change_loc < 0 % erosion happens
                for z = new_topz_loc:topz(px,py)
                    strata(px,py,z) = -totaltimestep;
                end
            end
            topz(px,py) = new_topz_loc;

            if boundflag(px,py) == TRUE
                sed_continue = FALSE;
            end
        end % --- end of one individual sediment parcel
    end % --- end of all sediment parcels

    % ---- topo diffusion after all sediment parcels
    for crossdiff = 1:N_crossdiff
        eta_diff = eta;
        for i = 2:L-1
            for j = 2:W-1
                if boundflag(i,j) == 0 && wall_flag(i,j) == 0
                    crossflux = 0;
                    for k = 1:Nnbr(i,j)
                        dxn(k) = nbr(i,j,k);
                        inbr = i+dxn_iwalk(dxn(k));
                        jnbr = j+dxn_jwalk(dxn(k));
                        if wall_flag(inbr,jnbr) == 0
                            crossflux_nb = dt/N_crossdiff*alpha*0.5*(qs(i,j)+qs(inbr,jnbr))*dx*(eta(inbr,jnbr)-eta(i,j))/dx;
                            crossflux = crossflux+crossflux_nb;

                            eta_diff(i,j) = eta_diff(i,j)+crossflux_nb/dx/dx;
                            zn = round((eta_diff(i,j)-z0)/dz);
                            zn = max(1,zn);
                            new_topz_loc = max(1,zn);
                            if crossflux_nb > 0
                                strata_age(i,j) = strata_age(inbr,jnbr);
                                for z = topz(i,j):new_topz_loc
                                    strata(i,j,z) = strata_age(i,j);
                                end
                            end
                            if crossflux_nb < 0
                                for z = new_topz_loc:topz(i,j)
                                    strata(i,j,z) = -totaltimestep;
                                end
                            end
                            topz(i,j) = new_topz_loc;

                        end
                    end

                end
            end
        end
        eta = eta_diff;
    end

    %  flooding/dry-wet correction
    for i = 1:L
        for j = 1:W
            if wet_flag(i,j) == 0 % locate dry nodes
                for k = 1:Nnbr(i,j)
                    dxn(k) = nbr(i,j,k);
                    inbr = i+dxn_iwalk(dxn(k));
                    jnbr = j+dxn_jwalk(dxn(k));
                    if wet_flag(inbr,jnbr) == 1 && H(inbr,jnbr)>eta(i,j)
                        H(i,j) = H(inbr,jnbr);
                    end
                end
            end
        end
    end
    h = H-eta; % no need to update velocities because they are not used in the following updates

    % upstream boundary condition - constant depth
    eta(px_start,py_start) = H(px_start,py_start)-h0;

    H_SL = H_SL+SLR*dt;

    R = (eta+hB)/(hB+0.5*S0*L*dx+H_SL);
    G = R;
    B = R;
    IMG(:,:,1) = R;
    IMG(:,:,2) = G;
    IMG(:,:,3) = B;
    for i = 1:L
        for j = 1:W
            if uw(i,j) >= u0*0.5 && qw(i,j) >= qw0*0.05
                kk = qw(i,j)/qw0;                
                IMG(i,j,1) = min(1,(0*kk+1*(1-kk))*IMG(i,j,1));
                IMG(i,j,2) = min(1,(0.72*kk+1*(1-kk))*IMG(i,j,2));
                IMG(i,j,3) = min(1,(1*kk+1*(1-kk))*IMG(i,j,3));
%                 IMG(i,j,3) = min(1,IMG(i,j,3)*(1+1.0*qw(i,j)/qw0));
%                 IMG(i,j,2) = min(1,IMG(i,j,2)*(1+0.5*qw(i,j)/qw0));
            end
        end
    end
    figure(1)
    imshow(IMG)
    figure(2)
    imagesc(H)
    axis equal
    colorbar
    pause(0.01)
    if mod(timestep,plotinterval) == 0
        fignum = fignum+1;
        fignum2 = fignum+plotnumber;
         handl = figure(1);
         saveas(handl,sprintf('WETIMG%d',fignum2),'png')
 %         save(sprintf('strata%d',fignum2),'strata')
    end
end
 
toc
 
handl = figure(2);
saveas(handl,sprintf('HFig%d',plotnumber),'png')
 

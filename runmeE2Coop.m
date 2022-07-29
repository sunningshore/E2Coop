% function [trajectAvoiding, trajectSet] = runmeE2Coop (velAmpIntru0, posIniIndiv)
%% test %%
clear;
clc;
plotFlag=1;
addpath('GVF');
addpath('snake2D');
%% Define intruders %%
velAmpIntru1=0:10; % Intruder speed 
n1=2:10; % Swarm size 
gamma1Set=.1:.1:.9; % Lamda_1 in paper  
varLen=length(velAmpIntru1);
for ind_l=1: 10 % Episode 
    for ind_ll=1: varLen
        gamma1 = gamma1Set(5); % Lamda_1 in paper  
        gamma2 = 1-gamma1; % Lamda_2 in paper  
        
        velAngIntru=[-135]; % Velocity direction of intruder
        posIniIntru=[180; 180]; % Position of intruder

        [~, intruNum]=size(posIniIntru);
        velAmpSwarm=10; % Velocity of swarm 
        velAmpIntru0=velAmpIntru1(ind_ll); % Velocity magnitude of intruder
        velAmpIntru = max(velAmpIntru0, velAmpSwarm);
        velAmpIntru = velAmpIntru * ones(intruNum, 1); 
        velVecIntru=[velAmpIntru.*cosd(velAngIntru), velAmpIntru.*sind(velAngIntru)];
        %% Define destination %%
        posIniDest=[250; 250]; % Swarm destination
        %% Define swarm %%
        n=n1(2); % Swarm size 
        % Arrange positions of swarm members 
        theta=linspace(-pi,pi,n+1);
        x0=95;
        y0=95;
        r=20;
        xi=r*cos(theta)+x0;
        yi=r*sin(theta)+y0;
        posIniIndivAlt=zeros(1, n);
        posIniIndiv=[xi(1:end-1); yi(1:end-1)];
        
        [~, swarmSize]=size(posIniIndiv);
        posIniSwarm=sum(posIniIndiv,2)./swarmSize;
        distSwarmIndiv = sqrt(sum((posIniIndiv-posIniSwarm).^2,1)); % Geometrical center of the swarm
        vecSwarm2Dest=posIniDest-posIniSwarm; % Offset vector towards the destination
        vecSwarm2DestUni=vecSwarm2Dest./norm(vecSwarm2Dest);
        posIniSwarm=posIniSwarm+vecSwarm2DestUni.*2*max(max(distSwarmIndiv), r); % Swarm center with offset
        %% Potential Fields %%
        fieldSize = 300; % Define size of the field (limit range of movement);
        posIniSwarmRel = round(posIniSwarm);
        posIniIntruRel = round(posIniIntru);
        groundIni = zeros(fieldSize);
        [groundRow, groundCol]=size(groundIni);
        fieldMotionIni = groundIni;
        %% Parameter Settings %%
        simulTime = 50;
        reslS=30; % Granularity of deltaS (length of arc)
        deltaSMax = max(velAmpIntru0, 5)*ones(1,1,swarmSize); % Maximum deltaS
        deltaSMin = max(velAmpIntru0, 5)*ones(1,1,swarmSize); % Manimum deltaS
        % deltaSMax = deltaSMax.*reshape((linspace(1,2,swarmSize)), 1,1,swarmSize); % Maximum deltaS
        % deltaSMin = deltaSMin.*reshape((linspace(1,2,swarmSize)), 1,1,swarmSize); % Manimum deltaS
        linSpcInd = linspace(0,1,reslS);
        deltaS = deltaSMin+linSpcInd'.*(deltaSMax-deltaSMin); % Dimension of arc length in PSO
                
        omega=0:1:359;
        omega=omega(:,:,ones(swarmSize, 1)); % Dimension of omega in PSO
        
        linSpcInd = linspace(0,1,350);
        kappa=ones(350, 1)*[-pi./(1*deltaSMax(:)')]+linSpcInd'.*(ones(350, 1)*(pi./(1*deltaSMax(:)')+pi./(1*deltaSMax(:)')));
        kappa=reshape(kappa, 1, 350, swarmSize); % Dimension of curvature in PSO
        
        vetVelMax=30*ones(1,1,swarmSize);
        linSpcInd = linspace(0,1,2*5+1);
        tau=-vetVelMax+linSpcInd'.*(vetVelMax+vetVelMax);
        % tau=0+linSpcInd'.*(vetVelMax-0);
        % tau=tau*0;
        
        Maxiter=100; % Iteratons of optimization in PSO
        N=100; % number of interpolation points of each arc
        %% Pre-trajectory Settings %%
        smoothField = groundIni;
        xPre = reshape(posIniIndiv(1,:),1,1,swarmSize);
        yPre = reshape(posIniIndiv(2,:),1,1,swarmSize);
        AltPre = reshape(posIniIndivAlt,1,1,swarmSize);
        kPre = zeros(1, 1, swarmSize);
        tPre = zeros(1, 1, swarmSize);
        omegaPre = zeros(1, 1, swarmSize);
        contourPlotCell=cell(simulTime, 1);
%         [arcxPre, arcyPre, x1Pre, y1Pre]=arcParameterization (xPre, yPre, kPre, omegaPre, 5*5, N);
        [arcxPre, arcyPre, arcAltPre, x1Pre, y1Pre]=arcParameterization3D (xPre, yPre, AltPre, kPre, tPre, omegaPre, 5*5, N);
        arcxPre=arcxPre-5*5;
        arczPre=interp2(smoothField, arcxPre, arcyPre); % Previous trajectories before avoidance
        x0=xPre;
        y0=yPre;
        Alt0=AltPre;
        z0=interp2(smoothField, x0, y0); % Starting points for trajectory planning
        
        arcPreOrigin=cat(2,arcxPre,arcyPre,arcAltPre);
        contourPlot=[0 0 0];
        trajectSet=[x0, y0, AltPre, zeros(size(x0)), z0];
        posSwarmSet=posIniSwarmRel;
        smoothFieldSet=smoothField;
        %% Plotting %%
        if plotFlag
            writerObj1=VideoWriter('AnimiFinal.avi');
            writerObj1.FrameRate = 1;
            open(writerObj1);
            
            fig2=figure(ind_ll);
            g3=mesh(gca, smoothField);
            axe2=get(fig2,'CurrentAxes');
            view (0,90);
            hold (axe2);
            g1=plot3(axe2, x0(:), y0(:), AltPre(:), 'r*');
            for ind_h=1: swarmSize
                eval(['h',num2str(ind_h),'=plot3(axe2, trajectSet(:,1,ind_h), trajectSet(:,2,ind_h), ones(size(trajectSet(:,3,ind_h)))*max(smoothField(:)), ''r'',''linewidth'',2); ']);
            end
            g4=plot3(contourPlot(:,1), contourPlot(:,2), contourPlot(:,3), 'w.','markersize',3);
            for ind_k=1: swarmSize
                eval(['k',num2str(ind_k),'=plot3(axe2, arcPreOrigin(:,1,ind_k), arcPreOrigin(:,1,ind_k), ones(size(arcPreOrigin(:,1,ind_k)))*max(smoothField(:)), ''r:'',''linewidth'',2); ']);
            end
            g2=plot3(posIniIntru(1,:), posIniIntru(2,:), interp2(smoothField, posIniIntru(1,:), posIniIntru(2,:)), 'ro');
            % g5=plot3(posIniDest(1,:), posIniDest(2,:), max(smoothField(:)), 'rs');
            for ind_p=1: swarmSize
                eval(['p',num2str(ind_p),'=plot3(axe2, arcPreOrigin(:,1,ind_p), arcPreOrigin(:,1,ind_p), ones(size(arcPreOrigin(:,1,ind_p)))*max(smoothField(:)), ''b:'',''linewidth'',2); ']);
            end
        end
        %% Optimization %%
        posSwarmRel=posIniSwarmRel; % Position of swarm
        posIntruRel=posIniIntruRel; % Position of intruder
        vecSwarm2Dest=posIniDest-posSwarmRel; % Vector from swarm center to destination
        vecSwarm2DestUni=vecSwarm2Dest./norm(vecSwarm2Dest);
        distSwarm2Dest=norm(vecSwarm2Dest);
        
        distIntruIndiv = sqrt(sum((posIniIndiv(:,:,ones(intruNum, 1))-reshape(posIniIntru, 2, 1, intruNum)).^2,1));
        distIntruIndiv = min(distIntruIndiv, [], 3); % Distance from swarm member to intruder
        
        vecIntruIndiv = reshape(posIniIntru, 2, 1, intruNum)-posIniIndiv(:,:,ones(intruNum, 1));
        vecIntruIndivUni = vecIntruIndiv./ vecnorm(vecIntruIndiv); % Vector from swarm member to intruder
        
        vecDestIndiv = posIniDest - posIniIndiv;
        vecDestIndivUni = vecDestIndiv./ vecnorm(vecDestIndiv); % Vector from swarm member to destination
        
        stateIndicator = dot(vecIntruIndivUni, vecDestIndivUni(:,:,ones(1, intruNum)));
        % stateIndicator is the dot product of vector UAV-intruder and vector
        % UAV-destination. It indicates when the UAVs finish avoiding. If the dot
        % product is positive, posibilities of collision exist. If the dot product
        % is nagative, the UAVs have passed the intruders, avoidance is deemed accomplished.
        stateIndicator (stateIndicator >= 0) = 0;
        stateIndicator (stateIndicator < 0) = -1;
        stateIndicator = sum(stateIndicator, 3);
        stateIndicator (stateIndicator == 0) = 1;
        
        vecIniSwarmIndiv =  posIniIndiv - posIniSwarmRel; % Records the initial relative
        % positions of swarm members. It's for formation restoration after
        % avoidance.
        
        posIndiv=posIniIndiv;
        posIndivAlt=posIniIndivAlt;
        swarmInd=1:swarmSize;
        guardDistSwarm=5;
        guardDist=10; % Distance around the intruder, within which no UAV is allowed.

        startDist=50+guardDist+velAmpIntru0; % Distance from the intruder, where UAVs start to perform avoidance.
        startDist=startDist(1, ones(1, swarmSize));
        deltaSTemp=deltaS(1,:,:);
        deltaSTemp=deltaSTemp(:);
        % First sort by velocity;
        [uniqS, uniqSInd]=unique(deltaSTemp);
        if length(uniqSInd) > 1
            sortMethod='ascend';
            [uniqSSort1, uniqSSortInd1]=sort(uniqS, sortMethod); % Slow ones avoid first.
            uniqSSortInd1=uniqSInd(uniqSSortInd1);
            uniqSdiff=abs(diff(uniqSSort1));
            if strcmp(sortMethod, 'descend')
                uniqSSortInd1=uniqSSortInd1(1: end-1);
            elseif strcmp(sortMethod, 'ascend')
                uniqSSortInd1=uniqSSortInd1(2: end);
            end
            uniqSSortInd1=uniqSSortInd1';
            % Then sort by dist2Obs;
            ind2sort=setdiff(swarmInd, uniqSSortInd1);
            [uniqSSort2, uniqSSortInd2]=sort(distIntruIndiv(ind2sort), 'ascend');
            uniqSSortInd2=ind2sort(uniqSSortInd2);
            if strcmp(sortMethod, 'ascend')
                uniqSSortInd=[uniqSSortInd2, uniqSSortInd1];
            elseif strcmp(sortMethod, 'descend')
                uniqSSortInd=[uniqSSortInd1, uniqSSortInd2];
            end
        else
            % Then sort by dist2Obs;
            [uniqSSort2, uniqSSortInd2]=sort(distIntruIndiv, 'ascend');
            uniqSSortInd=uniqSSortInd2;
        end
        
        uniqSLen=length(uniqSSortInd);
        guardDistv2v = 20; % V2V safety distance.
        
        buff=40;
        detectDist=startDist+buff; % Distance from the intruder, where UAVs detect the intruder.
        swarmIndStop = [];
        swarmIndAvoiding = [];
        swarmInd2Avoid = [];
        swarmCostSet = [];
        internCostSet = [];
        trajectAvoiding = [];
        velCenter = 0;
        velIntru = velAmpIntru0;
        fieldMotion = fieldMotionIni;
        vecIndivFormNorm = ones(1, swarmSize)*5+1;
        dist2intru = [];
        xPred=trajectSet(end,1,:);
        yPred=trajectSet(end,2,:);
        PSet = cell(swarmSize, 1);
        % for ind_t=1 %1: 3 %simulTime
        ind_t=0;
        while length(find(vecIndivFormNorm<5)) ~= swarmSize && ind_t < 40 
            ind_t=ind_t+1;
            %% Define Swarm Dir %%
            vecSwarm2Dest=posIniDest-posSwarmRel;
            distSwarm2Dest=norm(vecSwarm2Dest);
            vecSwarm2DestUni=vecSwarm2Dest./distSwarm2Dest; % Vector from swarm center to the destination
            
            if isempty(stateIndicator (stateIndicator>0))
                velCenter = 0; % If avoidance accomplished, stop moving swarm center and start formation restoration.
                %     elseif ~isempty(swarmInd2Avoid)
                %         velCenter = velCenter-velAmpIntru;
            else
                velCenter = mean(deltaSMax);
                %         velCenter = abs(velAmpIntru- velCenter);
            end
            
            posSwarmRel=posSwarmRel+vecSwarm2DestUni.*abs(velCenter); % Update swarm center movement
            %     posSwarmRel=posSwarmRel+vecSwarm2DestUni.*(velCenter);
            
            posSwarmForm = posSwarmRel;
            posIndivForm = posSwarmForm + vecIniSwarmIndiv;
            vecIndivForm = posIndivForm - posIndiv;
            vecIndivFormNorm = vecnorm(vecIndivForm);
            vecIndivFormUni = vecIndivForm./vecIndivFormNorm; % Formation restoration
            
            vecIntru=[cosd(velAngIntru); sind(velAngIntru)];
            posIntruRel = posIntruRel + vecIntru.*abs(velIntru); % Update intruder movement
            %% Construct potential field and contours %%
            obsFlag=1;
            [posCurrIntru, smoothField1] = testUAVMotion ...
                (fieldMotion, posIntruRel, velAmpIntru*10, guardDist, obsFlag); % Construct potential field
            
            %     smoothField1 = smoothField1./ max(smoothField1(:));
            
            swarmInd2Avoid = swarmInd(distIntruIndiv <= detectDist);
            if isempty(swarmInd2Avoid)
                swarmInd2Avoid = zeros(1,0);
            elseif size(swarmInd2Avoid, 1) > size(swarmInd2Avoid, 2)
                swarmInd2Avoid = swarmInd2Avoid';
            end
            
            ind_tMark = 1;
            if ~isempty(swarmInd2Avoid)
                % Accumulation field without interpolation
                obsFlag=1;
                %         guardDistSwarm=5;
                [posCurrIntru, smoothField0] = testUAVMotion (fieldMotion, ...
                    [posSwarmSet(:,ind_tMark:end)], [velAmpSwarm(1,ones(ind_t-ind_tMark+1,1))*10], guardDistSwarm, obsFlag);
                smoothField = smoothField1.*(ind_t-ind_tMark+1) + smoothField0;
                % Memory of potential field, from when the obstacle is deteced.
            else
                obsFlag=1;
                %         guardDistSwarm=5;
                [posCurrIntru, smoothField0] = testUAVMotion(fieldMotion, ...
                    [posSwarmRel], [velAmpSwarm*10], guardDistSwarm, obsFlag);
                smoothField = smoothField1 + smoothField0;
                ind_tMark = ind_t; % Mark the time obstacle is detected.
            end
            smoothFieldSet=cat(3, smoothFieldSet, smoothField);
            %% PSO %%
            if plotFlag
                tic;
            end
            
            z0=interp2(smoothField, x0, y0); % Update z0 immediately after updating smoothField.
            
            %     unique(z0(:)')
            
            intruMember=[];
            zForbid = zeros(intruNum, 1);
            for ind_i=1: intruNum % No contours drawn within guardDist.
                zForbidMat=smoothField([max(round(posIntruRel(2,ind_i)-guardDist), 1)...
                    :round(posIntruRel(2,ind_i)+guardDist)], ...
                    [max(round(posIntruRel(1,ind_i)-guardDist), 1)...
                    :round(posIntruRel(1,ind_i)+guardDist)]);
                if ~isempty(zForbidMat)
                    zForbid(ind_i)=min(zForbidMat(:));
                end
                distIntruIndiv3D=reshape(distIntruIndiv, 1, 1, swarmSize);
                % %         intruMember=cat(2, intruMember, swarmInd(distIntruIndiv3D <= guardDist+50 & z0>=zForbid(ind_i)));
                %         z0(distIntruIndiv3D <= guardDist+50 & z0 >= zForbid(ind_i)+0) = zForbid(ind_i);
            end
            
            binField=zeros(size(smoothField, 1), size(smoothField, 2), swarmSize);
            EextField=zeros(size(smoothField, 1), size(smoothField, 2), swarmSize);
            paraField=0;
            pxField=zeros(size(smoothField, 1), size(smoothField, 2), swarmSize);
            pyField=zeros(size(smoothField, 1), size(smoothField, 2), swarmSize);
            uField=zeros(size(smoothField, 1), size(smoothField, 2), swarmSize);
            vField=zeros(size(smoothField, 1), size(smoothField, 2), swarmSize);
            for ind_s=1: swarmSize
                [binField(:,:,ind_s), EextField(:,:,ind_s), pxField(:,:,ind_s), ...
                    pyField(:,:,ind_s), uField(:,:,ind_s), vField(:,:,ind_s)] = ...
                    ExternalEnergyField (smoothField, z0(ind_s));
            end
            % Generate binary field and external energy field.
            %% PSO Here %%
            arcx=zeros(N, 1, swarmSize);
            arcy=zeros(N, 1, swarmSize);
            arcz=zeros(N, 1, swarmSize);
            arcAlt=zeros(N, 1, swarmSize);
            %     arcAlt=arcAltPre;
            
            swarmIndAvoiding = swarmInd(distIntruIndiv <= startDist & stateIndicator >= 0);
            % Find UAVs avoiding obstacles
            if isempty(swarmIndAvoiding)
                swarmIndAvoiding=zeros(1,0);
            elseif size(swarmIndAvoiding, 1) > size(swarmIndAvoiding, 2)
                swarmIndAvoiding=swarmIndAvoiding';
            end
            
            swarmIndNonAvoiding = swarmInd(distIntruIndiv > startDist & stateIndicator >= 0);
            % Find UAVs not start avoidance yet
            swarmIndStop = swarmInd(stateIndicator < 0);
            
            swarmIndStop=union(swarmIndStop, swarmIndNonAvoiding);
            % Find UAVs finished avoidance
            
            if isempty(stateIndicator (stateIndicator>0))
                swarmIndStop = swarmInd;
                swarmIndAvoiding = [];
                swarmIndNonAvoiding = [];
            end
            
            if isempty(swarmIndAvoiding)
                % UAVs not start avoidance yet just go straight to destination.
                kPre2 = zeros(1, 1, length(swarmIndNonAvoiding));
                tPre2 = zeros(1, 1, length(swarmIndNonAvoiding));
                omegaPre2 = zeros(1, 1, length(swarmIndNonAvoiding));
                omegaPre2(1,1,:)=atand(vecSwarm2DestUni(2)./vecSwarm2DestUni(1));
                %         [arcx2, arcy2, x1Pre2, y1Pre2]=arcParameterization ...
                %             (x0(:,:,swarmIndNonAvoiding), y0(:,:,swarmIndNonAvoiding), kPre2, omegaPre2, deltaSMax(:,:,swarmIndNonAvoiding), N);
                [arcx2, arcy2, arcAlt2, x1Pre2, y1Pre2]=arcParameterization3D ...
                    (x0(:,:,swarmIndNonAvoiding), y0(:,:,swarmIndNonAvoiding), Alt0(:,:,swarmIndNonAvoiding), kPre2, tPre2, omegaPre2, deltaSMax(:,:,swarmIndNonAvoiding), N);
                
                arcx(:,1,swarmIndNonAvoiding)=arcx2;
                arcy(:,1,swarmIndNonAvoiding)=arcy2;
                arcAlt(:,1,swarmIndNonAvoiding)=arcAlt2;
                omegaPre (:,:,swarmIndNonAvoiding) = omegaPre2;
                
                % UAVs finihsed avoidance proceed to formation restoration.
                deltaSMax(:,:,swarmIndStop)=5;
                deltaSMin(:,:,swarmIndStop)=5;
                kPre2 = zeros(1, 1, length(swarmIndStop));
                tPre2 = zeros(1, 1, length(swarmIndStop));
                omegaPre2 = zeros(1, 1, length(swarmIndStop));
                %         omegaPre2(1,1,:)=atand(vecIndivForm(2,swarmIndStop)./vecIndivForm(1,swarmIndStop));
                angTemp=atand(vecIndivFormUni(2,swarmIndStop)./vecIndivFormUni(1,swarmIndStop));
                angTemp(vecIndivFormUni(1,swarmIndStop)<0)=180+angTemp(vecIndivFormUni(1,swarmIndStop)<0);
                omegaPre2(1,1,:)=angTemp;
                %         [arcx2, arcy2, x1Pre2, y1Pre2]=arcParameterization ...
                %             (x0(:,:,swarmIndStop), y0(:,:,swarmIndStop), kPre2, omegaPre2, deltaSMax(:,:,swarmIndStop), N);
                
                [arcx2, arcy2, arcAlt2, x1Pre2, y1Pre2]=arcParameterization3D ...
                    (x0(:,:,swarmIndStop), y0(:,:,swarmIndStop), Alt0(:,:,swarmIndStop), kPre2, tPre2, omegaPre2, deltaSMax(:,:,swarmIndStop), N);
                
                arcx(:,1,swarmIndStop)=arcx2;
                arcy(:,1,swarmIndStop)=arcy2;
                arcAlt(:,1,swarmIndStop)=arcAlt2;
                omegaPre (:,:,swarmIndStop) = omegaPre2;
                
                x0temp=x0(:,:,vecnorm(vecIndivForm)<=5);
                y0temp=y0(:,:,vecnorm(vecIndivForm)<=5);
                Alt0temp=Alt0(:,:,vecnorm(vecIndivForm)<=5);
                arcx(:,1, vecnorm(vecIndivForm)<=5)=x0temp(ones(N,1),:,:);
                arcy(:,1, vecnorm(vecIndivForm)<=5)=y0temp(ones(N,1),:,:);
                arcz(:,1, vecnorm(vecIndivForm)<=5)=ones(size(arcy(:,1, vecnorm(vecIndivForm)<=5)))*100;
                arcAlt(:,1, vecnorm(vecIndivForm)<=5)=Alt0temp(ones(N,1),:,:);
            else
                if ~isempty(intruMember)
                    arcx=arcxPre;
                    arcy=arcyPre;
                    arcAlt=arcAltPre;
                else
                    %% Prediction %%
                    iteration = 100;
                    
                    theta1 = linspace(1*pi/2, 7*pi/4, 10);
                    theta2 = linspace(-1*pi/2, 3*pi/4, 10);
                    d = 2*r;
                    xStart1 = d*cos(theta1)+sum(x0(:))/n;
                    xStart2 = d*cos(theta2)+posIntruRel(1);
                    yStart1 = d*sin(theta1)+sum(y0(:))/n;
                    yStart2 = d*sin(theta2)+posIntruRel(2);
                    xStart = [xStart1, xStart2, xStart1(1)];
                    yStart = [yStart1, yStart2, yStart1(1)];

                    trajectStart = [yStart', xStart'];
                    trajectStart = InterpolateContourPoints2D(trajectStart, N);
                    trajectStart = trajectStart(:, :, ones(1,n));
                    
                    for ind_y = 1: swarmSize
                        trajectStartTemp = PSet{ind_y};
                        if ~isempty(trajectStartTemp)
                            trajectStart(:,:,ind_y) = trajectStartTemp;
                        end
                    end
                    
                    stepPred = 20;
                    xPred = zeros(stepPred*N, 1, swarmSize);
                    yPred = zeros(stepPred*N, 1, swarmSize);
                    PSet = cell(swarmSize, 1);
                    
                    xPredStep = xPred(N+1:2*N,:,:);
                    yPredStep = yPred(N+1:2*N,:,:);
                    init_curv = sum(sqrt((xPredStep(3:end,:,:)-2*xPredStep(2:end-1,:,:)+xPredStep(1:end-2,:,:)).^2 ...
                        +(yPredStep(3:end,:,:)-2*yPredStep(2:end-1,:,:)+yPredStep(1:end-2,:,:)).^2), 1);
                    init_omega = atand((yPredStep(end,:,:) - y0)./(xPredStep(end,:,:) - x0));
                    
                    [xPred(:,:,swarmIndAvoiding), yPred(:,:,swarmIndAvoiding), P] = snakePredict (iteration, ...
                        binField(:,:,swarmIndAvoiding), trajectStart(:,:,swarmIndAvoiding), ...
                        x0(:,:,swarmIndAvoiding), y0(:,:,swarmIndAvoiding), posIniDest, ...
                        arcxPre(:,:,swarmIndAvoiding), arcyPre(:,:,swarmIndAvoiding), arczPre(:,:,swarmIndAvoiding), ...
                        trajectSet(:,:,swarmIndAvoiding));
                    
                    for ind_u = 1: length(swarmIndAvoiding)
                        PSet{swarmIndAvoiding(ind_u)} = P(:,:,ind_u);
                    end
                    %% De-conflict %%
                    AltAdjust = Alt0*0;
                    if length(swarmIndAvoiding) > 1
                        dist2Obs = min(sqrt((xPred-posIntruRel(1)).^2 + (yPred-posIntruRel(2)).^2), [], 1);
                        combo = nchoosek(swarmIndAvoiding, 2);
                        distV2V = sqrt((xPred(end,:,combo(:,1))-xPred(end,:,combo(:,2))).^2 ...
                            + (yPred(end,:,combo(:,1))-yPred(end,:,combo(:,2))).^2 ...
                            + (Alt0(end,:,combo(:,1))-Alt0(end,:,combo(:,2))).^2);
                        
                        squeeze(distV2V)
                        
                        thrV2V = 20;
                        dangerInd = combo(distV2V < thrV2V, :);
                        %% Alt PSO %%
                        %                 AltAdjust = Alt0;
                        if ~isempty(dangerInd)
                            t = tau(:,:,swarmIndAvoiding);
                            xEnd = xPred(end,:,swarmIndAvoiding);
                            yEnd = yPred(end,:,swarmIndAvoiding);
                            zEnd = z0(:,:,swarmIndAvoiding);
                            AltEnd = Alt0(:,:,swarmIndAvoiding);
                            Maxiter = 100;
                            [tBestInterp] = SpeciesPSOAlt (xEnd, yEnd, zEnd, ...
                                AltEnd, t, thrV2V, Maxiter);
                        else
                            tBestInterp = 0;
                        end
                    else
                        tBestInterp = 0;
                    end
                    AltAdjust(:,:,swarmIndAvoiding) = tBestInterp;
                    altTran=Alt0(ones(N,1),:,:);
                    for ind_j = swarmIndAvoiding
                        altStart = Alt0(ind_j);
                        altEnd = AltAdjust(ind_j)+Alt0(ind_j);
                        altTran(:,:,ind_j) = linspace(altStart, altEnd, N);
                    end
                    %% Horizontal PSO %%
                    [arcx1, arcy1, arcz1, kNew, omegaNew, deltaSNew, BestCosts, swarmVecPre, swarmVecCurr, swarmCost] = ...
                        SpeciesPSO (gamma1, gamma2, swarmIndAvoiding, posIniIndiv(:,  swarmIndAvoiding), posIniSwarmRel, ...
                        posSwarmRel, x0(:,:, swarmIndAvoiding), y0(:,:, swarmIndAvoiding), z0(:,:, swarmIndAvoiding), ...
                        zForbid, kappa(:,:, swarmIndAvoiding), omega(:,:, swarmIndAvoiding), deltaS(:,:, swarmIndAvoiding), ...
                        arcxPre(:,:, swarmIndAvoiding), arcyPre(:,:, swarmIndAvoiding), arczPre(:,:, swarmIndAvoiding), ...
                        smoothField, EextField, Maxiter, N, init_curv, init_omega);
                    %             arcAlt1=arcz1;
                    arcAlt1=Alt0(ones(N,1),:,swarmIndAvoiding);
                    
                    arcx(:,1, swarmIndAvoiding)=arcx1;
                    arcy(:,1, swarmIndAvoiding)=arcy1;
                    arcz(:,1, swarmIndAvoiding)=arcz1;
                    arcAlt(:,1, swarmIndAvoiding)=altTran(:,:,swarmIndAvoiding);
                    
                    deltaSMax(:,:,swarmIndStop)=5;
                    deltaSMin(:,:,swarmIndStop)=5;
                    kPre2 = zeros(1, 1, length(swarmIndStop));
                    tPre2 = zeros(1, 1, length(swarmIndStop));
                    omegaPre2 = zeros(1, 1, length(swarmIndStop));
                    %             omegaPre2(1,1,:)=atand(vecIndivForm(2,swarmIndStop)./vecIndivForm(1,swarmIndStop));
                    angTemp=atand(vecIndivFormUni(2,swarmIndStop)./vecIndivFormUni(1,swarmIndStop));
                    angTemp(vecIndivFormUni(1,swarmIndStop)<0)=180+angTemp(vecIndivFormUni(1,swarmIndStop)<0);
                    omegaPre2(1,1,:)=angTemp;
                    %             [arcx2, arcy2, x1Pre2, y1Pre2]=arcParameterization ...
                    %                 (x0(:,:,swarmIndStop), y0(:,:,swarmIndStop), kPre2, omegaPre2, deltaSMax(:,:,swarmIndStop), N);
                    
                    [arcx2, arcy2, arcAlt2, x1Pre2, y1Pre2]=arcParameterization3D ...
                        (x0(:,:,swarmIndStop), y0(:,:,swarmIndStop), Alt0(:,:,swarmIndStop), kPre2, tPre2, omegaPre2, deltaSMax(:,:,swarmIndStop), N);
                    
                    arcx(:,1,swarmIndStop)=arcx2;
                    arcy(:,1,swarmIndStop)=arcy2;
                    arcz(:,1,swarmIndStop)=ones(size(arcy2))*100;
                    arcAlt(:,1,swarmIndStop)=arcAlt2;
                    omegaPre (:,:,swarmIndStop) = omegaPre2;
                    
                    kPre2 = zeros(1, 1, length(swarmIndNonAvoiding));
                    tPre2 = zeros(1, 1, length(swarmIndNonAvoiding));
                    omegaPre2 = zeros(1, 1, length(swarmIndNonAvoiding));
                    omegaPre2(1,1,:)=atand(vecSwarm2DestUni(2)./vecSwarm2DestUni(1));
                    %             [arcx2, arcy2, x1Pre, y1Pre]=arcParameterization ...
                    %                 (x0(:,:,swarmIndNonAvoiding), y0(:,:,swarmIndNonAvoiding), kPre2, omegaPre2, deltaSMax(:,:,swarmIndNonAvoiding), N);
                    
                    [arcx2, arcy2, arcAlt2, x1Pre, y1Pre]=arcParameterization3D ...
                        (x0(:,:,swarmIndNonAvoiding), y0(:,:,swarmIndNonAvoiding), Alt0(:,:,swarmIndNonAvoiding), ...
                        kPre2, tPre2, omegaPre2, deltaSMax(:,:,swarmIndNonAvoiding), N);
                    
                    arcx(:,1,swarmIndNonAvoiding)=arcx2;
                    arcy(:,1,swarmIndNonAvoiding)=arcy2;
                    arcz(:,1,swarmIndNonAvoiding)=ones(size(arcy2))*500;
                    arcAlt(:,1,swarmIndNonAvoiding)=arcAlt2;
                    omegaPre (:,:,swarmIndNonAvoiding) = omegaPre2;
                    
                    x0temp=x0(:,:,vecnorm(vecIndivForm)<=5);
                    y0temp=y0(:,:,vecnorm(vecIndivForm)<=5);
                    Alt0temp=Alt0(:,:,vecnorm(vecIndivForm)<=5);
                    arcx(:,1, vecnorm(vecIndivForm)<=5)=x0temp(ones(N,1),:,:);
                    arcy(:,1, vecnorm(vecIndivForm)<=5)=y0temp(ones(N,1),:,:);
                    arcz(:,1, vecnorm(vecIndivForm)<=5)=ones(size(arcy(:,1, vecnorm(vecIndivForm)<=5)))*100;
                    arcAlt(:,1, vecnorm(vecIndivForm)<=5)=Alt0temp(ones(N,1),:,:);
                end
            end
            
            timeCost=toc;
            %     trajectSet(end,4,:)=timeCost(:,:,ones(1, swarmSize));
            trajectTemp=cat(2, arcx, arcy, arcAlt, timeCost(ones(size(arcz))), arcz);
            trajectSet=cat(1,trajectSet,trajectTemp);
            trajectAvoiding=cat(1, trajectAvoiding, [trajectTemp, ind_t(ones(1,N),1,ones(1,swarmSize))]);
            %     if plotFlag
            %         fprintf('\n');
            %     end
            
            arcxPre=trajectSet(end-N+1:end, 1, :);
            arcyPre=trajectSet(end-N+1:end, 2, :);
            arcAltPre=trajectSet(end-N+1:end, 3, :);
            x0=arcxPre(end,:,:);
            y0=arcyPre(end,:,:);
            Alt0=arcAltPre(end,:,:);
            
            posIndiv = [x0(:)'; y0(:)'];
            posIndivAlt = Alt0(:)';
            
            dist2intruTemp=vecnorm(posIntruRel-posIndiv, 2);
            dist2intru=cat(1, dist2intru, dist2intruTemp);
            
            distIntruIndiv = sqrt(sum((posIndiv(:,:,ones(intruNum, 1))-reshape(posIntruRel, 2, 1, intruNum)).^2,1));
            distIntruIndiv = min(distIntruIndiv, [], 3);
            distIndivDest = sqrt(sum((posIndiv-posIniDest).^2,1));
            
            vecIntruIndiv = reshape(posIntruRel, 2, 1, intruNum) - posIndiv(:,:,ones(intruNum, 1));
            vecIntruIndivUni = vecIntruIndiv./ vecnorm(vecIntruIndiv);
            %     vecDestIndiv = posIniDest - posIndiv;
            vecDestIndiv = [350; 350] - posIniIndiv;
            vecDestIndivUni = vecDestIndiv./ vecnorm(vecDestIndiv);
            
            stateIndicator=dot(vecIntruIndivUni, vecDestIndivUni(:,:,ones(1, intruNum)));
            stateIndicator (stateIndicator >= 0) = 0;
            stateIndicator (stateIndicator < 0) = -1;
            stateIndicator = sum(stateIndicator, 3);
            stateIndicator (stateIndicator == 0) = 1;
            
            if isempty(stateIndicator (stateIndicator>0))
                stateIndicator = -ones(1, swarmSize);
            end
            if plotFlag
                %% Plotting %%
                pause(.01);
                
                [contourPlot, sectorNum] = ContourPlot (smoothField, z0);
                contourPlotCell{ind_t}=contourPlot;
                g4.XData=contourPlot(:,1);
                g4.YData=contourPlot(:,2);
                g4.ZData=contourPlot(:,3);
                %     smoothField(smoothField>min(zForbid))=0;
                %     g3.ZData=smoothField;
                %     g5.ZData=max(smoothField(:));
                g1.XData=x0;
                g1.YData=y0;
                g1.ZData=Alt0+max(smoothField(:));
                g2.XData=posIntruRel(1,:);
                g2.YData=posIntruRel(2,:);
                g2.ZData=interp2(smoothField, posIntruRel(1,:), posIntruRel(2,:));
                for ind_h=1: swarmSize
                    eval(['h',num2str(ind_h),'.XData=trajectSet(:,1,ind_h);']);
                    eval(['h',num2str(ind_h),'.YData=trajectSet(:,2,ind_h);']);
                    eval(['h',num2str(ind_h),'.ZData=trajectSet(:,3,ind_h);']);
                    %         eval(['h',num2str(ind_h),'.ZData=ones(size(trajectSet(:,3,ind_h)))*max(smoothField(:));']);
                end
                for ind_k=1: swarmSize
                    eval(['k',num2str(ind_k),'.XData=arcPreOrigin(:,1,ind_k);']);
                    eval(['k',num2str(ind_k),'.YData=arcPreOrigin(:,2,ind_k);']);
                    eval(['k',num2str(ind_k),'.ZData=arcPreOrigin(:,3,ind_k);']);
                    %         eval(['k',num2str(ind_k),'.ZData=ones(size(arcPreOrigin(:,3,ind_k)))*max(smoothField(:));']);
                end
                for ind_p=1: swarmSize
                    eval(['p',num2str(ind_p),'.XData=xPred(:,1,ind_p);']);
                    eval(['p',num2str(ind_p),'.YData=yPred(:,1,ind_p);']);
                    eval(['p',num2str(ind_p),'.ZData=ones(size(xPred(:,1,ind_p)))*max(smoothField(:));']);
                end
                
                eval(['title(axe2,''Iteration ',num2str(ind_t),''');' ]);
                
                posSwarmSet=cat(2, posSwarmSet, posSwarmRel);
                
                frame = getframe(fig2);
                writeVideo(writerObj1,frame);
            end
        end
        endX=arcx(end,:,:);
        endY=arcy(end,:,:);
        endZ=arcz(end,:,:);
        arcAlt(end,:,:)=posIniIndivAlt;
        endAlt=arcAlt(end,:,:);
        % trajectTemp=cat(2, endX(ones(1, N),:,:), endY(ones(1, N),:,:), endAlt(ones(1, N),:,:), timeCost(ones(size(arcz))));
        trajectTemp=cat(2, arcx(end,:,:), arcy(end,:,:), arcAlt(end,:,:), timeCost(ones(size(arcAlt(end,:,:)))), arcz(end,:,:));
        trajectSet=cat(1,trajectSet,trajectTemp);
        for ind_h=1: swarmSize
            eval(['h',num2str(ind_h),'.XData=trajectSet(:,1,ind_h);']);
            eval(['h',num2str(ind_h),'.YData=trajectSet(:,2,ind_h);']);
            eval(['h',num2str(ind_h),'.ZData=trajectSet(:,3,ind_h);']);
            %         eval(['h',num2str(ind_h),'.ZData=ones(size(trajectSet(:,3,ind_h)))*max(smoothField(:));']);
        end
        for ind_k=1: swarmSize
            eval(['k',num2str(ind_k),'.XData=arcPreOrigin(:,1,ind_k);']);
            eval(['k',num2str(ind_k),'.YData=arcPreOrigin(:,2,ind_k);']);
            eval(['k',num2str(ind_k),'.ZData=arcPreOrigin(:,3,ind_k);']);
            %         eval(['k',num2str(ind_k),'.ZData=ones(size(arcPreOrigin(:,3,ind_k)))*max(smoothField(:));']);
        end
        pause (1);
        if plotFlag
            hold (axe2);
            close(writerObj1);
        end
    end
end

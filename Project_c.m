%% Introduction to Computational Materials Engineering
% Final Project 8.24 - task c
% Written by Ofir Friedman 300483153
% Written in Matlab R2013a
% ofirfr@post.bgu.ac.il

%% initial function
% recives first parameters pf the simulation
% sets the world size accordingly and runs the simulation.
function Project_c(btnDimentionsValue,btnValParticles,btnValInitialTemp)
  
 constSigma = 3.4*10^-10;
 
% Setting basic parameters and triangular lattice:
    TotalTime=2*10^-10;
    worldsizeX=constSigma*8 * btnValParticles/8+constSigma*5/2;
    worldsizeY=worldsizeX*(sqrt(3)/2);

    if btnDimentionsValue==1
        TotalTime=10^-9;
        worldsizeX=constSigma*8*btnValParticles+constSigma*5/2;
    end 
 Run_Model(btnValParticles,btnDimentionsValue,TotalTime,worldsizeX,worldsizeY, btnValInitialTemp);
   
end

%% Model Generator
function Run_Model(numParticles,Dimentions,time,WorldSizeX,WorldSizeY, InitialTemperature)


     %Variables
    global constSigma;
    global constE;
    constE = 1.65*10^-21;
    constSigma = 3.4*10^-10;
  
    WorldSize=[WorldSizeX WorldSizeY];
    M(1:numParticles,1)=6.67*10^-26;
    
    if Dimentions==1
        TimeLap=(constSigma*sqrt((M(1,1)/constE)))/40;
    else
        TimeLap=(constSigma*sqrt((M(1,1)/constE)))/40;
    end
    
    t=0:TimeLap:time-TimeLap;
    steps=floor(time/TimeLap);
    ParticlesPositions=zeros(numParticles,Dimentions,steps);
    F=zeros(numParticles,Dimentions,steps);
    
    %The key new command and function that distinguishes this task from the
    %original molecular dynamics

    V(:,:,1)=SetInitialVelocities(InitialTemperature,numParticles,Dimentions);
    
    % Instantaneous temperature
    temperature=zeros(steps,1);
    temperature(1)=CalculateTemperature(numParticles,Dimentions,V,1);
        
    % 1=Ek, 2=Ep
    Etot= zeros(2,steps);
    
    %Set Starting Positions
    if Dimentions==1
                
        ParticlesPositions(1,1,1)=constSigma*8;
        for k=2:numParticles
         ParticlesPositions(k,1,1)=ParticlesPositions(k-1,1,1)+constSigma*8;
        end
    end
    
    if Dimentions~=1
        
       ParticlesPositions(1,1,1)=constSigma*8*sqrt(3)/2;
       ParticlesPositions(1,2,1)=constSigma*8*sqrt(3)/2;

       flagStep=0;
       floarLevel=1;
       flagReset=0;
        for k=2:numParticles

            if mod(k,8)==0 && flagStep==0
                
                ParticlesPositions(k,1,1)=ParticlesPositions(k-1,1,1)+constSigma*8;
                ParticlesPositions(k,2,1)=constSigma*8*(sqrt(3)/2)*floarLevel;
                flagStep=1;
                floarLevel=floarLevel+1;
                flagReset=1;
                
            elseif mod(k,8)==0 && flagStep~=0
            
                ParticlesPositions(k,1,1)=ParticlesPositions(k-1,1,1)+constSigma*8;
                ParticlesPositions(k,2,1)=constSigma*8*(sqrt(3)/2)*floarLevel;
                flagStep=0;
                floarLevel=floarLevel+1;
                flagReset=1;
            
            else
                
            if flagReset==1 && flagStep==0
                ParticlesPositions(k,1,1)=ParticlesPositions(1,1,1);
                flagReset=0;
                
            elseif flagReset==1 && flagStep==1
                 ParticlesPositions(k,1,1)=ParticlesPositions(1,1,1)+constSigma*4;
                flagReset=0;               
            else
                ParticlesPositions(k,1,1)=ParticlesPositions(k-1,1,1)+constSigma*8;
            end
            
            ParticlesPositions(k,2,1)=constSigma*8*(sqrt(3)/2)*floarLevel;

            end
        end

    end
    
    
    hCalc = waitbar(0,'Calculating...');
   
%    Generate Steps
     for j=1:Dimentions
            for k=1:numParticles

                F(k,j,1)=LennardJonesForceCalculation(k,1,ParticlesPositions,WorldSize,numParticles,Dimentions,j);
                
            end
            
     end
        
    % Calculate initial Energies 
   [Etot(1,1), Etot(2,1)] = CalcE(Dimentions,numParticles,M,V,1,ParticlesPositions,WorldSize);
   
   %Main VERLET EOM
    for i=2:steps
        
        for j=1:Dimentions
            for k=1:numParticles

                ParticlesPositions(k,j,i) = VerletAlgorithm(ParticlesPositions(k,j,i-1),V(k,j,i-1),TimeLap,WorldSize,j,F(k,j,i-1)/M(k));

            end
        end

        for j=1:Dimentions
            for k=1:numParticles

                F(k,j,i)=LennardJonesForceCalculation(k,i,ParticlesPositions,WorldSize,numParticles,Dimentions,j);

                V(k,j,i)=V(k,j,i-1)+0.5*((F(k,j,i)/M(k))+F(k,j,i-1)/M(k))*TimeLap;

            end
        end
          
            waitbar(i/steps ,hCalc,sprintf('Currect Step: %d , Total Steps: %d',i, steps));

            V(:,:,i)=RescaleVelocities(InitialTemperature,V,i,numParticles,Dimentions);

            [Etot(1,i), Etot(2,i)] = CalcE(Dimentions,numParticles,M,V,i,ParticlesPositions,WorldSize);

            temperature(i)=CalculateTemperature(numParticles,Dimentions,V,i);

    end
    
    close(hCalc);
   
    % Plot
    %Movment
    if Dimentions==1
        
        figure(2);
         set(2,'Position',[370 150 500 500]);
        
         for i=1:numParticles

            plot(t,squeeze(ParticlesPositions(i,1,:)),'Color',[rand, rand, rand]);
            hold on;
            xlabel('t');
            ylabel('X(t)');
            grid on;
            ylim([0 WorldSizeX]);
    
         end
        
    else
        
        figure(2);
        set(2,'Position',[370 150 500 500]);
        
         for i=1:numParticles
            
          
            plot(squeeze(ParticlesPositions(i,1,:)),squeeze(ParticlesPositions(i,2,:)),'Marker','.','Color','r');
            hold on;
            plot(squeeze(ParticlesPositions(i,1,1)),squeeze(ParticlesPositions(i,2,1)),'Marker','.','Color','b');
            hold on;
            plot(squeeze(ParticlesPositions(i,1,steps)),squeeze(ParticlesPositions(i,2,steps)),'Marker','.','Color','g');
            hold on;
            xlabel('X(t)');
            ylabel('Y(t)');
        

            answerRed = (['route']);
            answerBlue = (['Start']);
            answerGreen = (['End']);
            legend(answerRed,answerBlue,answerGreen);
            grid on;
            xlim([0 WorldSizeX]);
            ylim([0 WorldSizeY]);
    
         end
    
    end
    
  
    for i=1:steps
        
        Esum(i)=Etot(1,i)+Etot(2,i);
        
    end
    
    %energies
    figure(3);
    
    subplot(3,1,1);
    plot(t,Etot(1,:));
    answerEk = (['Ek']);
    legend(answerEk);
     xlabel('t');
    ylabel('E');
    grid on;
    
    subplot(3,1,2);
    plot(t,Etot(2,:),'color','red');
    answerEp = (['Ep']);
    legend(answerEp);
    xlabel('t');
    ylabel('E');
    grid on;
    
    subplot(3,1,3);
    plot(t,Esum,'color','green');
    answerEsum = (['E-Total']);
    legend(answerEsum);
    xlabel('t');
    ylabel('E');
    grid on;
    
    %Plot temperature change in time
    
    figure (4);
    plot(t,temperature(:));
    xlabel('t');
    ylabel('T(t)');
    grid on;

     %export data to excle
     
    Table_temperature = 'temperature.xlsx';
    xlswrite(Table_temperature,temperature(:));

    Table_Velocities_x = 'Velocities_x.xlsx';
    xlswrite(Table_Velocities_x,squeeze(V(:,1,:)));

    Table_Masses_x = 'Masses_x.xlsx';
    xlswrite(Table_Masses_x,M(:,1));

    Table_Forces_x = 'Forces_x.xlsx';
    xlswrite(Table_Forces_x,squeeze(F(:,1,:)));

    Table_Positions_x = 'Positions_x.xlsx';
    xlswrite(Table_Positions_x,squeeze(ParticlesPositions(:,1,:)));

    Table_Energy_x = 'Energy_x.xlsx';
    xlswrite(Table_Energy_x,Etot(:,:));

    if Dimentions~=1
        Table_Velocities_y = 'Velocities_y.xlsx';
        xlswrite(Table_Velocities_y,squeeze(V(:,2,:)));

        Table_Forces_y = 'Forces_y.xlsx';
        xlswrite(Table_Forces_y,squeeze(F(:,2,:)));

        Table_Positions_y = 'Positions_y.xlsx';
        xlswrite(Table_Positions_y,squeeze(ParticlesPositions(:,2,:)));
    end
    
end


%% Check the current temperature
% compare it to the initial temperature
% and rescale all velocities
function balancedVelocity = RescaleVelocities(InitialTemperature,V,currentStep,numParticles,Dimentions)

    Velocity=zeros(numParticles,Dimentions);

    CurrentTemperature = CalculateTemperature(numParticles,Dimentions,V,currentStep);

    RescaleVelocityFactor = InitialTemperature/CurrentTemperature;

    for i=1:numParticles

        Velocity(i,1)= V(i,1,currentStep)*RescaleVelocityFactor;

        if Dimentions>1
         Velocity(i,2)= V(i,2,currentStep)*RescaleVelocityFactor;
        end

    end
    
 balancedVelocity=Velocity;
 
end


%% Calculate the current temperature
% relaing on the particles current speeds
function currentT = CalculateTemperature(numParticles,Dimentions,Velocity, currentStep)

    MassVelocity=0;

    for i=1:numParticles

        
        if Dimentions==1
         MassVelocity= MassVelocity +  Velocity(i,currentStep)* Velocity(i,currentStep);
        
        elseif Dimentions>1
             MassVelocity= MassVelocity +  Velocity(i,1,currentStep)* Velocity(i,1,currentStep);
             MassVelocity= MassVelocity + Velocity(i,2,currentStep)* Velocity(i,2,currentStep);
        end

    end
    
    if Dimentions==1
        currentT =  MassVelocity(1)/(numParticles*Dimentions*2);

    elseif Dimentions>1
        currentT =  MassVelocity(1)/(numParticles*Dimentions);
    end

end

%%Recieve  initial temperature and set initial velocities to the particles 
%initial velocities are based on the equation:
% Ek*T(t)=sum(M*V^2)/N*d
%M=mass,V=velocity, N=total particles, d=amount of dimentions
function balancedVelocity = SetInitialVelocities(InitialTemperature, numParticles,Dimentions)

    Velocity=zeros(numParticles,Dimentions);
    VxTotal=0;
    VyTotal=0;

    % assign randon initial velocities
    for i=1:numParticles
        
        Velocity(i,1) = rand() - 0.5;
        
        if Dimentions>1
            Velocity(i,2) = rand() - 0.5;
        end
        
        VxTotal=VxTotal+Velocity(i,1);
        
         if Dimentions>1
        VyTotal=VyTotal+ Velocity(i,2);
         end
        
    end
    
    %zero center of mass momentum
    VxCentered=VxTotal/numParticles;
    
     if Dimentions>1
        VyCentered=VyTotal/numParticles;
     end
     
    for i=1:numParticles
        
        Velocity(i,1)= Velocity(i,1) - VxCentered;
        
        if Dimentions>1
         Velocity(i,2)= Velocity(i,2) - VyCentered;
        end
         
    end
    
    Vsum=0;
    
    %calculate sum of velocities
    for i=1:numParticles

        Vsum=Vsum +  Velocity(i,1)* Velocity(i,1);
        
         if Dimentions>1
             Vsum=Vsum+ Velocity(i,2)* Velocity(i,2);
         end

    end
    
    %find total kinetic energy and set the rescale factor
   Ek=0.5*Vsum/numParticles;
   RescaleVelocityFactor = sqrt(InitialTemperature/Ek);
   
   %update the velocities
    for i=1:numParticles
        
        Velocity(i,1)= Velocity(i,1)*RescaleVelocityFactor;
        
        if Dimentions>1
            Velocity(i,2)= Velocity(i,2)*RescaleVelocityFactor;
        end
        
    end
     
    balancedVelocity=Velocity;
    
end

%% Calculate total Ek
% and Ep in the system
function [EKr, EPr] = CalcE(Dimentions,numParticles,m,v,stepNumber,ParticlesPositions,WorldSize)

    Ecalc=zeros(Dimentions,numParticles);
    Ep=zeros(1,numParticles);
    EKtot=0;
    Eptot=0;
    Estep=1;
    
    myESP=1.65*10^-21;
    mySIG=3.4*10^-10;
    
    
         for i=1:Dimentions
               for j=1:numParticles

                    Ecalc(i,j)=0.5*(m(j)*(v(j,i,stepNumber)^2));

               end
         end
         
          for i=1:Dimentions
               for j=1:1
                    for k=1:numParticles
                        if k~=j

                                Directions=vectorDirection(k,j, stepNumber,ParticlesPositions,WorldSize,Dimentions);


                            if Dimentions==1

                                r=CalcDistance(Directions+2,0,k,j,stepNumber,ParticlesPositions,WorldSize,Dimentions);

                            else

                                r=CalcDistance(Directions(1)+2,Directions(2)+2,j,k,stepNumber,ParticlesPositions,WorldSize,Dimentions);

                            end

                            Ep(j)=4*myESP*((mySIG/r)^12-(mySIG/r)^6);

                        end
                    end

                     Eptot = Eptot+Ep(j);

               end
         end
         

   
    if Dimentions==1

        for i=1:numParticles
                
             EKtot=EKtot+ Ecalc(1,i);
              
        end
    
    end
     
     if Dimentions==2
         
         for i=1:numParticles
             
             Estep=0;
             
           for j=1:Dimentions
               
                Estep=Estep+Ecalc(j,i);
                
           end
           
           EKtot=EKtot+Estep;
           
          end
         
     end
     
  
     EKr = EKtot;
     EPr = Eptot;
     
     
end

%% Calculate next position with Verlet Algorithm
function newPosition = VerletAlgorithm(LastPosition,CurrentVelocity,TimeLap,WorldSize,CurrentDimention,oldAcceleration)


    calculatednewPosition=LastPosition+CurrentVelocity*TimeLap+0.5*oldAcceleration*TimeLap^2;
    
    
    if calculatednewPosition<WorldSize(CurrentDimention) && calculatednewPosition>0
        
        newPosition = calculatednewPosition;
    
    elseif calculatednewPosition>WorldSize(CurrentDimention)
        
          calculatednewPosition= mod(calculatednewPosition,WorldSize(CurrentDimention));
        
    else
        
         calculatednewPosition=WorldSize(CurrentDimention)-mod(WorldSize(CurrentDimention),calculatednewPosition);
         
    end
    
    newPosition = calculatednewPosition;
    
end

%% Calculate the total forces applied on a particle from all other particles in the box
function clostestForceTotal = LennardJonesForceCalculation(particleNumber, stepNumber,ParticlesPositions,WorldSize,numParticles,Dimentions,currentDomain)

    totalForce=0;
    
    % direction size is as nubmer of dimentions, 1 for positive direction and
    % for negative directions the value is -1 per dimention. 
    % 0 for original spot.
    
    Directions=zeros(Dimentions);
 
    for k=1:numParticles
        if k~=particleNumber
         
                Directions=vectorDirection(k, particleNumber, stepNumber,ParticlesPositions,WorldSize,Dimentions);
                totalForce=totalForce+canculateSingleForce(k, particleNumber,Directions,stepNumber,ParticlesPositions,WorldSize,Dimentions,currentDomain);

        end
    end
        
    
        clostestForceTotal=totalForce;
        
  

end

%% Find the closest direction between two particles
% This function sets the initial search parameter
% and runs associated functions
function directions = vectorDirection(sourceParticle,targetParticle,step,ParticlesPositions,WorldSize,Dimentions)

    % building arrays of flags for the required dimentions set
    %1D:
    if Dimentions==1
        
        % 3 options: 1,0,-1 and set of flags to note the closest option
       distances= zeros(3);
       
       for i=1:3
           distances(i)=i-2;
       end
        
    end
    
    
    %2D:
     if Dimentions==2
        
        % 9 options: 1,0,-1 for both X and Y axis
       distances= zeros(3,3,2);
       
       for i=1:3
           for j=1:3
            distances(i,j,1)=i-2;
            distances(i,j,2)=j-2;
           end
       end
        
     end    

     
     % Call Distance calculators and return the closest direction
    if Dimentions==1
        closestFlagLocation = recursionFlagCloseDistance(1,0,sourceParticle,targetParticle,step,ParticlesPositions,WorldSize,Dimentions);
      %  disp(closestFlagLocation);
        directions=distances(closestFlagLocation);
        
    end
    
    if Dimentions==2
        closestFlagLocation = recursionFlagCloseDistance(1,1,sourceParticle,targetParticle,step,ParticlesPositions,WorldSize,Dimentions);
        directions=[distances(closestFlagLocation(1),closestFlagLocation(2),1) distances(closestFlagLocation(1),closestFlagLocation(2),2)];
    end
    
    
    
end

%% Recursive function which return the closest directions; 
function flaggedDistance = recursionFlagCloseDistance(PositionX,PositionY,sourceParticle,targetParticle,step,ParticlesPositions,WorldSize,Dimentions)
 
global endFlag;
endFlag=0;

    %1D:
    if Dimentions==1

        if PositionX<3 && endFlag==0
            
            closestPosition = recursionFlagCloseDistance(PositionX+1,0,sourceParticle,targetParticle,step,ParticlesPositions,WorldSize,Dimentions);
            
        end
        
        if PositionX==3
        
            endFlag=1;
            flaggedDistance=3;
            
        end
        
         if PositionX<3 && endFlag==1
             
                if  CalcDistance(PositionX,0,sourceParticle,targetParticle,step,ParticlesPositions,WorldSize,Dimentions)< CalcDistance(closestPosition,0,sourceParticle,targetParticle,step,ParticlesPositions,WorldSize,Dimentions)
            
                     flaggedDistance=PositionX;
                    
                
                elseif PositionX==0 && CalcDistance(PositionX,0,sourceParticle,targetParticle,step,ParticlesPositions,WorldSize,Dimentions)== CalcDistance(closestPosition,0,sourceParticle,targetParticle,step,ParticlesPositions,WorldSize,Dimentions)
                    
                     flaggedDistance=PositionX;
                    
                else
                    
                    flaggedDistance=closestPosition;
             
                end
             
             
         end
         
    end
    
    
    %2D:
    if Dimentions==2
        
        %opens all cells while end isnt reached
        if PositionX<=3 &&  PositionY~=3 && endFlag==0

            closestPosition = recursionFlagCloseDistance(PositionX,PositionY+1,sourceParticle,targetParticle,step,ParticlesPositions,WorldSize,Dimentions);
   
        end
        
        if PositionX<3 &&  PositionY==3 && endFlag==0
            
                 closestPosition = recursionFlagCloseDistance(PositionX+1,1,sourceParticle,targetParticle,step,ParticlesPositions,WorldSize,Dimentions);
    
        end

        %return of the last cell
        if PositionX==3  &&  PositionY==3 

            endFlag=1;
            flaggedDistance=[3 3];
             return;
        end

        % closing and comparing values
        if  endFlag==1

            if  CalcDistance(PositionX,PositionY,sourceParticle,targetParticle,step,ParticlesPositions,WorldSize,Dimentions) < CalcDistance(closestPosition(1),closestPosition(2),sourceParticle,targetParticle,step,ParticlesPositions,WorldSize,Dimentions)

                 flaggedDistance=[PositionX PositionY];


            elseif PositionX==0 && PositionY == 0 && CalcDistance(PositionX,PositionY,sourceParticle,targetParticle,step,ParticlesPositions,WorldSize,Dimentions) == CalcDistance(closestPosition(1),closestPosition(2),sourceParticle,targetParticle,step,ParticlesPositions,WorldSize,Dimentions)

                 flaggedDistance=[PositionX PositionY];

                
            else

                flaggedDistance=closestPosition;

            end


        end
         
     end
     
    
end

 %% Calculate Distance between two particles
function distance = CalcDistance(PositionX,PositionY,sourceParticle,targetParticle,step,ParticlesPositions,WorldSize,Dimentions)


  %1D:
    if Dimentions==1

        if (PositionX-2)~=0
            
            distance = abs((PositionX-2)*WorldSize(1)+ParticlesPositions(targetParticle,1,step)-ParticlesPositions(sourceParticle,1,step));
        
        else
            
            distance = abs(ParticlesPositions(targetParticle,1,step)-ParticlesPositions(sourceParticle,1,step));
            
        end
    end
    
    
    %2D:
     if Dimentions==2
 
         if (PositionX-2)~=0&&(PositionY-2)~=0
             
            distanceX=abs((PositionX-2)*WorldSize(1)+ParticlesPositions(targetParticle,1,step)-ParticlesPositions(sourceParticle,1,step));
            distanceY=abs((PositionY-2)*WorldSize(2)+ParticlesPositions(targetParticle,2,step)-ParticlesPositions(sourceParticle,2,step));

         elseif (PositionX-2)~=0&&(PositionY-2)==0
         
            distanceX=abs((PositionX-2)*WorldSize(1)+ParticlesPositions(targetParticle,1,step)-ParticlesPositions(sourceParticle,1,step));
            distanceY=abs(ParticlesPositions(targetParticle,2,step)-ParticlesPositions(sourceParticle,2,step));
                
        elseif (PositionX-2)==0&&(PositionY-2)~=0
         
            distanceX=abs(ParticlesPositions(targetParticle,1,step)-ParticlesPositions(sourceParticle,1,step));
            distanceY=abs((PositionY-2)*WorldSize(2)+ParticlesPositions(targetParticle,2,step)-ParticlesPositions(sourceParticle,2,step));

         else
             
            distanceX=abs(ParticlesPositions(targetParticle,1,step)-ParticlesPositions(sourceParticle,1,step));
            distanceY=abs(ParticlesPositions(targetParticle,2,step)-ParticlesPositions(sourceParticle,2,step));
            
         end
         
         if distanceX==0 && distanceY~=0
             
             distance=abs(distanceY);
             
         elseif distanceX~=0 && distanceY==0
             
             distance=abs(distanceX);
             
         else
             
            distance=sqrt(distanceX^2+distanceY^2);
            
         end
         
     end

end

%% Calculate Force between two particles given calculation direction
function singleForceVector = canculateSingleForce(sourceParticle,targetParticle,Directions,step,ParticlesPositions,WorldSize,Dimentions,currentDomain)

    constE = 1.65*10^-21;
    constSigma = 3.4*10^-10;

        if length(Directions)==1
            if Directions~=0

                r = Directions*WorldSize(1)+ParticlesPositions(targetParticle,1,step)-ParticlesPositions(sourceParticle,1,step);

            else

                r = ParticlesPositions(targetParticle,1,step)-ParticlesPositions(sourceParticle,1,step);

            end

            singleForceVector=((24*constE)/(r))*(2*(((constSigma)/(r))^(12))-((constSigma)/(r))^(6));
            
        else
            
             if Directions(1)~=0&&Directions(2)~=0

               distanceX=Directions(1)*WorldSize(1)+ParticlesPositions(targetParticle,1,step)-ParticlesPositions(sourceParticle,1,step);
               distanceY=Directions(2)*WorldSize(2)+ParticlesPositions(targetParticle,2,step)-ParticlesPositions(sourceParticle,2,step);

            elseif Directions(1)~=0&&Directions(2)==0

                distanceX=Directions(1)*WorldSize(1)+ParticlesPositions(targetParticle,1,step)-ParticlesPositions(sourceParticle,1,step);
               distanceY=ParticlesPositions(targetParticle,2,step)-ParticlesPositions(sourceParticle,2,step);

             elseif Directions(1)==0&&Directions(2)~=0

               distanceX=ParticlesPositions(targetParticle,1,step)-ParticlesPositions(sourceParticle,1,step);
               distanceY=Directions(2)*WorldSize(2)+ParticlesPositions(targetParticle,2,step)-ParticlesPositions(sourceParticle,2,step);

            else

              distanceX=ParticlesPositions(targetParticle,1,step)-ParticlesPositions(sourceParticle,1,step);
              distanceY=ParticlesPositions(targetParticle,2,step)-ParticlesPositions(sourceParticle,2,step);

            end
            
            if distanceX==0 && distanceY~=0

                r=abs(distanceY);

                if currentDomain==1

                    singleForceSize=((24*constE)/(r))*(2*(((constSigma)/(r))^(12))-((constSigma)/(r))^(6));

                    singleForceVector=0;

                else

                    singleForceSize=((24*constE)/(r))*(2*(((constSigma)/(r))^(12))-((constSigma)/(r))^(6));

                    singleForceVector=singleForceSize*r;

                end

            elseif distanceX~=0 && distanceY==0

                r=abs(distanceX);

                if currentDomain==1

                    singleForceSize=((24*constE)/(r))*(2*(((constSigma)/(r))^(12))-((constSigma)/(r))^(6));

                    singleForceVector=singleForceSize*r;

                else

                    singleForceSize=((24*constE)/(r))*(2*(((constSigma)/(r))^(12))-((constSigma)/(r))^(6));

                    singleForceVector=0;

                end

            else

                r=sqrt(distanceX^2+distanceY^2);

                if currentDomain==1

                    singleForceSize=((24*constE)/(r))*(2*(((constSigma)/(r))^(12))-((constSigma)/(r))^(6));

                    singleForceVector=singleForceSize*(distanceX/r);

                else

                    singleForceSize=((24*constE)/(r))*(2*(((constSigma)/(r))^(12))-((constSigma)/(r))^(6));

                    singleForceVector=singleForceSize*(distanceY/r);

                end

            end
        end
    
end

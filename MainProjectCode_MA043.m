%% Clean up
clc
clear 
close all
%% Input and input restrictions 
[crane_rgb] = imread('Crane.JPG'); %Reads the crane's image into the workspace. (Notes: This built-in function wasn't covered in class)
imshow(crane_rgb)%Displays the crane's image to the user

disp('Please refer to the crane model shown and enter the values below:') %Displaying this image minimizes data entry mistakes. 

M1=input('Please enter mass 1 of the main truck body in Mg:'); %Mg
M2=input('Please enter mass 2 of the boom in Mg:'); %Mg
M3=input('Please enter mass 3 of the hanging weight in Mg:'); %Mg
W1=M1*9.81*10^3; %N , multiplies the mass by (10^3) to turn it into Kg then by the mass due to gravity (9.8) to turn into Newtons.
W2=M2*9.81*10^3; %N
W3=M3*9.81*10^3; %N

L1=input('Please enter a vecotr of lengths for L1 in meters such that the vector contains five elements:'); %m , allows the user to enter five lengths for L1 in order to test them agains five lengths for L2
L2=input('Please enter a vecotr of lengths for L2 in meters such that the vector contains five elements:'); %m , allows the user to enter five lengths for L2 in order to test them agains five lengths for L1

while length(L1)~=length(L2)||length(L1)~=5 %This while loop will keep prompting the user to enter L1 and L2 until both L1 and L2 contain five elements each. 
    fprintf('Pleas enter the lengths L1 and L2 such that they each have five elements.\n')
    L1=input('Please enter a vercotr of lengths for L1 in meters such that the vector contains five elements::');%m
    L2=input('Please enter a vercotr of lengths for L2 in meters such that the vector contains five elements:');%m
end

%for i=1:length(L1) %This for loop runs through each element of L1 and L2 to check them against each other through the nested while loop.
  % while L2(i)/L1(i)~=0.96 %This while loop will keep prompting the user to enter L1 and L2 until the condition is satisfied. This condition is to ensure that the lengths to the center of gravity are entered correctly by the user. Note that the accuracy has to be up to 9 decimal places. 
       % fprintf('Pleas enter the lengths L1 and L2 such that L2/L1=0.96 for each two corresponding lengths L1 and L2.\n')
       % L1=input('Please enter a vercotr of lengths for L1:');%m
       % L2=input('Please enter a vercotr of lengths for L2:');%m
  % end
%end

L3=input('Please enter length L3 in meters:'); %m
L4=input('Please enter length L4 in meters:'); %m
L5=input('Please enter length L5 in meters:'); %m

while L5~=(1/4)*(L3+L4+L5)%This while loop will keep prompting the user to enter L5 until the condition is satisfied. This loop is to ensure that the length L5 to the center of gravity is entered correctly by the user.
    fprintf('Pleas enter length 5 such that the center of gravity of the truck body is 1/4 the length of AB from B.\n')
    L5=input('Please enter length 5 between outrigger B and the center of gravity of the truck body:'); %m
end 

L7=input('Please enter length L7 in meters:'); %m
theta=input('Please enter the desired theta (in degrees) at which the critical angle phi should be calculated:'); %degrees , this will be the theta used to calculte the critical angles phiB and phiD
desiredPhi=input('Please enter the desired phi (in degrees) at which the critical angle phi should be calculated:'); %degrees , this will be the phi used to calculate the critical angle theta. This phi can be set to zero to compare the results with the answers from the statics textbook.

disp('Please allow 30 seconds for the code to run')%This is to notify the user that the program takes a while to display the results. Compared to the conventional ways of finding such results, the program is still quite efficient and convinient to use despite its long run time.

%% Calculations

criticalThetas=zeros(1,5);%Pre-allocates memory space for critical thetas to save time for efficiency of code. 
criticalThetas=thetaCritical(W1,W2,W3,L1,L2,L3,L4,L5,desiredPhi); %degrees, calls the function thetaCritical to calculate critical thetas. 
criticalPhisD=zeros(1,5);%Pre-allocates memory space for critical phis to save time for efficiency of code. 
criticalPhisB=zeros(1,5);%Pre-allocates memory space for critical phis to save time for efficiency of code. 
[criticalPhisD,criticalPhisB]=criticalPhis(W1,W2,W3,L1,L2,L3,L4,L5,L7,theta); %degrees, calls the function criticalPhis to calculate both critical angles phiB and phiD which are on opposite sides of each other.

%The spacing operations below create a vector of seven elements from zero
%to each critical theta and critical phi corresponding to each combination
%of L1 and L2 to be used in the calculation of the reactions at outriggers
%A,B,and D. 
thetas1=linspace(0,criticalThetas(1),7);%degrees
thetas2=linspace(0,criticalThetas(2),7);%degrees
thetas3=linspace(0,criticalThetas(3),7);%degrees
thetas4=linspace(0,criticalThetas(4),7);%degrees
thetas5=linspace(0,criticalThetas(5),7);%degrees

phisB1=linspace(0,criticalPhisB(1),7);%degrees
phisB2=linspace(0,criticalPhisB(2),7);%degrees
phisB3=linspace(0,criticalPhisB(3),7);%degrees
phisB4=linspace(0,criticalPhisB(4),7);%degrees
phisB5=linspace(0,criticalPhisB(5),7);%degrees

phisD1=linspace(0,criticalPhisD(1),7);%degrees
phisD2=linspace(0,criticalPhisD(1),7);%degrees
phisD3=linspace(0,criticalPhisD(1),7);%degrees
phisD4=linspace(0,criticalPhisD(1),7);%degrees
phisD5=linspace(0,criticalPhisD(1),7);%degrees

%Each line below calls for a specific function that calculates the
%reactions at outriggers A, B, and D at each combination of L1 and L2 
[NA1,NA2,NA3,NA4,NA5]=NAreaction(W1,W2,W3,L1,L2,L3,L4,L5,L7,desiredPhi,thetas1,thetas2,thetas3,thetas4,thetas5);%N
[NB1,NB2,NB3,NB4,NB5] = NBreaction(W1,W2,W3,L1,L2,L3,L4,L5,L7,phisB1,phisB2,phisB3,phisB4,phisB5,theta);%N
[ND1,ND2,ND3,ND4,ND5] = NDreaction(W1,W2,W3,L1,L2,L3,L4,L5,L7,phisD1,phisD2,phisD3,phisD4,phisD5,theta);%N


%% Outputs

%The fprintf command below prints the three critical angles calculated for
%each combination of L1 and L2 enetered by the user. 
fprintf('The critical angles for each combination of the lengths L1 and L2 is as follows: .\n')
for i=1:length(L1)
fprintf('For L1=%g meters and L2=%g meters, the critical angle theta equals %g degrees, the critical angle phi on the left side equals %g degrees, and the critical angle phi on the right side equals %g degrees \n',L1(i),L2(i),criticalThetas(i),criticalPhisD(i),criticalPhisB(i))
end 

[minTheta,idx]=min(criticalThetas); %finds the minimum value and index of the angle. The index will be used to link the angle to the combination of L1 and L2 that resulted in it.
[minPhi,indx]=min(criticalPhisD); %finds the minimum value and index of the angle. The index will be used to link the angle to the combination of L1 and L2 that resulted in it.
[minPhiB,indxB]=min(criticalPhisB); %finds the minimum value and index of the angle. The index will be used to link the angle to the combination of L1 and L2 that resulted in it.
[maxTheta,index]=max(criticalThetas); %finds the maximum value and index of the angle. The index will be used to link the angle to the combination of L1 and L2 that resulted in it.
[maxPhiB,idexB]=max(criticalPhisB); %finds the maximum value and index of the angle. The index will be used to link the angle to the combination of L1 and L2 that resulted in it.
[maxPhi,idex]=max(criticalPhisD); %finds the maximum value and index of the angle. The index will be used to link the angle to the combination of L1 and L2 that resulted in it.

fprintf('The minimum critical angle theta occurs at lengths L1=%g meters and L2=%g meters.\n',L1(idx),L2(idx))%Uses the index found above to display the minimum critical angle along with its associated L1 and L2.
fprintf('The minimum critical angle phi on the left side occurs at lengths L1=%g meters and L2=%g meters.\n',L1(indx),L2(indx))%Uses the index found above to display the minimum critical angle along with its associated L1 and L2.
fprintf('The minimum critical angle phi on the right side occurs at lengths L1=%g meters and L2=%g meters.\n',L1(indxB),L2(indxB))%Uses the index found above to display the minimum critical angle along with its associated L1 and L2.
fprintf('The maximum critical angle theta occurs at lengths L1=%g meters and L2=%g meters.\n',L1(index),L2(index))%Uses the index found above to display the maximum critical angle along with its associated L1 and L2.
fprintf('The maximum critical angle phi on the left side occurs at lengths L1=%g meters and L2=%g meters.\n',L1(idex),L2(idex))%Uses the index found above to display the maximum critical angle along with its associated L1 and L2.
fprintf('The maximum critical angle phi on the right side occurs at lengths L1=%g meters and L2=%g meters.\n',L1(idexB),L2(idexB))%Uses the index found above to display the maximum critical angle along with its associated L1 and L2.

%% Plots 
%Plots are an important part of the output, and thus, have a dedicated section. 

%Plot #1 , this plot plots the critical angles phiB which are on the right
%against the reaction at outrigger B, which is on the same side. This helps
%the user determine whether or not the outrigger B will be able to
%withstand the force inflicted on it as a result of the range of angles
%phiB from zero to the critical angle for each combination of L1 and L2. 
subplot(3,1,1)
plot(phisB1,NB1)
hold on
plot(phisB2,NB2)
hold on
plot(phisB3,NB3)
hold on
plot(phisB4,NB4)
hold on
plot(phisB5,NB5)
hold off
axis tight
xlabel('Angle \phi_B (degrees)')
ylabel('Reaction at outrigger B (Newtons)')
legend('At 1st length combination','At 2nd length combination','At 3rd length combination','At 4th length combination','At 5th length combination','Location','bestoutside')%Because the location of the data on the plot will vary according to the data entered by the user, the best place for the legend is outside the axes of the plot.
title('Angle \phi_B VS. Reaction at outrigger B')

%Plot #2 , this plot plots the critical angles phiB which are on the left
%against the reaction at outrigger D, which is on the same side. This helps
%the user determine whether or not the outrigger D will be able to
%withstand the force inflicted on it as a result of the range of angles
%phiD from zero to the critical angle for each combination of L1 and L2. 

subplot(3,1,2)
plot(phisD1,ND1)
hold on
plot(phisD2,ND2)
hold on
plot(phisD3,ND3)
hold on
plot(phisD4,ND4)
hold on
plot(phisD5,ND5)
hold off
axis tight
xlabel('Angle \phi_D (degrees)')
ylabel('Reaction at outrigger D (Newtons)')
legend('At 1st length combination','At 2nd length combination','At 3rd length combination','At 4th length combination','At 5th length combination','Location','bestoutside')%Because the location of the data on the plot will vary according to the data entered by the user, the best place for the legend is outside the axes of the plot.
title('Angle \phi_D VS. Reaction at outrigger D')

%Plot #3, this plot plots the reaction at outrigger A against the angle
%theta. This plot helps the user determine the value of the suitable angle
%theta at each combination of L1 and L2 before the crane tips when the
%reaction at outrigger A equals zero. 

subplot(3,1,3)
plot(thetas1,NA1)
hold on
plot(thetas2,NA2)
hold on
plot(thetas3,NA3)
hold on
plot(thetas4,NA4)
hold on
plot(thetas5,NA5)
hold off
axis tight
xlabel('Angle \Theta (degrees)')
ylabel('Reaction at outrigger A (Newtons)')
legend('At 1st length combination','At 2nd length combination','At 3rd length combination','At 4th length combination','At 5th length combination','Location','bestoutside')%Because the location of the data on the plot will vary according to the data entered by the user, the best place for the legend is outside the axes of the plot.
title('Angle \Theta VS. Reaction at outrigger A')

%% Functions 
%Functions are an important part of the code as they make the main code
%more concise and faster to run. 

%Function #1: thetaCritical 

%This function calculates the critical theta at a pre-defined angle phi (entered by the user) for each L1 and L2 combination
%by looping through their vectors and applying each value pair to the
%equation inside the loop.This function uses the built-in function (vpasolve) to
%numerically solve the equilibrium equation and find an answer (critical
%theta) at each condition. The resulting angles are in degrees.

%Notes: vpasolve is a built-in function that was not discussed in class. 
function[criticalThetas]=thetaCritical(W1,W2,W3,L1,L2,L3,L4,L5,desiredPhi)
for i=1:length(L1)
    syms criticalTheta
    criticalThetas(i)=vpasolve(-W1*(L3+L4)-W2*(L3+(L2(i)*sind(criticalTheta)*cosd(desiredPhi)))-W3*(L3+((L1(i)+L2(i))*sind(criticalTheta)*cosd(desiredPhi)))+((W1+W2+W3)*(L3+L4+L5))==0,criticalTheta); %degrees
end
end

%Function #2: NA

%This function calculates the reaction at outrigger A at each L1 and L2
%combination for each vector of thetas ranging from zero to the critical
%theta at each L1 and L2 combination. This function uses the built-in
%function vpasolve to find the resutls. The resulting reactions are in
%Newtons. 

function [NA1,NA2,NA3,NA4,NA5]=NAreaction(W1,W2,W3,L1,L2,L3,L4,L5,L7,desiredPhi,thetas1,thetas2,thetas3,thetas4,thetas5)
for i=1:length(thetas1)
   syms NA
   NA1(i)=vpasolve(-W1*(L3+L4)-W2*(L3+(L2(1)*sind(thetas1(i))*cosd(desiredPhi)))-W3*(L3+((L1(1)+L2(1))*sind(thetas1(i))*cosd(desiredPhi)))+((W1+W2+W3-NA)*(L3+L4+L5))==0,NA); %N
end

for i=1:length(thetas2)
   syms NA
   NA2(i)=vpasolve(-W1*(L3+L4)-W2*(L3+(L2(2)*sind(thetas2(i))*cosd(desiredPhi)))-W3*(L3+((L1(2)+L2(2))*sind(thetas2(i))*cosd(desiredPhi)))+((W1+W2+W3-NA)*(L3+L4+L5))==0,NA); %N
end

for i=1:length(thetas3)
   syms NA
   NA3(i)=vpasolve(-W1*(L3+L4)-W2*(L3+(L2(3)*sind(thetas3(i))*cosd(desiredPhi)))-W3*(L3+((L1(3)+L2(3))*sind(thetas3(i))*cosd(desiredPhi)))+((W1+W2+W3-NA)*(L3+L4+L5))==0,NA); %N
end

for i=1:length(thetas4)
   syms NA
   NA4(i)=vpasolve(-W1*(L3+L4)-W2*(L3+(L2(4)*sind(thetas4(i))*cosd(desiredPhi)))-W3*(L3+((L1(4)+L2(4))*sind(thetas4(i))*cosd(desiredPhi)))+((W1+W2+W3-NA)*(L3+L4+L5))==0,NA); %N
end

for i=1:length(thetas5)
   syms NA
   NA5(i)=vpasolve(-W1*(L3+L4)-W2*(L3+(L2(5)*sind(thetas5(i))*cosd(desiredPhi)))-W3*(L3+((L1(5)+L2(5))*sind(thetas5(i))*cosd(desiredPhi)))+((W1+W2+W3-NA)*(L3+L4+L5))==0,NA); %N
end
end

%Function #3: criticalPhis

%This function calculates the critical phis at both sides of the y-axis at a pre-defined angle theta (entered by the user) for each L1 and L2 combination
%by looping through their vectors and applying each value pair to the
%equation inside the loop.This function uses the built-in function (vpasolve) to
%numerically solve the equilibrium equation and find an answer (critical
%phi) at each condition. The resulting angles are in degrees.

function [phisD,phisB] = criticalPhis(W1,W2,W3,L1,L2,L3,L4,L5,L7,theta)
for i=1:length(L1)
syms phiD
phisD(i)=vpasolve(((-W1*(L3+L4))-(W2*(L3+(L2(i)*sind(theta)*cosd(phiD))))-(W3*(L3+(L1(i)+L2(i))*sind(theta)*cosd(phiD)))+(((-sind(phiD)*(L3+L4+L5)*(W2*L2(i)*sind(theta)+W3*L1(i)*sind(theta)+W3*L2(i)*sind(theta)))/L7)*(L3+L4+L5)))==0,phiD); %degrees
end 

for i=1:length(L1)
    syms phiB
    phisB(i)=vpasolve(((-W1*(L3+L4))-(W2*(L3+(L2(i)*sind(theta)*cosd(phiB))))-(W3*(L3+(L1(i)+L2(i))*sind(theta)*cosd(phiB)))+(((sind(phiB)*(L3+L4+L5)*(W2*L2(i)*sind(theta)+W3*L1(i)*sind(theta)+W3*L2(i)*sind(theta)))/L7)*(L3+L4+L5)))==0,phiB); %degrees
end
end

%Function #4: NB

%This function calculates the reaction at outrigger B at each L1 and L2
%combination for each vector of thetas ranging from zero to the critical
%theta at each L1 and L2 combination. This function uses the built-in
%function vpasolve to find the resutls. The resulting reactions are in
%Newtons. 

function [NB1,NB2,NB3,NB4,NB5] = NBreaction(W1,W2,W3,L1,L2,L3,L4,L5,L7,phisB1,phisB2,phisB3,phisB4,phisB5,theta)
for i=1:length(phisB1)
    syms NB
    NB1(i)=vpasolve((W2*L2(1)*sind(theta)*sind(phisB1(i)))+(W3*(L1(1)+L2(1))*sind(theta)*sind(phisB1(i)))-(NB*L7)+(L7*((((W1*(L3+L4)+(W2*(L3+L2(1)*sind(theta)*cosd(phisB1(i))))+(W3*(L3+(L1(1)+L2(1))*sind(theta)*cosd(phisB1(i))))))/(L3+L4+L5))-NB))==0,NB);%N
end
for i=1:length(phisB2)
    syms NB
    NB2(i)=vpasolve((W2*L2(2)*sind(theta)*sind(phisB2(i)))+(W3*(L1(2)+L2(2))*sind(theta)*sind(phisB2(i)))-(NB*L7)+(L7*((((W1*(L3+L4)+(W2*(L3+L2(2)*sind(theta)*cosd(phisB2(i))))+(W3*(L3+(L1(2)+L2(2))*sind(theta)*cosd(phisB2(i))))))/(L3+L4+L5))-NB))==0,NB);%N
end
for i=1:length(phisB3)
    syms NB
    NB3(i)=vpasolve((W2*L2(3)*sind(theta)*sind(phisB3(i)))+(W3*(L1(3)+L2(3))*sind(theta)*sind(phisB3(i)))-(NB*L7)+(L7*((((W1*(L3+L4)+(W2*(L3+L2(3)*sind(theta)*cosd(phisB3(i))))+(W3*(L3+(L1(3)+L2(3))*sind(theta)*cosd(phisB3(i))))))/(L3+L4+L5))-NB))==0,NB);%N
end
for i=1:length(phisB4)
    syms NB
    NB4(i)=vpasolve((W2*L2(4)*sind(theta)*sind(phisB4(i)))+(W3*(L1(4)+L2(4))*sind(theta)*sind(phisB4(i)))-(NB*L7)+(L7*((((W1*(L3+L4)+(W2*(L3+L2(4)*sind(theta)*cosd(phisB4(i))))+(W3*(L3+(L1(4)+L2(4))*sind(theta)*cosd(phisB4(i))))))/(L3+L4+L5))-NB))==0,NB);%N
end
for i=1:length(phisB5)
    syms NB
    NB5(i)=vpasolve((W2*L2(5)*sind(theta)*sind(phisB5(i)))+(W3*(L1(5)+L2(5))*sind(theta)*sind(phisB5(i)))-(NB*L7)+(L7*((((W1*(L3+L4)+(W2*(L3+L2(5)*sind(theta)*cosd(phisB5(i))))+(W3*(L3+(L1(5)+L2(5))*sind(theta)*cosd(phisB5(i))))))/(L3+L4+L5))-NB))==0,NB);%N
end
end

%Function #5: ND

%This function calculates the reaction at outrigger D at each L1 and L2
%combination for each vector of thetas ranging from zero to the critical
%theta at each L1 and L2 combination. This function uses the built-in
%function vpasolve to find the resutls. The resulting reactions are in
%Newtons. 

function [ND1,ND2,ND3,ND4,ND5] = NDreaction(W1,W2,W3,L1,L2,L3,L4,L5,L7,phisD1,phisD2,phisD3,phisD4,phisD5,theta)
for i=1:length(phisD1)
    syms ND
    ND1(i)=vpasolve((W2*L2(1)*sind(theta)*sind(phisD1(i)))+(W3*(L1(1)+L2(1))*sind(theta)*sind(phisD1(i)))+(ND*L7)-(L7*((((W1*(L3+L4)+(W2*(L3+L2(1)*sind(theta)*cosd(phisD1(i))))+(W3*(L3+(L1(1)+L2(1))*sind(theta)*cosd(phisD1(i))))))/(L3+L4+L5))-ND))==0,ND);%N
end
for i=1:length(phisD2)
    syms ND
    ND2(i)=vpasolve((W2*L2(2)*sind(theta)*sind(phisD2(i)))+(W3*(L1(2)+L2(2))*sind(theta)*sind(phisD2(i)))+(ND*L7)-(L7*((((W1*(L3+L4)+(W2*(L3+L2(2)*sind(theta)*cosd(phisD2(i))))+(W3*(L3+(L1(2)+L2(2))*sind(theta)*cosd(phisD2(i))))))/(L3+L4+L5))-ND))==0,ND);%N
end
for i=1:length(phisD3)
    syms ND
    ND3(i)=vpasolve((W2*L2(3)*sind(theta)*sind(phisD3(i)))+(W3*(L1(3)+L2(3))*sind(theta)*sind(phisD3(i)))+(ND*L7)-(L7*((((W1*(L3+L4)+(W2*(L3+L2(3)*sind(theta)*cosd(phisD3(i))))+(W3*(L3+(L1(3)+L2(3))*sind(theta)*cosd(phisD3(i))))))/(L3+L4+L5))-ND))==0,ND);%N
end
for i=1:length(phisD4)
    syms ND
    ND4(i)=vpasolve((W2*L2(4)*sind(theta)*sind(phisD4(i)))+(W3*(L1(4)+L2(4))*sind(theta)*sind(phisD4(i)))+(ND*L7)-(L7*((((W1*(L3+L4)+(W2*(L3+L2(4)*sind(theta)*cosd(phisD4(i))))+(W3*(L3+(L1(4)+L2(4))*sind(theta)*cosd(phisD4(i))))))/(L3+L4+L5))-ND))==0,ND);%N
end
for i=1:length(phisD5)
    syms ND
    ND5(i)=vpasolve((W2*L2(5)*sind(theta)*sind(phisD5(i)))+(W3*(L1(5)+L2(5))*sind(theta)*sind(phisD5(i)))+(ND*L7)-(L7*((((W1*(L3+L4)+(W2*(L3+L2(5)*sind(theta)*cosd(phisD5(i))))+(W3*(L3+(L1(5)+L2(5))*sind(theta)*cosd(phisD5(i))))))/(L3+L4+L5))-ND))==0,ND);%N
end
end






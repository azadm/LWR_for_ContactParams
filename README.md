# LWR_for_ContactParams
Locally Weighted Regression (LWR) algorithm for learning contact parameters of compliant contacts.
Running function [time, P, V, F] = data_2_state_force in MATLAB extracts the position, velocity and force of the contact points on the left or right foot of the robot.   
Running function main.m for different mat files (lines 12 to 14) compares predicted values of force (from LWR algorithm) with the actual force from sensors for the recorded data.

function [time, P, V, F] = data_2_state_force( m )
% data_2_state_force calculates position, velocity and forces at the
% contact points of the right foot.
% This function extracts the position and velocity of the contact points
% from the 53 DoF readings.  time is a vector of time (n by 1), where n is
% the number of samples.  Assuming that the robot's right foot has m
% contact points, P is a 3 times m by n matrix including the positions of
% the contact points and V is a 3 times m by n matrix containing the
% velocities of those points.  F is a 6 by n vector including the total
% force and torque acting on the robot's right foot.  This is estimated
% from F/T sensors. Here, it is assumed that the base link is the left
% foot.

if nargin < 1
    % defining 4 contact points on the foot
    m = 5;
    r(:,1) = [0.05;0.02;0];
    r(:,2) = [-0.05;0.02;0];
    r(:,3) = [0.05;-0.02;0];
    r(:,4) = [-0.05;-0.02;0];
end

% load recorded joint data
load('Joint_data');
% extract the required data from the whole data (32 DoF to use in IK)
[time, matrix32, labels32] = convert_data_from_53_to_32 ( matrix, labels );
% calculate the robot's states
[q, qd, qdd] = extract_states ( time, matrix32 );


% set up dynComp class for dynamics calculations
dynComp = iDynTree.DynamicsComputations();
dynComp.loadRobotModelFromFile('icub.urdf');
dof = dynComp.getNrOfDegreesOfFreedom();

gravity = iDynTree.SpatialAcc();
gravity.zero();
gravity.setVal(2, -9.81);

% % define the robot's state vectors
qj = iDynTree.VectorDynSize(dof);
qdj = iDynTree.VectorDynSize(dof);
qddj = iDynTree.VectorDynSize(dof);

% get frame id s
l_sole_id = dynComp.getFrameIndex('l_sole');
r_sole_id = dynComp.getFrameIndex('r_sole');
l_foot_id = dynComp.getFrameIndex('l_foot');
r_foot_id = dynComp.getFrameIndex('r_foot');

% set foating base link
dynComp.setFloatingBase('l_foot');
fl_base = dynComp.getFloatingBase();


% time_index = 10;
% for i = time_index
for i = 1:size(q,2)
    Q = q(:,i);
    Qd = qd(:,i);
    Qdd = qdd(:,i);
    % set the values of q and qdot
    qj.fromMatlab(Q);
    qdj.fromMatlab(Qd);
    qddj.fromMatlab(Qdd);
    
    
    % set robot states
    dynComp.setRobotState(qj,qdj,qddj,gravity);
    
    % left foot position and rotation matrix
    R = dynComp.getRelativeTransform(l_foot_id,l_sole_id);
    p_lfoot = R.getPosition().toMatlab();
    R_lfoot = R.getRotation().toMatlab();
    
    % right foot position and rotation matrix
    R = dynComp.getRelativeTransform(l_foot_id,r_sole_id);
    p_rfoot = R.getPosition().toMatlab();
    R_rfoot = R.getRotation().toMatlab();

    
    % left foot twist
    w = dynComp.getFrameTwist(l_sole_id);
    w_lfoot = w.toMatlab();
  
    % right foot twist
    w = dynComp.getFrameTwist(r_sole_id);
    w_rfoot = w.toMatlab();
    
    % Jacobian calculations
    lfoot_jac = iDynTree.MatrixDynSize(6,6+dof);
    rfoot_jac = iDynTree.MatrixDynSize(6,6+dof);
    jac_calc1 = dynComp.getFrameJacobian(l_sole_id, lfoot_jac);
    jac_calc2 = dynComp.getFrameJacobian(r_sole_id, rfoot_jac);
    
    % convert Jacobians from idyntree to matlab
    J_lfoot = lfoot_jac.toMatlab();
    J_rfoot = rfoot_jac.toMatlab();
    
    % position and velocity matrices
    P(1:3,i) = p_rfoot;
    vel = J_rfoot*[zeros(6,1); Qd];
    V(1:3,i) = vel(1:3);
    omega = skew(vel(4:6));
    for j = 1:(m-1)
        P(3*j+1:3*j+3,i) = p_rfoot + R_rfoot*r(:,j);
        V(3*j+1:3*j+3,i) = vel(1:3) + omega*R_rfoot*r(:,j);
    end
    % twist matrix
    W(:,i) = w_rfoot;
end


% load recorded force/torque data
load('FT_data');

time_force = matrix(:,1);
F_l = matrix(:,14:16)';
T_l = matrix(:,17:19)';
F_r = matrix(:,20:22)';
T_r = matrix(:,23:25)';

F = [F_r; T_r];

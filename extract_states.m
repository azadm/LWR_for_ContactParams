function [q, qd, qdd] = extract_states ( time, matrix )

% q and qdot
for i = 1:size(matrix,1)
    [SG0, SG1 ] = sgolayDiff(matrix(i,:), time);
    q(i,:) = SG0;
    qd(i,:) = SG1;
end

% qddot
for i = 1:size(matrix,1)
    [SG0, SG1 ] = sgolayDiff(qd(i,:), time);
    qdd(i,:) = SG1;
end
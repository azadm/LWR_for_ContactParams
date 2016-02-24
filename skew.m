function Skew_Mat = skew( Vector )

Skew_Mat = [0, -Vector(3), Vector(2);
    Vector(3), 0, -Vector(1);
    -Vector(2), Vector(1), 0];


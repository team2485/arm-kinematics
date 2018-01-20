syms t m1 m2 g d1 d2 J1 J2 L1 T1 T2 Fx Fy theta1 theta2 w1 w2 a1 a2 real;
theta1Function(t) = theta1 + w1*t + a1*t^2/2;
theta2Function(t) = theta2 + w2*t + a2*t^2/2;
x1 = d1 * [cos(theta1Function) sin(theta1Function) 0];
x2 = d2 * [cos(theta2Function) sin(theta2Function) 0];
R = L1 * [cos(theta1Function) sin(theta1Function) 0];
W1 = m1 * g * [0 -1 0];
W2 = m2 * g * [0 -1 0];
F = [Fx Fy 0];
zHat = [0 0 1];
Rdotdot = diff(R, 2);
x2dotdot = diff(x2, 2);


eq1 = J2 * a2 == dot(cross(x2(0), W2), zHat) + T2;
eq2 = J1 * a1 == dot(cross(x1(0), W1) - cross(R(0), F), zHat) + T1 - T2;
eq3 = m2 * (Rdotdot(0) + x2dotdot(0)) == F + W2;
[T1sln T2sln Fxsln Fysln] = solve([eq1 eq2 eq3], [T1 T2 Fx Fy]);
simplify(T1sln)
simplify(T2sln)
simplify (Fxsln)
simplify(Fysln)





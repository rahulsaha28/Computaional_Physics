#     Assignment 04 : Solving eigenvalue problem using Sturm-Liouville problems
#     Name: Rahul Saha
#     Reg. No: 2014132028
#     Session : 2014-2015
#     python version: 3.6

import matplotlib.pyplot as plt;
import math;


#  here I use
# c2 = (2p_i + h(p_i)' )
# c1 = 4p_i + 2h^2 q_i
# c0 = ( h(p_i)-2p_i )
# d = 2h^2*s_i

#          c1 ∗ u_i − c0 ∗ u_i−1 + d
# u_i+1 = ----------------------------
#                    c2


# required fun p(x)
def p(x):
    return 1 - x * x;


# required fun p'(x)
def fp(x):
    return -2 * x;


# required fun q(x)
def q(x, l):
    return l * (l + 1);


# required fun s(x)
def s(x):
    return 0;


# this function find the eigenfunction of corresponding eigenvalue l
def stufun(l):
    x = 0;
    h = 0.01;
    pos = [x, h];
    # given boundary condition is u(0) = 0, u(1) = 1
    ui = [0, h];  # (0 is the ui-1 value ) , (h  is the ui value)
    i = 0;
    while x < 1:
        #         all coefficient are here
        c2 = h * fp(x) + 2 * p(x);
        c1 = 4 * p(x) - 2 * h * h * q(x, l);
        c0 = 2 * p(x) - h * fp(x);
        d = 2 * h * h * s(x);

        w = (c1 * ui[i + 1] - (c0 * ui[i]) + d) / c2
        ui.append(w);
        x = x + h;
        pos.append(x)  # update the position
        i = i + 1;
    return [pos, ui];


# this function find the root using secand method of the function f(l) =  Ul(1) - U(1)
def secant(l):
    dx = 0.01;
    error = 1 * math.exp(-6);
    l1 = l + dx;
    while abs(dx) > error:
        #       here we represent f(x)
        time, y1 = stufun(l);
        #     hrer we represent f(x+h)
        time, y2 = stufun(l1);

        fx = y1[-1] - 1;
        fx1 = y2[-1] - 1;
        #       secand method
        d = fx1 - fx;
        l2 = l1 - fx1 * (l1 - l) / d;
        l = l1;
        l1 = l2;
        dx = l1 - l;
    return l1;


# now assign the eigen value of l=1;
l = 1;

a = secant(l);

result = stufun(a);
plt.plot(result[0], result[1]);
plt.xlabel("value of x");
plt.ylabel("U(x)");
plt.legend(['Eigenvalue is ' + str(round(a))]);
plt.show();




#     Assignment 03 : Solving seconder order differential equation using RK method 12 order
#     Name: Rahul Saha
#     Reg. No: 2014132028
#     Session : 2014-2015
#     python version: 3.6


# Runge–Kutta methods is  implicit or explicit iterative method use for find the solution curve of a differential EQ
# if I  use the method at first I convert the given differential EQ in the form like below
#
#    y' = f(x,y)  [ y(x0) = y0 which is given initial condition ]
#
#                                  s
# Then the rule is  y_t+h = y_t + h*Σa_i*k_i  ----------------------(1)
#                                  i
# where
#                       s
#       k_i = f(y_t + h*ΣB_ij*k_j, t_n + α_i*h)         where s represent which order RK method
#                      j=1
#
# By use the constant α_i and B_ij as the value (from https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
#
#     α_i            B_ij
#     α_1 = 0        B_21 = 1/2
#     α_2 = 1/2      B_32 = 1/2
#     α_3 = 1/2      B_43 = ,,         --------------> for n order RK method
#     α_4 = ,,       .........
#    .........       .........
#    .........       B_n(n-1) = 0
#    a_n = 1
#
#  now expend the both site of EQ (1) using Taylor's series we get
#  LHS is
#                              h^2   d             h^3   d^2
#  y_t+h = y_t + h*f(y_t, t) + ---- --f(y_t, t) + ----  -----f(y_t, t) + ........ ----------->(2)
#                              2!  dx             3!    dx^2
#
# RHS is
# I can not represent this here
# now compairing the two side of each coefficient I get
#
#     a + b + c + d = 1
#   (1/2)*b +(1/2)c +d = 1/2
#   (1/4)c + (1/2)d = 1/6
#            (1/4)d = 1/24
#
#
import matplotlib.pyplot as plt
import math;


# this function take an uppertrainglur matrix and solve the linear system of equation and find x,y,z....
def unknownSolve(a):
    result = [];  # where z , y , x are append

    x = 0;
    for i in range(len(a) - 1, -1, -1):  # accessing last row
        index = 0;
        for j in range(len(a[0]) - 1, -1, -1):  # accessing last colum

            # if last element of then  colun add to x
            if j == len(a[0]) - 1:
                x = -(a[i][j]);
            # if find a element like a33 where i == j divide by that element and add result list;
            elif i == j:
                x = x / a[i][j];
                result.append(x);
                break;  # I don't run this for loop any more
            else:
                x = x - (a[i][j] * result[index]);
                index = index + 1;
    result.reverse();
    return result;


# this function is use to find constant a, b, c, d, e, f ..... for a given value of n where n is order
def f(n):
    a = [];
    c = 1;
    k = 0;
    for i in range(n):  # row accessing loop
        a.append([]);
        for j in range(n + 1):  # column accessing loop
            # in this last column (where j=n ) I wont to put 1/2! , 1/3!, 1/4! ......
            if j == n:
                a[i].append(-1 / c);
            #             in this row (where i=0) I want to put 1, 1, 1, 1, 1, 1........
            elif i == 0:
                a[i].append(1);

            elif i == 1 and j < 3 and j > 0:
                a[i].append(1 / 2);

            #             in this row I set 0 when j<i
            elif j < i:
                a[i].append(0);

            #             But what about the other condition like j==i  or j>i
            elif j >= i:
                if j == i:
                    a[i].append(1 / 4);
                    k = j + 1
                elif k == j:
                    a[i].append(1 / 2);
                else:
                    a[i].append(1);
        c = c * (i + 2);

    return unknownSolve(a);


# this function set orbitrary αi and return the value
def findAlp(n):
    alp = [];

    for i in range(n):
        #         set α1 as the value 0
        if i == 0:
            alp.append(0);
        #         set αn as the value 1
        elif i == n - 1:
            alp.append(1);
        #         otherwise set the 1/2
        else:
            alp.append(1 / 2);
    return alp;


# this function set orbitrary βij and return the value
def findBeta(n):
    beta = [];
    for i in range(n):
        beta.append([]);
        for j in range(n):
            #             set the value at β[n][n-1] = 1
            if i == n - 1 and j == n - 2:
                beta[i].append(1);
            #             set the value at β[i][i-1]= 1/2 so j==(i-1)
            elif j == i - 1:
                beta[i].append(1 / 2);
            #             set β[][] = 0 otherwise
            else:
                beta[i].append(0);
    return beta;


# this is the main function which take fig of the function, initial value of x0 and y0 and h and order n return k1, k2, k3 ....
#  return coefficient
def findK(fig, y1, y2, h, t, n):
    alpha = findAlp(n);
    beta = findBeta(n);

    k = [];
    kk = [0];
    for i in range(n):  # accessing the row of Beta[i][j]
        cons = 0;
        for j in range(i):  # accessing the col of Beta[i][j]
            cons = cons + beta[i][j] * kk[j];  # ---------------->B[i][j]*k[j]

        k.append(fig(y1 + h * cons, y2 + h * cons, t + h * alpha[i]));  # ------------->k_i
        kk = k;

    #     print(k);
    a = f(n);
    con = 0;
    for j in range(n):
        con = con + a[j] * k[j];  # ------------------------> Σa[j]*k[i]

    return h * con;  # ---------------->h*Σa[j]*k[i]


# -----------------------------------solving second order diffrential equation using Rk12--------------------
#
# θ'' + qθ + sinθ = bcos(wt)
#
#  where q = 1.15, b = 0.9, w = 2/3
#  initial value------------>θ(0) = 0 so θ = 0 when t = 0
#                            θ'(0) = 1 so θ' = 1 when t = 0


# this is the function that return the θ = x
def positionfun(y1, y2, t0):
    return y2;


# this is the second function that return θ' = v
def velocityfun(y1, y2, t0):
    return -1.15 * y2 - math.sin(y1) + 0.9 * math.cos((2 / 3) * t0);


# creating graph
# where y1 is the initial θ,  y2 is the initial θ', h is the interval, t0 is the initial time, n is the oreder of RK method
def createCurve(y1, y2, h, t0, n):
    y_first = [y1];
    y_second = [y2];
    time = [t0];
    t = t0;
    i = 0;

    while t < 50:
        y_first.append(y_first[i] + findK(positionfun, y_first[i], y_second[i], h, t, n))
        y_second.append(y_second[i] + findK(velocityfun, y_first[i], y_second[i], h, t, n))
        t = t + h;
        time.append(t);
        i = i + 1;
    #this is the place  where plt is use to visulize the graph
    plt.subplot(121)
    plt.plot(time, y_first, 'ro');
    plt.plot(time, y_second, 'bo');
    plt.xlabel('Time')
    plt.ylabel('Magnitude')
    plt.legend(['position', 'velocity'])

    plt.subplot(122)
    plt.plot(y_first, y_second);
    plt.xlabel('position')
    plt.ylabel('velocity')
    plt.legend(['velocity vs position'])
    plt.savefig('r2.pdf')

    plt.show();

createCurve(0, 1, 0.1, 0, 12);

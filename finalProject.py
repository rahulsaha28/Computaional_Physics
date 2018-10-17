
'''
  using   python  v3.4
'''

# importing important lib...
import matplotlib.pyplot as plt;
import math;
from mpl_toolkits.mplot3d import Axes3D


# The  differential equation which control the motion of the pendulum are given below
'''
      g                        cos(ϑ)
      -- * sin(ϑ) +  ϑ'' -   ---------- * h^2 = 0
      l                       (sin(ϑ))^3
                                         
  which get from the L = T - V 
         
     d   d(L)      d(L)
    ---( ---- ) - ------ = 0
    dt    dϑ'       dϑ      
    
    where 
         
         constant = h = m * l^2 * (sin(ϑ))^2 * ϕ'
         
     
                     cos(ϑ)               g
     ϑ'' = y2' =  -------------* h^2  -  ---- * sin(ϑ)  
                   (sin(ϑ))^3             l
                   
                   
     for this process using Rangu kutta method -------->
     ϑ' = y2
     ϑ = y1              
                                  
'''

'''
for adding drag force in the ϑ direction the fig of the equation look like

                    cos(ϑ)                g
     ϑ'' = y2' =  -------------* h^2  -  ---- * sin(ϑ) - b* ϑ'  where b is the drag force
                   (sin(ϑ))^3             l
                                      

'''

'''
for varying mass the equation is look like

h = (m - r*t) * l^2 * (sin(ϑ))^2 * ϕ'   where r is the rate of change (for mass)

here h is not constant in time. But constant in moment

'''

'''
for canonical pendulum  ϑ = constant
so the value of  ϕ' is given by
                         g
            ϕ' = sqrt(---------)
                       l*cos(ϑ)
                       
'''

# now the constant are here
# g = 9.81;
# l = 1;
# r = 0.6;
# b = 0;
# h = .5;
# m = .1;



def pendulum_graph(g, l, m, b, control1, control2):

    h = .4;
    # this fun is for find y1 = ϑ

    def g1(y1, y2, t):
        return y2;

    # this fun is for find y2 = ϑ'

    def g2(y1, y2, t):
        return (h * h * math.cos(y1) / (math.sin(y1)) ** 3) - (g / l) * math.sin(y1) - b * y2;

    # this is the fun for z1 = ϕ
    def g3(y1, y2, t):
        return h / (math.sin(y1)) ** 2;

    # find coefficient
    def coeffi(fig, y1, y2, t, tau):
        k1 = tau * fig(y1, y2, t);
        k2 = tau * fig(y1 + k1 / 2, y2 + k1 / 2, t + tau / 2);
        k3 = tau * fig(y1 + k2 / 2, y2 + k2 / 2, t + tau / 2);
        k4 = tau * fig(y1 + k3, y2 + k3, t + tau);

        return (k1 + 2 * k2 + 2 * k3 + k4) / 6;

    def R_kutta():
        t = 0;
        tau = 0.01;
        time = [t];
        y1 = [math.pi / 6];
        y2 = [math.pi/6];
        z1 = [math.pi/6];
        z2 = [];
        m_y1 = [];
        m_z1 = [];
        k_y1 = [];
        k_z1 = [];
        total_k = [];
        U = [];
        Total_E = [];

        z2.append( h/(math.sin(y1[-1])* math.sin(y1[-1])) );
        m_y1.append(m*l*l*y2[-1]);
        m_z1.append(m*l*l*math.sin(y1[-1])*math.sin(y1[-1])*z2[-1]);
        k_y1.append((m_y1[-1]*m_y1[-1])/(2*m));
        k_z1.append((m_z1[-1] * m_z1[-1]) / (2 * m));
        total_k.append(k_y1[-1]+k_z1[-1]);
        U.append(-(m*g*l*math.cos(y1[-1])));
        Total_E.append(total_k[-1]+U[-1]);

        x1 = [];
        x2 = [];
        x3 = [];

        while t < 25:

            z1.append(z1[-1] + coeffi(g3, y1[-1], y2[-1], t, tau));

            if control1 == "n":
                y1.append(y1[-1] + coeffi(g1, y1[-1], y2[-1], t, tau));
            else:
                y1.append(y1[-1]);
            y2.append(y2[-1] + coeffi(g2, y1[-1], y2[-1], t, tau));

            z2.append(h / (math.sin(y1[-1]) * math.sin(y1[-1])));

            m_y1.append(m * l * l * y2[-1]);
            m_z1.append(m * l * l * math.sin(y1[-1]) * math.sin(y1[-1]) * z2[-1]);

            k_y1.append((m_y1[-1] * m_y1[-1]) / (2 * m));
            k_z1.append((m_z1[-1] * m_z1[-1]) / (2 * m));
            total_k.append(k_y1[-1] + k_z1[-1]);
            U.append(-(m * g * l * math.cos(y1[-1])));
            Total_E.append(total_k[-1] + U[-1]);

            t = t + tau;
            time.append(t);
            # position of the boob
            x1.append(l * math.sin(y1[-1]) * math.cos(z1[-1]));
            x2.append(l * math.sin(y1[-1]) * math.sin(z1[-1]));
            x3.append(-(l * math.cos(y1[-1])));

        return [[time, y1, y2, z1, z2], [x1, x2, x3], [m_y1, m_z1], [k_y1, k_z1, U, total_k, Total_E]];

    return R_kutta();


# representing the 3D graph
def graph3d(resultX, resultY, resultZ, title, xlabel, ylabel, zlabel):
    fig = plt.figure();
    ax = plt.gca(projection="3d");
    plt.plot(resultX, resultY, resultZ);
    plt.title(title);
    plt.xlabel(xlabel);
    plt.ylabel(ylabel);
    plt.show();

# representing the 2D graph
def graph2d(x, y, title, xlabel, ylabel):
    plt.plot(x, y);
    plt.title(title);
    plt.xlabel(xlabel);
    plt.ylabel(ylabel);
    plt.show();


# --------------------------------------------user interface-------------------------
def myProject():
    print("Welcome to the spherical pendulum graph visualization.");
    g = float(input("Enter the value of g (gravitational acceleration): "));
    l = float(input("Enter the length of the pendulum: "));
    m = float(input("Enter the mass of the bob: "));
    b = float(input("Enter the drag coefficient ( along ϑ): "));


    result = pendulum_graph(g, l, m, b, "n", "n");

    # result = pendulum_graph(g, l, m, b, 'n');

    title1 = "Graph for g = "+str(g)+", l = "+str(l)+", m= "+str(m)+", drag coefficient = "+str(b);
    graph3d(result[1][0], result[1][1], result[1][2], title1, "x position", "y position", "z position");

    i =0;
    # graph for time VS ϑ
    graph2d(result[0][0], result[0][1], title1, "Time", "value of ϑ");

    # graph for time VS ϑ'
    graph2d(result[0][0], result[0][2], title1, "Time", " velocity along ϑ");

    # graph for time VS ϕ
    graph2d(result[0][0], result[0][3], title1, "Time", "value of ϕ");

    # graph for time VS ϕ'
    graph2d(result[0][0], result[0][4], title1, "Time", "velocity along ϕ");

    # graph for ϑ VS ϑ'
    graph2d(result[0][1], result[0][2], title1, "value of ϑ", "velocity along ϑ");

    # graph for ϑ VS ϕ
    graph2d(result[0][1], result[0][3], title1, "value of ϑ", "value of ϕ");

    # graph for ϑ' VS ϕ'
    graph2d(result[0][2], result[0][4], title1, "velocity along ϑ", "velocity along ϕ");

    # graph for ϑ VS ϕ'
    graph2d(result[0][1], result[0][4], title1, "value of ϑ", "velocity along ϕ");

    # graph for ϕ VS ϑ'
    graph2d(result[0][3], result[0][2], title1, "value of ϕ", "velocity along ϑ");

    # graph for ϕ VS ϕ'
    graph2d(result[0][3], result[0][4], title1, "value of ϕ", "velocity along ϕ");

    # graph for time VS momentum along ϑ
    graph2d(result[0][0], result[2][0], title1, "time", "momentum along ϑ");

    # graph for time VS momentum along ϕ
    graph2d(result[0][0], result[2][1], title1, "time", "momentum along ϕ");

    # graph for time VS momentum ϕ
    graph2d(result[2][0], result[2][1], title1, "momentum along ϑ", "momentum along ϕ");

    # graph for time VS k along ϑ
    graph2d(result[0][0], result[3][0], title1, "Time", "kinetic energy along ϑ");

    # graph for time VS k along ϕ
    graph2d(result[0][0], result[3][1], title1, "Time", "kinetic energy along ϕ");

    # graph for k along ϑ VS k along ϕ
    graph2d(result[3][0], result[3][1], title1, "kinetic energy along ϑ", "kinetic energy along ϕ");

    # graph for time VS potential energy
    graph2d(result[0][0], result[3][2], title1, "Time", "potential energy");

    # graph for time VS Total kinetic energy
    graph2d(result[0][0], result[3][3], title1, "Time", "Total kinetic energy");

    # graph for potential energy VS Total kinetic energy
    graph2d(result[3][2], result[3][3], title1, "potential energy", "Total kinetic energy");

    # graph for time VS Total kinetic energy
    graph2d(result[0][0], result[3][4], title1, "Time", "Total energy");

    # graph for ϑ VS k along ϑ
    graph2d(result[0][1], result[3][0], title1, "value of ϑ", "k along ϑ");

    # graph for ϕ VS k along ϕ
    graph2d(result[0][2], result[3][1], title1, "value of ϕ", "k along ϕ");

    # graph for ϕ VS k along ϑ
    graph2d(result[0][2], result[3][0], title1, "value of ϕ", "k along ϑ");

    # graph for ϑ VS k along ϕ
    graph2d(result[0][1], result[3][1], title1, "value of ϑ", "k along ϕ");





    control1 = input("Is theta will be constant?: [y/n]");


    result = pendulum_graph(g, l, m, b, control1, "n");
    title = "Graph for conical pendulum where g = "+str(g)+", l = "+str(l)+", m= "+str(m)+", drag coefficient = "+str(b);
    graph3d(result[1][0], result[1][1], result[1][2], title, "x position", "y position", "z position");






myProject();

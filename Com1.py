# using python language
# python v3.6
# must install matplotlib for visulazing  graph

import matplotlib.pyplot as plt;
import math;


# Side force coefficient Cs
Cs = 0.4;
Cd = 0.5;

# This function is the y vs x
# y = 0.0077*Sc*x^2
#
# where,
#     y is the Bowling crease
#     x is the Bounce point
#     Cs is the side force coefficient (0.3)

def fx(x):
    return .0077*Cs*x*x;






x1 = [];
x2 = [];
y1 = [];
y2 = [];


i = 1.22;
i1 = 0;
while i<19.9:
    x1.append(i);
    i = i + .1;

while i1<4.12:
    x2.append(i1);
    i1 = i1 + .1;

for i in range(len(x1)):
    y1.append(fx(x1[i]));

for i in range(len(x2)):
    y2.append(fx(x2[i]));

print(x2, y2);


# before bouncing
plt.subplot(121);
plt.plot(y1,x1, "ro");
plt.xlabel("Bowling Crease");
plt.ylabel("Bounce point( Deviation before bouncing)")
plt.legend("Point of Bouncing")
plt.grid(True);

# after bouncing
plt.subplot(122);
plt.plot(x2,y2, "bo");
plt.xlabel("Stamp to interaction point");
plt.ylabel("Bounce point( Deviation after bouncing)")
plt.legend("Point of Bouncing")
plt.grid(True);
plt.show();





These modules create a family of curve types which inherit from the FinCurve class. They all provide a discount factor function which can be used to present value a future cash flow. 

Included in the classes is a credit risky curve which builds a survival probability curve on top of a risk-free curve They also rely heavily on the FinInterpolate module which provides fast interpolation.

It also includes a number of parametric curves that can be used to fit yield curves such as Nelson-Siegel.
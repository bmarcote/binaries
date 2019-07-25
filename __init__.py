# import abc

# class EmissionModels(object):
#     """ This abstract class defines what an emission model must (at least) define.
#     """
#     __metaclass__ = abc.ABCMeta

#     @abc.abstractmethod
#     def compute(self, x, *args):
#         """ Returns the flux emitted at the frequency x for the given parameters.
#         Input parameters
#         -------------------

#         Output value
#         ____________________

#         f : float
#             Flux emitted at frequency x for the model, in units of ####.
#         """
#         return NotImplemented
    
#     @abc.abstractmethod
#     def derivative(self, x, *args):
#         return NotImplemented

#     @abc.abstractmethod
#     def compute_error(self, x, *args):
#         return NotImplemented

#     @abc.abstractmethod
#     def maximum(self, x, *args):
#         return NotImplemented



# The implementation is just class(EmissionModels)



# from numpy import *

# # Basic models without physical components.
# def ssa(x, a, b):
#     return a*x**(2.5)*(1-exp(-b*x**(-3)))
# def ssa_free(x, a, b, c):
#     #return a*(x**(-0.5)/(x**(-3) + x**(-2.1)))*(1-exp(-b*(x**(-3) + x**(-2.1))))
#     return c*(1-exp(-a*x**(-3)-b*x**(-2)))/(a*x**(-2.5)+b*x**(-1.5))
# def razin(x, a, b):
#     return a*x**(-0.5)*exp(-b/x)
# def ssa_razin(x, a, b, c):
#     return a*x**(2.5)*( exp(-c/x) - exp(-b*x**(-3)-c/x) )
# def ssa_free_razin(x, a, b, c, d):
#     #return a*(x**(-0.5)/(x**(-3) + x**(-2.1)))*(exp(-c/x)-exp(-b*(x**(-3) + x**(-2.1)) - c/x))
#     return ssa_free(x, a, b, c) * exp(-d/x)





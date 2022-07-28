import numpy as np
from scipy import integrate

def NP_rhs(t, x, a, b, c, d, e):

    rhs = np.zeros(2)
    rhs[0] = a - b*x[0]*x[1] - e*x[0]
    rhs[1] = c*b*x[0]*x[1] - d*x[1]
    
    return rhs

def np_solver(t_bounds, t_eval):

    def solver(x0, args):

        def neg_population(tval, yval, *args):
            val = +1. if np.all(yval) >=0 else -1.
            return val

        neg_population.terminal = True

        result = integrate.solve_ivp(NP_rhs, t_bounds, x0, t_eval=t_eval, args=args, events=neg_population, method='RK45', atol=1e-12)

        if len(result.t_events[0]) > 0:
            assert False, "Negative population detected"

        return result

    return solver


def plot_labels(variable_interactions, variable_texts):
    """
    Computes a list of labels for plotting based on variable interactions and given variable texts.
    """

    labels = []
    for q in range(len(variable_interactions)):
        temp = list(variable_interactions[q])
        label_str = "$I = \\{$"
        nonein = True
        if 0 in temp:
            label_str += variable_texts[0]
            #label_str += "$B$"
            nonein = False

        if 1 in temp:
            if nonein:
                label_str += variable_texts[1]
                #label_str += "$N(\\tau_B)$"
                nonein = False
            else:
                label_str += ", " + variable_texts[1]
                #label_str += ", $N(\\tau_B)$"

        if 2 in temp:
            if nonein:
                label_str += variable_texts[2]
                #label_str += "$P(\\tau_B)$"
                nonein = False
            else:
                label_str += ", " + variable_texts[2]
                #label_str += ", $P(\\tau_B)$"

        label_str += "$\\}$"
        labels.append(label_str)

    return labels

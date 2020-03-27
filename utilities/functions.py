import numpy as np
import math
from scipy.interpolate import interp1d


class Function:
    def __init__(self, x_, y_):
        self.x = x_
        self.y = y_

    def __add__(self, other):
        return FunctionSum(self) + other

    def __mul__(self, other):
        return FunctionMul(self) * other

    def __rmul__(self, other):
        return FunctionMul(self) * other

    def __sum__(self):
        return Function(self.x, sum(self.y))

    def calculate(self, x):
        return self.calculate(x)


class GaussFunction(Function):
    def __init__(self, mu_, sigma_, x_=range(0, 4000)):
        self.mu = mu_
        self.sigma = sigma_
        super().__init__(x_, self.calculate_gauss(x_))

    def calculate_gauss(self, x_):
        y = 1 / (self.sigma * math.sqrt(2 * math.pi)) * np.exp(
                - (np.array(x_) - self.mu) ** 2 / 2 / self.sigma ** 2)
        return y

    def calculate(self, x):
        return self.calculate_gauss(x)


class LorenzFunction(Function):
    def __init__(self, x0_, gamma_, x_):
        self.x0 = x0_
        self.gamma = gamma_
        super().__init__(x_, self.calculate_lorenz(x_))

    def calculate_lorenz(self, x_):
        y = []
        for x_i in x_:
            y.append(1 / math.pi * (
                        self.gamma / ((x_i - self.x0) ** 2 + self.gamma ** 2)))
        return y

    def calculate(self, x):
        return self.calculate_lorenz(x)


class FunctionSum(Function):
    def __init__(self, func):
        self.func_list = [func]
        super().__init__(func.x, func.y)

    def __add__(self, other):
        if isinstance(other, Function):
            self.func_list.append(other)
            for i, y_i in enumerate(other.y):
                self.y[i] += y_i
        else:
            raise NotImplemented
        return self

    def calculate(self, x):
        sum_ = 0
        for function in self.func_list:
            sum_ += function.calculate(x)
        return sum_


class FunctionMul(Function):
    def __init__(self, func):
        self.func = func
        self.number = 1
        super().__init__(func.x, func.y)

    def __mul__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            self.number *= other
            for i, y_i in enumerate(self.y):
                self.y[i] *= other
        else:
            raise NotImplemented
        return self

    def calculate(self, x):
        return self.number * self.func.calculate(x)


class Spectra(Function):
    def __init__(self, x, y):
        super().__init__(x, y)

    def calculate(self, x):
        return interp1d(self.x, self.y)


def gauss_func(x, x_true_, mu_):
    function = 0 * GaussFunction(mu_=1, sigma_=1, x_=x_true_)
    for i in range(1, len(mu_) + 1):
        function += x[i - 1] * GaussFunction(mu_=mu_[i - 1], sigma_=x[i - 1 + len(mu_)], x_=x_true_)
    return function


def error(x, args):
    y_true = args[0]
    function = gauss_func(x, args[1], args[2])
    total_error = sum((np.array(y_true) - np.array(function.y))**2)
    return total_error


# if __name__ == '__main__':


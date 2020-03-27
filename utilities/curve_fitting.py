from utilities import functions as func
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import minimize

# ectracting peacks position
# algorithms were taken from https://www.sciencedirect.com/science/article/pii/S1386142504006316

class GaussFittingCurve:
    def __init__(self, data_x, data_y):
        self.x = np.array(data_x)
        self.f_real = np.array(data_y)
        self.a = None
        self.sigma = None
        self.mu = self.f_real
        self.i = [i for i, v in enumerate(self.x)]
        self.j = [j for j, v in enumerate(self.x)]
        self.t_max = None
        self.max_a = self.f_real.copy()
        self.peak_pos = None

    def f(self, j):
        sum_ = 0
        for i_i in self.i:
            sum_ += self.a[i_i] * np.exp(-(i_i - j) ** 2 / 2 / self.sigma[i_i] ** 2)
        return sum_

    def error(self):
        sum_ = 0
        for j in self.j:
            sum_ += (self.f(j) - self.f_real[j]) ** 2
        return 1/2 * sum_

    def d_error_d_a(self):
        result = []
        for i in self.i:
            sum_ = 0
            for j in self.j:
                sum_ += (self.f(j) - self.f_real[j]) * \
                        np.exp(-(i - j) ** 2 / 2 / self.sigma[i] ** 2)
            result.append(sum_)
        return np.array(result)

    def d_error_d_sigma(self):
        result = []
        for i in self.i:
            sum_ = 0
            for j in self.j:
                sum_ += (self.f(j) - self.f_real[j]) * \
                        np.exp(-(i - j) ** 2 / 2 / self.sigma[i] ** 2) * \
                        self.a[i] * (i - j) ** 2 / self.sigma[i] ** 3
            result.append(sum_)
        return np.array(result)

    def calculate_steps_number(self):
        # step 1
        self.a = np.random.rand(len(self.j))
        self.sigma = np.random.rand(len(self.j))
        p = 0.1
        alpha = 0.3
        delta_a = np.zeros(len(self.j))
        delta_sigma = np.zeros(len(self.j))
        t = 1
        error_max = 1

        # step 2 - 5
        while self.error() > error_max:
            delta_a = - p * self.d_error_d_a() + alpha * delta_a
            delta_sigma = - p * self.d_error_d_sigma() + alpha * delta_sigma
            self.a += delta_a
            self.sigma += delta_sigma
            t += 1

        # step 6
        self.t_max = 3 * t
        print('nit:', t)
        print('error:', self.error())

    def extract_peak_position(self, r):
        # step 1
        self.a = np.random.rand(len(self.j))
        self.sigma = np.random.rand(len(self.j))
        p = 0.1
        alpha = 0.3
        delta_a = np.zeros(len(self.j))
        delta_sigma = np.zeros(len(self.j))
        t = 1

        # step 2 - 5
        while t <= 0.4 * self.t_max:
            delta_a = - p * self.d_error_d_a() + alpha * delta_a
            delta_sigma = - p * self.d_error_d_sigma() + alpha * delta_sigma
            i_mins = []
            i_maxs = []
            for i in self.i:
                i_min = max(int(i - r * np.exp(-(1 - t/self.t_max))), 0)
                i_max = min(int(i + r * np.exp(-(1 - t/self.t_max))), len(self.i) - 1)
                # if i_min == 0:
                #     i_max = r
                i_mins.append(i_min)
                i_maxs.append(i_max)
                self.max_a[i] = max(self.f_real[i_min: i_max + 1])
                if self.f_real[i] < self.max_a[i]:
                    self.a[i] = 0
                    self.max_a[i] = 0
                else:
                    self.a[i] = self.max_a[i]
            t += 1
        self.peak_pos = [i for i, x in enumerate(self.max_a) if x != 0]
        mus = [self.x[i] for i in self.peak_pos]
        return mus

    def f2(self, j):
        sum_ = 0
        for i in self.i:
            sum_ += self.a[i] * math.exp(-(self.x[self.peak_pos[i]] - self.x[j]) ** 2 / 2 / self.sigma[i] ** 2)
        return sum_

    def error_2(self):
        sum_ = 0
        for j in self.j:
            sum_ += (self.f2(j) - self.f_real[j]) ** 2
        return 1/2 * sum_

    def d_error_d_a_2(self):
        result = []
        for i in self.i:
            sum_ = 0
            for j in self.j:
                sum_ += (self.f(j) - self.f_real[j]) * \
                        math.exp(-(self.x[j] - self.x[self.peak_pos[i]]) ** 2 / 2 / self.sigma[i] ** 2)
            result.append(sum_)
        return np.array(result)

    def d_error_d_sigma_2(self):
        result = []
        for i in self.i:
            sum_ = 0
            for j in self.j:
                sum_ += (self.f(j) - self.f_real[j]) * \
                        math.exp(-(self.x[j] - self.x[self.peak_pos[i]]) ** 2 / 2 / (self.sigma[i] ** 2)) * \
                        self.a[i] * (self.x[j] - self.x[self.peak_pos[i]]) ** 2 / self.sigma[i] ** 3
            result.append(sum_)
        return np.array(result)

    def calculate_sigma(self):
        self.peak_pos = [i for i, x in enumerate(self.max_a) if x != 0]
        self.a = np.random.rand(len(self.peak_pos))
        self.sigma = np.random.rand(len(self.a))
        self.i = [i for i, v in enumerate(self.a)]
        p = 0.1
        alpha = 0.3
        delta_a = np.zeros(len(self.a))
        delta_sigma = np.zeros(len(self.a))
        t = 1
        error_max = 0.1

        while self.error() > error_max:
            delta_a = - p * self.d_error_d_a_2() + alpha * delta_a
            delta_sigma = - p * self.d_error_d_sigma_2() + alpha * delta_sigma
            self.a += delta_a
            self.sigma += delta_sigma
            t += 1

        print('number of iterrations:', t)
        print('error:', self.error_2())
        print(self.sigma, self.a)


if __name__ == '__main__':
    # start_time = time()
    x_true = np.arange(0, 100, 2)
    data = 20 * func.GaussFunction(mu_=10, sigma_=20, x_=x_true) + 15 * func.GaussFunction(mu_=50, sigma_=20, x_=x_true) + 10 * func.GaussFunction(mu_=70, sigma_=20, x_=x_true)
    gfc = GaussFittingCurve(data.x, data.y)
    # start_calculations_time = time()
    gfc.calculate_steps_number()
    mu = gfc.extract_peak_position(5)
    print(mu)
    func_number = len(mu)
    args_ = [data.y, data.x, mu]
    point_1 = np.array([1] * 2 * func_number)
    res = minimize(func.error, args=args_, x0=point_1)
    print(res)
    plt.plot(data.x, data.y)
    plt.plot(data.x, func.gauss_func(res.x, x_true, mu).y)
    plt.show()




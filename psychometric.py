
from scipy.optimize import least_squares
from sklearn.base import BaseEstimator, RegressorMixin
import numpy as np
import pandas as pd

class PsychometricFit(BaseEstimator, RegressorMixin):
    def __init__(self):
        self.slope_ = None
        self.threshold_ = None
        self.popt_ = None
        self.r_squared_ = None

    def _sigmoid(self, x, x0, k, L=1, b=0):
        return L / (1 + np.exp(-k * (x - x0))) + b

    def _residuals(self, params, x, y):
        return y - self._sigmoid(x, *params)

    def fit(self, X, y):
        initial_guess = [np.median(X), 1.0, 1.0, 0.0]
        bounds = (
            [0, 0, 0, -1],
            [1, 500, 1, 1]
        )
        result = least_squares(self._residuals, initial_guess, bounds=bounds, args=(X, y))
        self.popt_ = result.x
        self.threshold_ = result.x[0]
        self.slope_ = result.x[1]

        residuals = y - self._sigmoid(X, *result.x)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y - np.mean(y))**2)
        self.r_squared_ = 1 - (ss_res / ss_tot)
        return self

    def predict(self, X):
        return self._sigmoid(X, *self.popt_)

    @property
    def slope(self):
        return self.slope_

    @property
    def threshold(self):
        return self.threshold_

    @property
    def r_squared(self):
        return self.r_squared_

    def summary_dataframe(self):
        data = {
            'slope': [self.slope_],
            'threshold': [self.threshold_],
            'x0': [self.popt_[0]],
            'k': [self.popt_[1]],
            'L': [self.popt_[2]],
            'b': [self.popt_[3]],
            'r_squared': [self.r_squared_]
        }
        return pd.DataFrame(data)
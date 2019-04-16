# -*- coding: utf-8 -*-


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


# Lambda types
COULOMBIC_LAMBDA = "Coulombic"
STERIC_LAMBDA = "Steric"
DUAL_LAMBDA = "Dual"
LAMBDA_TYPES = [COULOMBIC_LAMBDA, STERIC_LAMBDA, DUAL_LAMBDA]


class Lambda(object):
    def __init__(self, value, index=0, previous_lambda=None, next_lambda=None,
                 lambda_type=DUAL_LAMBDA):
        try:
            value = float(value)
        except ValueError:
            raise ValueError("Lambda Error: invalid value for a lambda: " +
                             "\'{}\'".format(value))
        if (value > 1) or (value < 0):
            raise ValueError("Lambda Error: invalid value for a lambda: " +
                             "\'{}\'. ".format(value) +
                             "It has to be between 0 and 1")

        if (lambda_type not in LAMBDA_TYPES):
            raise ValueError("Lambda Error: invalid lambda type: " +
                             "\'{}\'. ".format(lambda_type))

        self.__value = value
        self.__index = index

        self.__previous_lambda = previous_lambda
        if (previous_lambda is not None):
            previous_lambda.set_next_lambda(self)

        self.__next_lambda = next_lambda
        self.__lambda_type = lambda_type

    @property
    def value(self):
        return self.__value

    @property
    def index(self):
        return self.__index

    @property
    def previous_lambda(self):
        return self.__previous_lambda

    @property
    def next_lambda(self):
        return self.__next_lambda

    @property
    def type(self):
        return self.__lambda_type

    def __str__(self):
        return str(self.type) + ' Lambda: ' + str(self.value)

    def set_next_lambda(self, next_lambda):
        self.__next_lambda = next_lambda

    def set_previous_lambda(self, previous_lambda):
        self.__previous_lambda = previous_lambda

    def get_delta(self):
        if (self.index == 0):
            return self.value
        else:
            return abs(self.value - self.previous_lambda.value) / 2.


class IterateOverLambdas(object):
    def __init__(self, lambdas, lambda_type=DUAL_LAMBDA):
        self._lambdas = lambdas
        self._lambda_type = lambda_type
        self._current = 0
        self._previous_lambda = None

    def __iter__(self):
        return self

    def __next__(self):
        if (self._current == len(self._lambdas)):
            raise StopIteration
        else:
            self._current += 1
            new_lambda = Lambda(self._lambdas[self._current - 1],
                                self._current - 1,
                                previous_lambda=self._previous_lambda,
                                lambda_type=self._lambda_type)
            self._previous_lambda = new_lambda
            return new_lambda

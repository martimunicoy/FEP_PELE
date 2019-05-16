# -*- coding: utf-8 -*-


# Python imports
import sys
import matplotlib.pyplot as plt


# FEP_PELE imports
from .Calculators import calculateMean
from .Calculators import calculateStandardDeviationOfMean


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


class dEDistributionPlot(object):
    def __init__(self, values, averages, first_plot=0):
        self.lambdas = list(values.keys())
        self.values = list(values.values())
        self.averages = list(averages.values())
        self.total_plots = len(values)
        if (first_plot < self.total_plots):
            self.current_plot = first_plot
        else:
            self.current_plot = 0

        # Initiate plot
        self.fig = plt.figure()
        self.fig.canvas.mpl_connect('key_press_event', self._key_pressed)
        self.ax = self.fig.add_subplot(111)

        # Initiate empty args and kwargs
        self.args = []
        self.kwargs = {}

    def plotHistogram(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        average = calculateMean(self.averages[self.current_plot])
        stdev = calculateStandardDeviationOfMean(
            self.averages[self.current_plot])

        self._setTitle(average, stdev)
        self._setAxisLabels()
        self.ax.axvspan(average - stdev, average + stdev, alpha=0.5,
                        color='red')
        self.ax.axvline(x=average, linewidth=2, color='r', linestyle="--")
        self.ax.hist(self.values[self.current_plot], *self.args,
                     **self.kwargs)

        plt.show()

    def _setTitle(self, average, stdev):
        plt.suptitle(r'$\Delta E$ distribution', fontsize=12)
        plt.title("({} of {})".format(self.current_plot + 1,
                                      self.total_plots), fontsize=10)
        title = self.ax.text(0.85, 0.80, "",
                             bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 3},
                             transform=self.ax.transAxes, ha="center")
        average_box = self.ax.text(
            0.20, 0.85, "", bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 3},
            transform=self.ax.transAxes, ha="center", color='red')

        lambda_shift = self.lambdas[self.current_plot][:2]
        lambda_type = self.lambdas[self.current_plot][2]

        title.set_text(lambda_type + '\n' +
                       r'$\lambda_{0} =$' +
                       "{}".format(lambda_shift[0]) + '\n' +
                       r'$\lambda_{1} =$' +
                       "{}".format(lambda_shift[1]))
        average_box.set_text(r'$\Delta E = $' +
                             "{: .2f}".format(round(average, 2)) +
                             r' $kcal/mol$' + '\n' +
                             r'$\sigma = $' +
                             "{:.3f}".format(round(stdev, 3)))

    def _setAxisLabels(self):
        plt.xlabel(r'$\Delta E$')
        plt.ylabel('Frequency')

    def _key_pressed(self, event):
        if (event.key == "left"):
            self._previous(event)
        elif (event.key == "right"):
            self._next(event)

    def _previous(self, event):
        if (self.current_plot > 0):
            self.current_plot -= 1
            average = calculateMean(self.averages[self.current_plot])
            stdev = calculateStandardDeviationOfMean(
                self.averages[self.current_plot])

            self.ax.clear()

            self._setTitle(average, stdev)
            self._setAxisLabels()
            self.ax.axvspan(average - stdev, average + stdev, alpha=0.5,
                            color='red')
            self.ax.axvline(x=average, linewidth=2, color='r', linestyle="--")
            self.ax.hist(self.values[self.current_plot], *self.args,
                         **self.kwargs)

            plt.draw()
        else:
            sys.stdout.write('\a')
            sys.stdout.flush()

    def _next(self, event):
        if (self.current_plot < self.total_plots - 1):
            self.current_plot += 1
            average = calculateMean(self.averages[self.current_plot])
            stdev = calculateStandardDeviationOfMean(
                self.averages[self.current_plot])

            self.ax.clear()

            self._setTitle(average, stdev)
            self._setAxisLabels()
            self.ax.axvspan(average - stdev, average + stdev, alpha=0.5,
                            color='red')
            self.ax.axvline(x=average, linewidth=2, color='r', linestyle="--")
            self.ax.hist(self.values[self.current_plot], *self.args,
                         **self.kwargs)

            plt.draw()
        else:
            sys.stdout.write('\a')
            sys.stdout.flush()


class dEVariationPlot(object):
    def __init__(self, values, first_plot=0):
        self.lambdas = list(values.keys())
        self.values = list(values.values())
        self.total_plots = len(values)
        if (first_plot < self.total_plots):
            self.current_plot = first_plot
        else:
            self.current_plot = 0

        # Initiate plot
        self.fig = plt.figure()
        self.fig.canvas.mpl_connect('key_press_event', self._key_pressed)
        self.ax = self.fig.add_subplot(111)

        # Initiate empty args and kwargs
        self.args = []
        self.kwargs = {}

    def plotGraph(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

        self._setTitle()
        self._setAxisLabels()
        self.ax.plot(self.values[self.current_plot], *self.args,
                     **self.kwargs)

        plt.show()

    def _setTitle(self):
        plt.suptitle(r'$\Delta E$ variation', fontsize=12)
        plt.title("({} of {})".format(self.current_plot + 1,
                                      self.total_plots), fontsize=10)
        title = self.ax.text(0.85, 0.80, "",
                             bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 3},
                             transform=self.ax.transAxes, ha="center")

        lambda_shift = self.lambdas[self.current_plot][:2]
        lambda_type = self.lambdas[self.current_plot][2]

        title.set_text(lambda_type + '\n' +
                       r'$\lambda_{0} =$' +
                       "{}".format(lambda_shift[0]) + '\n' +
                       r'$\lambda_{1} =$' +
                       "{}".format(lambda_shift[1]))

    def _setAxisLabels(self):
        plt.ylabel('PELE accepted steps')
        plt.ylabel(r'$\Delta E$')

    def _key_pressed(self, event):
        if (event.key == "left"):
            self._previous(event)
        elif (event.key == "right"):
            self._next(event)

    def _previous(self, event):
        if (self.current_plot > 0):
            self.current_plot -= 1

            self.ax.clear()

            self._setTitle()
            self._setAxisLabels()

            self.ax.plot(self.values[self.current_plot], *self.args,
                         **self.kwargs)

            plt.draw()
        else:
            sys.stdout.write('\a')
            sys.stdout.flush()

    def _next(self, event):
        if (self.current_plot < self.total_plots - 1):
            self.current_plot += 1

            self.ax.clear()

            self._setTitle()
            self._setAxisLabels()

            self.ax.plot(self.values[self.current_plot], *self.args,
                         **self.kwargs)

            plt.draw()

            plt.draw()
        else:
            sys.stdout.write('\a')
            sys.stdout.flush()

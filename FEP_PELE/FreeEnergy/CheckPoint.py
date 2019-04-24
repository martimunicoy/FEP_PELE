# -*- coding: utf-8 -*-


# Python imports
import sys
import pickle


# FEP_PELE imports
from FEP_PELE.Utils.InOut import isThereAFile


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


class CheckPoint(object):
    def __init__(self, path, settings):
        self._path = path
        self._settings = settings
        self._restart_available = False
        self._completed_commands = {}
        self._initialized = False
        self._data = None

    @property
    def restart_available(self):
        return self._restart_available

    @property
    def initialized(self):
        return self._initialized

    def initialize(self):
        self._initialized = True

        if (not self._settings.restart):
            self._initiateCheckPointFile()

        if (not isThereAFile(self._path)):
            self._initiateCheckPointFile()

        try:
            self._loadCheckPointFile()
        except EOFError:
            self._initiateCheckPointFile()
            self._loadCheckPointFile()

    def _initiateCheckPointFile(self):
        with open(self._path, 'wb') as file:
            data = {'Settings': str(self._settings)}
            pickle.dump(data, file)

    def _loadCheckPointFile(self):
        with open(self._path, 'rb') as file:
            data = pickle.load(file)

            if (data['Settings'] == str(self._settings)):
                self._restart_available = True
                self._data = data
            else:
                raise RuntimeError("settings in the current checkpoint " +
                                   "do not match with settings from the " +
                                   "input file. Restart events will not " +
                                   "take place.")

    def save(self, checkPointData):
        if (not self.initialized):
            print("CheckPoint Error: called save() before initializing" +
                  "the checkpoint")
            sys.exit(1)

        command_to_save, checkpoint = checkPointData

        for command, checkpoints in self._data.items():
            if (command == command_to_save):
                checkpoints.add(checkpoint)
                break
        else:
            self._data[command_to_save] = set([checkpoint, ])

        with open(self._path, 'wb') as file:
            pickle.dump(self._data, file)

    def check(self, checkPointData):
        if (not self.initialized):
            print("CheckPoint Error: called check() before initializing " +
                  "the checkpoint")
            sys.exit(1)
        if (not self.restart_available):
            return False

        if (not self._settings.restart):
            return False

        command_to_save, checkpoint = checkPointData

        with open(self._path, 'rb') as file:
            data = pickle.load(file)

        for command, checkpoints in data.items():
            if (command == command_to_save):
                return checkpoint in checkpoints

        return False

# -*- coding: utf-8 -*-


# Python imports
import os


# FEP_Pele imports
from . import Constants as co
from FEP_PELE.Utils.InOut import checkPath, checkFile
from FEP_PELE.Utils.InOut import getFileFromPath
from FEP_PELE.Tools.StringTools import asPath
from FEP_PELE.Tools.StringTools import asBool


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


# Class definitions
class Settings(object):
    def __init__(self):
        # Setting default values
        self.__general_path = asPath(os.getcwd())
        self.__serial_pele = co.DEF_SERIAL_PELE
        self.__mpi_pele = co.DEF_MPI_PELE
        self.__initial_template = co.DEF_INITIAL_TEMPLATE
        self.__final_template = co.DEF_FINAL_TEMPLATE
        self.__atom_links = co.DEF_ATOM_LINKS
        self.__lambdas = co.DEF_LAMBDAS
        self.__number_of_processors = co.DEF_NUMBER_OF_PROCESSORDS
        self.__commands = co.DEF_COMMANDS
        self.__min_control_file = co.DEF_MIN_CONTROL_FILE
        self.__sim_control_file = co.DEF_SIM_CONTROL_FILE
        self.__pp_control_file = co.DEF_PP_CONTROL_FILE
        self.__sp_control_file = co.DEF_SP_CONTROL_FILE
        self.__safety_check = co.DEF_SAFETY_CHECK
        self.__minimization_folder = asPath(co.DEF_MIN_FOLDER)
        self.__simulation_folder = asPath(co.DEF_SIM_FOLDER)
        self.__calculation_folder = asPath(co.DEF_CAL_FOLDER)

        # Other
        self.__default_lambdas = True

    @property
    def general_path(self):
        return self.__general_path

    @property
    def serial_pele(self):
        return self.__serial_pele

    @property
    def mpi_pele(self):
        return self.__mpi_pele

    @property
    def initial_template(self):
        return self.__initial_template

    @property
    def final_template(self):
        return self.__final_template

    @property
    def atom_links(self):
        return self.__atom_links

    @property
    def lambdas(self):
        return sorted(self.__lambdas)

    @property
    def number_of_processors(self):
        return self.__number_of_processors

    @property
    def command_names(self):
        return self.__commands

    @property
    def min_control_file(self):
        return self.__min_control_file

    @property
    def sim_control_file(self):
        return self.__sim_control_file

    @property
    def pp_control_file(self):
        return self.__pp_control_file

    @property
    def sp_control_file(self):
        return self.__sp_control_file

    @property
    def safety_check(self):
        return self.__safety_check

    @property
    def final_template_name(self):
        return getFileFromPath(self.__final_template)

    @property
    def minimization_folder(self):
        return self.__minimization_folder

    @property
    def minimization_path(self):
        return self.general_path + self.minimization_folder

    @property
    def simulation_folder(self):
        return self.__simulation_folder

    @property
    def simulation_path(self):
        return self.general_path + self.simulation_folder

    @property
    def calculation_folder(self):
        return self.__calculation_folder

    @property
    def calculation_path(self):
        return self.general_path + self.calculation_folder

    def set(self, key, value):
        if (key == co.CONTROL_FILE_DICT["GENERAL_PATH"]):
            value = self._getSingleValue(key, value)
            self._checkPath(key, value)
            self.__general_path = asPath(str(value))

        elif (key == co.CONTROL_FILE_DICT["SERIAL_PELE_PATH"]):
            value = self._getSingleValue(key, value)
            self._checkFile(key, value)
            self.__serial_pele = str(value)

        elif (key == co.CONTROL_FILE_DICT["MPI_PELE_PATH"]):
            value = self._getSingleValue(key, value)
            self._checkFile(key, value)
            self.__mpi_pele = str(value)

        elif (key == co.CONTROL_FILE_DICT["INITIAL_TEMPLATE"]):
            value = self._getSingleValue(key, value)
            value = self._checkFile(key, value)
            self.__initial_template = str(value)

        elif (key == co.CONTROL_FILE_DICT["FINAL_TEMPLATE"]):
            value = self._getSingleValue(key, value)
            value = self._checkFile(key, value)
            self.__final_template = str(value)

        elif (key == co.CONTROL_FILE_DICT["ATOM_LINK"]):
            self._checkAtomLink(key, value)
            self.__atom_links.append((str(value[0]), str(value[2])))

        elif (key == co.CONTROL_FILE_DICT["LAMBDAS"]):
            value = self._getCommaSeparatedList(key, value)
            self._checkListOfLambdas(key, value)
            if (self.__default_lambdas):
                self.__lambdas = set()
                self.__default_lambdas = False

            for lambda_chunk in value:
                for lambda_value in lambda_chunk.split(','):
                    self.__lambdas.add(float(lambda_value))

        elif (key == co.CONTROL_FILE_DICT["NUMBER_OF_PROCESSORS"]):
            value = self._getSingleValue(key, value)
            self._checkPositiveInteger(key, value)
            self.__number_of_processors = int(value)

        elif (key == co.CONTROL_FILE_DICT["COMMANDS"]):
            value = self._getCommaSeparatedList(key, value)
            self._checkCommandNames(key, value)
            for command_chunk in value:
                for command in command_chunk.split(','):
                    self.__commands.append(str(command))

        elif (key == co.CONTROL_FILE_DICT["MIN_CONTROL_FILE"]):
            value = self._getSingleValue(key, value)
            value = self._checkFile(key, value)
            self.__min_control_file = str(value)

        elif (key == co.CONTROL_FILE_DICT["SIM_CONTROL_FILE"]):
            value = self._getSingleValue(key, value)
            value = self._checkFile(key, value)
            self.__sim_control_file = str(value)

        elif (key == co.CONTROL_FILE_DICT["PP_CONTROL_FILE"]):
            value = self._getSingleValue(key, value)
            value = self._checkFile(key, value)
            self.__pp_control_file = str(value)

        elif (key == co.CONTROL_FILE_DICT["SP_CONTROL_FILE"]):
            value = self._getSingleValue(key, value)
            value = self._checkFile(key, value)
            self.__sp_control_file = str(value)

        if (key == co.CONTROL_FILE_DICT["MIN_FOLDER"]):
            value = self._getSingleValue(key, value)
            self.__minimization_folder = asPath(str(value))

        if (key == co.CONTROL_FILE_DICT["SIM_FOLDER"]):
            value = self._getSingleValue(key, value)
            self.__simulation_folder = asPath(str(value))

        if (key == co.CONTROL_FILE_DICT["CAL_FOLDER"]):
            value = self._getSingleValue(key, value)
            self.__calculation_folder = asPath(str(value))

        elif (key == co.CONTROL_FILE_DICT["SAFETY_CHECK"]):
            value = self._getSingleValue(key, value)
            value = self._checkBool(key, value)
            self.__safety_check = value

    def _getSingleValue(self, key, value):
        if (len(value) != 1):
            raise NameError('Expected a single value in line with key' +
                            '{}. '.format(key) + 'Obtained {}'.format(value))

        return value[0]

    def _getCommaSeparatedList(self, key, value):
        list_of_chunks = [i.replace(',', ' ') for i in value]

        list_of_strings = []
        for chunk in list_of_chunks:
            for item in chunk.split():
                list_of_strings.append(item)

        return list_of_strings

    def _checkBool(self, key, value):
        try:
            value = asBool(value)
        except ValueError as exception:
            print("Error while setting \'{}\',\'{}\'".format(key, value) +
                  ": " + str(exception))
            exit(1)

        return value

    def _checkPath(self, key, value):
        try:
            checkPath(asPath(str(value)))
        except NameError as exception:
            print("Error while setting \'{}\',\'{}\'".format(key, value) +
                  ": " + str(exception))
            exit(1)

    def _checkFile(self, key, value):
        try:
            checkFile(str(value))
        except NameError:
            try:
                checkFile(self.general_path + str(value))
            except NameError as exception:
                print("Error while setting \'{}\',\'{}\'".format(key, value) +
                      ": " + str(exception))
                exit(1)
            return self.general_path + str(value)

        return value

    def _checkAtomLink(self, key, value):
        okay = True
        message = ""

        if (len(value) != 3):
            okay = False
            message += "Incorrect number of fields detected. "

        else:
            atom1 = value[0]
            atom2 = value[2]
            if (len(atom1) != 4) or (len(atom2) != 4):
                okay = False
                message += "Incorrect atom name, it must have 4 characters. "

        if (not okay):
            print("Error while setting \'{}\',\'{}\'".format(key, value) +
                  ": " + message)
            exit(1)

    def _checkListOfLambdas(self, key, value):
        okay = True
        message = ""

        for lambda_value in value:
            try:
                lambda_value = float(lambda_value)
            except ValueError:
                okay = False
                message += "Lambda values need to be floats. "
                break

            if (lambda_value > 1) or (lambda_value < 0):
                okay = False
                message += "Lambda value out of range [0-1]. "
                break

        if (not okay):
            print("Error while setting \'{}\',\'{}\'".format(key, value) +
                  ": " + message)
            exit(1)

    def _checkPositiveInteger(self, key, value):
        okay = True
        message = ""

        try:
            value = int(value)
        except ValueError:
            okay = False
            message += "Input value is not an integer. "

        if (value < 1):
            okay = False
            message += "Input value is not a positive integer. "

        if (not okay):
            print("Error while setting \'{}\',\'{}\'".format(key, value) +
                  ": " + message)
            exit(1)

    def _checkCommandNames(self, key, value):
        okay = True
        message = ""

        for command in value:
            if (command not in co.COMMAND_NAMES_LIST):
                okay = False
                message += "Command not recogniced. "

        if (not okay):
            print("Error while setting \'{}\',\'{}\'".format(key, value) +
                  ": " + message)
            exit(1)

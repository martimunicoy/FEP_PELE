import Constants as co
from ..Utils.InOut import checkPath


class Settings(object):
    def __init__(self):
        # Setting default values
        self.__general_path = co.DEF_GENERAL_PATH
        self.__initial_template = co.DEF_INITIAL_TEMPLATE
        self.__final_template = co.DEF_FINAL_TEMPLATE
        self.__atom_links = co.DEF_ATOM_LINKS
        self.__lambdas = co.DEF_LAMBDAS
        self.__number_of_processors = co.DEF_NUMBER_OF_PROCESSORDS
        self.__commands = co.DEF_COMMANDS

    @property
    def general_path(self):
        return self.__general_path

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
        return self.__lambdas

    @property
    def number_of_processors(self):
        return self.__number_of_processors

    @property
    def commands(self):
        return self.__commands

    def set(self, key, value):
        if (key == co.CONTROL_FILE_DICT["GENERAL_PATH"]):
            self._checkPath(key, value)
            self.__general_path = value

        elif (key == co.CONTROL_FILE_DICT["INITIAL_TEMPLATE"]):
            self._checkPath(key, value)
            self.__initial_template = value

        elif (key == co.CONTROL_FILE_DICT["FINAL_TEMPLATE"]):
            self._checkPath(key, value)
            self.__final_template = value

        elif (key == co.CONTROL_FILE_DICT["ATOM_LINK"]):
            self._checkAtomLink(key, value)
            self.__atom_links.append((value[0], value[2]))

        elif (key == co.CONTROL_FILE_DICT["LAMBDAS"]):
            self._checkListOfLambdas(key, value)
            for lambda_value in value:
                self.__lambdas.append(float(lambda_value))

        elif (key == co.CONTROL_FILE_DICT["LAMBDAS"]):
            self._checkPositiveInteger(key, value)
            for lambda_value in value:
                self.__lambdas.append(float(lambda_value))

    def _checkPath(self, key, value):
        try:
            checkPath(value)
        except NameError as exception:
            print("Error while setting \'{}\',\'{}\'".format(key, value) +
                  ": " + exception.message)
            exit(1)

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
            message += "Number of processors needs to be an integer. "

        if (value < 1):
            okay = False
            message += "Number of processors can only be a positive integer. "

        if (not okay):
            print("Error while setting \'{}\',\'{}\'".format(key, value) +
                  ": " + message)
            exit(1)

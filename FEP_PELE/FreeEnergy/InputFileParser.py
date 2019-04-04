from . import Constants as co
from .Settings import Settings
from FEP_PELE.Utils.InOut import isThereAFile


class InputFileParser(object):
    def __init__(self, path):
        if (not isThereAFile(path)):
            raise IOError('FEP_PELE input file not found')

        self.path = path

    def createSettings(self):
        settings = Settings()

        with open(self.path, 'r') as input_file:
            # If GeneralPath is defined, it needs to be set first
            for line in input_file:
                line = line.strip()
                if (line.startswith('#')) or (line == ''):
                    continue
                key, value = parseLine(line)
                if (key != co.CONTROL_FILE_DICT["GENERAL_PATH"]):
                    continue
                settings.set(key, value)

            # Go to first line again
            input_file.seek(0)

            # Set the rest of the settings
            for line in input_file:
                line = line.strip()
                if (line.startswith('#')) or (line == ''):
                    continue
                key, value = parseLine(line)
                if (key == co.CONTROL_FILE_DICT["GENERAL_PATH"]):
                    continue
                settings.set(key, value)

        return settings


def parseLine(line):
    fields = line.split()

    okay = True

    if (len(fields) < 2):
        okay = False

    key = fields[0].strip(':')

    if (key not in co.INPUT_FILE_KEYS):
        okay = False

    if (not okay):
        raise NameError('Line \'{}\' has an invalid format'.format(line))

    return key, fields[1:]

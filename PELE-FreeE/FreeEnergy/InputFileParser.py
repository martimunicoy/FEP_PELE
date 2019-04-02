import Constants as co
from .Settings import Settings
from ..Utils.InOut import isThereAPath


class InputFileParser(object):
    def __init__(self, path):
        if (not isThereAPath(path)):
            raise IOError('PELE-FreeE input file not found')

        self.path = path

    def createSettings(self):
        settings = Settings()

        with open(self.path, 'r') as input_file:
            for line in input_file:
                if (line.startswith('#')):
                    continue
                key, value = parseLine(line)
                settings.set(key, value)

        return settings


def parseLine(line):
    line = line.strip()
    fields = line.split()

    okay = True

    if (fields < 2):
        okay = False

    if (fields[0] not in co.INPUT_FILE_KEYS):
        okay = False

    if (not okay):
        raise NameError('Line \'{}\' has an invalid format'.format(line))

    return fields[0], fields[1:]

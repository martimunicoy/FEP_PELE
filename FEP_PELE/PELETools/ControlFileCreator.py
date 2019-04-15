# -*- coding: utf-8 -*-


# Python imports
import sys


# FEP_PELE imports
from FEP_PELE.PELETools import PELEConstants as pele_co

from FEP_PELE.Utils.InOut import checkFile
from FEP_PELE.Utils.InOut import checkPath
from FEP_PELE.Utils.InOut import getPathFromFile


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


# Private constants
CONTROL_FILE_FLAGS = ["$LICENSE_PATH$",
                      "$LOG_PATH$",
                      "$REPORT_PATH$",
                      "$TRAJECTORY_PATH$",
                      "$SOLVENT_TYPE$",
                      "$INPUT_PDB_NAME$",
                      "$SEED$"]

FROM_FLAGS_TO_NAMES = {
    CONTROL_FILE_FLAGS[0]: pele_co.CONTROL_FILE_FLAG_NAMES[0],
    CONTROL_FILE_FLAGS[1]: pele_co.CONTROL_FILE_FLAG_NAMES[1],
    CONTROL_FILE_FLAGS[2]: pele_co.CONTROL_FILE_FLAG_NAMES[2],
    CONTROL_FILE_FLAGS[3]: pele_co.CONTROL_FILE_FLAG_NAMES[3],
    CONTROL_FILE_FLAGS[4]: pele_co.CONTROL_FILE_FLAG_NAMES[4],
    CONTROL_FILE_FLAGS[5]: pele_co.CONTROL_FILE_FLAG_NAMES[5],
    CONTROL_FILE_FLAGS[6]: pele_co.CONTROL_FILE_FLAG_NAMES[6]}

FROM_NAMES_TO_FLAGS = {
    pele_co.CONTROL_FILE_FLAG_NAMES[0]: CONTROL_FILE_FLAGS[0],
    pele_co.CONTROL_FILE_FLAG_NAMES[1]: CONTROL_FILE_FLAGS[1],
    pele_co.CONTROL_FILE_FLAG_NAMES[2]: CONTROL_FILE_FLAGS[2],
    pele_co.CONTROL_FILE_FLAG_NAMES[3]: CONTROL_FILE_FLAGS[3],
    pele_co.CONTROL_FILE_FLAG_NAMES[4]: CONTROL_FILE_FLAGS[4],
    pele_co.CONTROL_FILE_FLAG_NAMES[5]: CONTROL_FILE_FLAGS[5],
    pele_co.CONTROL_FILE_FLAG_NAMES[6]: CONTROL_FILE_FLAGS[6]}


# Class definitions
class ControlFileFromTemplateCreator(object):
    def __init__(self, template_path):
        try:
            checkFile(template_path)
        except NameError:
            print("Error: invalid path to Control File template: " +
                  "{}".format(template_path))
            sys.exit(1)

        self._template_path = template_path

        self._text = self.read()

        self._flags_to_replace = self.findFlags()

    def read(self):
        with open(self._template_path, 'r') as file:
            return file.read()

    def findFlags(self):
        flags_to_replace = set()

        for line in self._text.split('\n'):
            for flag in CONTROL_FILE_FLAGS:
                if flag in line:
                    flags_to_replace.add(FROM_FLAGS_TO_NAMES[flag])

        return flags_to_replace

    def replaceFlag(self, flag, value):
        if (flag not in pele_co.CONTROL_FILE_FLAG_NAMES):
            print("  - Warning: unknown flag {}".format(flag))

        new_text = ""

        for line in self._text.split('\n'):
            new_text += line.replace(FROM_NAMES_TO_FLAGS[flag], str(value))
            new_text += '\n'

        self._text = new_text

        self._flags_to_replace = self.findFlags()

    def write(self, output_path):
        try:
            checkPath(getPathFromFile(output_path))
        except NameError:
            print("Error: invalid output path to write customized Control " +
                  "File: {}".format(output_path))
            sys.exit(1)

        if (len(self._flags_to_replace) != 0):
            print("Error: there are still flags that were not replaced in " +
                  "the original templatized Control File: " +
                  "{}".format(self._flags_to_replace))
            sys.exit(1)

        with open(output_path, 'w') as file:
            file.write(self._text)

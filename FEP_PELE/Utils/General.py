# -*- coding: utf-8 -*-


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


# Function definitions
def fromDictValuesToList(input_dict):
    output_list = []
    for key, value in input_dict.items():
        if type(value) is dict:
            value = fromDictValuesToList(value)
        for element in value:
            output_list.append(element)
    return output_list

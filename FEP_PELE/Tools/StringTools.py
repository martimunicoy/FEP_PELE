# -*- coding: utf-8 -*-


# Python imports
import re


# FEP_PELE imports
from . import Constants as co


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


def natural_sort(l):
    def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
        return [int(text) if text.isdigit() else text.lower()
                for text in _nsre.split(s)]

    return sorted(l, key=natural_sort_key)


def asPath(path):
    if (path[-1] != '/'):
        path += '/'
    return path


def asBool(value):
    if (value.lower() in co.TRUE_STRING_BOOLEANS):
        return True
    elif (value.lower() in co.FALSE_STRING_BOOLEANS):
        return False
    else:
        raise ValueError("could not convert {} to boolean. ".format(value) +
                         "To define true booleans, you can use " +
                         "{}. ".format(co.TRUE_STRING_BOOLEANS) +
                         "To define false booleans, you can use " +
                         "{}. ".format(co.FALSE_STRING_BOOLEANS))

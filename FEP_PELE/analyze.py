# -*- coding: utf-8 -*-


# Python imports
import argparse


# FEP_PELE imports
from FreeEnergy.InputFileParser import InputFileParser
from FreeEnergy.CommandTypes.EmptyCommand import EmptyCommand
from FreeEnergy.Analysis import FEPAnalysis

from FEP_PELE.Utils.InOut import printCommandTitle


# Function definitions
def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', metavar='PATH', type=str, nargs=1,
                        help='Path to input file')
    parser.add_argument('-d', '--divisions', metavar='INTEGER', type=int,
                        default=1, help='Number of divisions to calculate ' +
                        'standard deviation')

    args = parser.parse_args()

    path_to_input_file = args.input_file[0]
    divisions = args.divisions

    return path_to_input_file, divisions


def main():
    path_to_input_file, divisions = parseArguments()

    inputFileParser = InputFileParser(path_to_input_file)
    settings = inputFileParser.createSettings()

    printCommandTitle("FEP-PELE Analysis Script")

    com = EmptyCommand(settings, path=settings.calculation_path)

    print(" - Retrieving lambda folders from {}".format(com.path))

    lambda_folders = com._getLambdaFolders()

    print(" - Analyzing results")

    analysis = FEPAnalysis.FEPAnalysis(lambda_folders,
                                       com.settings.sampling_method,
                                       divisions=divisions)

    print(" - Plotting energetic histogram")

    analysis.plotHistogram()

    analysis.printResults()


# Set an executable behaviour
if __name__ == '__main__':
    main()

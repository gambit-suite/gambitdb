# a class which takes in 2 spreadsheets, assembly_metadata_spreadsheet and gambit_results_file, then compares them with pandas
# so that we can work out how good gambit is at recalling species. The column they are joined on is the 
# accession number which is the first column in both spreadsheets. The species column from the assembly_metadata_spreadsheet needs to be compared to the equivalent predicted.name column in the gambit_results_file.
#The spreadsheets are in the following formats.
#assembly_metadata_spreadsheet file format:
# uuid,species_taxid,assembly_accession,species
# GCA_000277025.1,1,GCA_000277025.1,Legionella pneumophila
# GCA_000277065.1,1,GCA_000277065.1,Legionella pneumophila
# GCA_000300315.1,2,GCA_000300315.1,Coxiella burnetii
# GCA_000359545.5,2,GCA_000359545.5,Coxiella burnetii
#
#gambit_results_file file format:
# query,predicted.name,predicted.rank,predicted.ncbi_id,predicted.threshold,closest.distance,closest.description,next.name,next.rank,next.ncbi_id,next.threshold
# GCA_000277025.1,Legionella pneumophila,species,1,1.0,0.0,,,,,
# GCA_000277065.1,Legionella pneumophila,species,1,1.0,0.0,,,,,
# GCA_000300315.1,Coxiella burnetii,species,2,0.9998,0.0,,,,,
# GCA_000359545.5,Coxiella burnetii,species,2,0.9998,0.0,,,,,

import pandas as pd
import logging

class DatabaseRecall:
    def __init__(self, assembly_metadata_spreadsheet, gambit_results_file, output_filename, debug, verbose):
        self.assembly_metadata_spreadsheet = assembly_metadata_spreadsheet
        self.gambit_results_file = gambit_results_file
        self.output_filename = output_filename
        self.debug = debug
        self.verbose = verbose

        self.logger = logging.getLogger('DatabaseRecall')
        self.logger.debug('DatabaseRecall.__init__')

    def compare_results(self):
        self.logger.debug('DatabaseRecall.compare_results')
        # Read in the assembly metadata spreadsheet
        assembly_metadata = pd.read_csv(self.assembly_metadata_spreadsheet)
        # Read in the gambit results file
        gambit_results = pd.read_csv(self.gambit_results_file)

        # Join the two spreadsheets on the accession number
        joined = pd.merge(assembly_metadata, gambit_results, left_on='assembly_accession', right_on='query')
        # Compare the species column from the assembly_metadata spreadsheet to the predicted.name column in the gambit_results file
        joined['correct'] = joined['species'] == joined['predicted.name']

        correct = 0
        incorrect = 0
        # Count the number of correct and incorrect predictions
        # check if joined['correct'] contains a true value
        if joined['correct'].any():
            correct = joined['correct'].value_counts()[True]
        
        if not joined['correct'].any():
            incorrect = joined['correct'].value_counts()[False]

        # Print the number of correct and incorrect predictions
        print('Correct predictions: ' + str(correct))
        print('Incorrect predictions: ' + str(incorrect))
        # Calculate the percentage of correct predictions
        percentage_correct = correct/(correct+incorrect)*100
        # Print the percentage of correct predictions
        print('Percentage correct: ' + str(percentage_correct) + '%')

        # Print the number of correct and incorrect predictions to a file
        with open(self.output_filename, 'w') as f:
            f.write('Correct predictions: ' + str(correct) + '\n')
            f.write('Incorrect predictions: ' + str(incorrect) + '\n')
            f.write('Percentage correct: ' + str(percentage_correct) + '%\n')

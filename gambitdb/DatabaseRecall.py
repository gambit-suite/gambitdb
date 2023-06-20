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
    """
    Compares two spreadsheets to work out how good Gambit is at recalling species.
    """
    def __init__(self, assembly_metadata_spreadsheet, gambit_results_file, output_filename, debug, verbose):
        """
    Initializes the DatabaseRecall class.
    Args:
      assembly_metadata_spreadsheet (str): The path to the assembly metadata spreadsheet.
      gambit_results_file (str): The path to the Gambit results file.
      output_filename (str): The path to the output file.
      debug (bool): Whether to enable debug logging.
      verbose (bool): Whether to enable verbose logging.
    Side Effects:
      Initializes the logger.
    """
        self.assembly_metadata_spreadsheet = assembly_metadata_spreadsheet
        self.gambit_results_file = gambit_results_file
        self.output_filename = output_filename
        self.debug = debug
        self.verbose = verbose

        self.logger = logging.getLogger('DatabaseRecall')
        self.logger.debug('DatabaseRecall.__init__')

    def compare_results(self):
        """
    Compares the species column from the assembly_metadata_spreadsheet to the predicted.name column in the gambit_results_file.
    Args:
      None
    Returns:
      None
    Side Effects:
      Prints the number of correct and incorrect predictions.
      Writes the output to the output_filename.
    Examples:
      >>> compare_results()
      Indentical predictions: 3
      Percentage identical: 75.0%
    """
        self.logger.debug('DatabaseRecall.compare_results')
        # Read in the assembly metadata spreadsheet
        assembly_metadata = pd.read_csv(self.assembly_metadata_spreadsheet)
        # Read in the gambit results file
        gambit_results = pd.read_csv(self.gambit_results_file)
        num_samples = gambit_results.shape[0]

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
        print('Identical predictions: ' + str(correct))
        # Calculate the percentage of correct predictions
        percentage_correct = correct/(num_samples)*100
        # Print the percentage of correct predictions
        print('Percentage identical: ' + str(percentage_correct) + '%')

        output_df = joined[['species', 'predicted.name', 'assembly_accession']]
        output_df.to_csv(self.output_filename, index=False)

        # select rows where correct is True
        incorrect_df = joined[joined['correct'] == False]
        incorrect_df = incorrect_df[['species', 'predicted.name', 'assembly_accession']]
        # sort by species
        incorrect_df = incorrect_df.sort_values(by=['species'])
        incorrect_df.to_csv(self.output_filename + '.differences.csv', index=False)
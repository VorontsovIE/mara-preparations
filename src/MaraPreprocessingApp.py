#!/usr/bin/env python
import os
import sys
import json
import re
import subprocess
import argparse
import shutil
import logging
import math

BASE_DIR = os.getcwd()

# All available parameters
USER_INPUT = {
    'tss_clusters_file',
}

SOURCE_DATA = {
    'genome_file',
    'motif_dir'
}

REQUIRED_PARAMETERS = {
    'flank_5',
    'flank_3'
}

OPTIONAL_PARAMETERS = {
    'custom_motifs',
    'scoring_mode',
    'num_processes',
    'log_level',
    'stage',
    'motif'
}

AVAILABLE_FLANK_PAIRS = {
    '250u 0u',
    '250u 0d', # Handle both 0u and 0d
    '250u 10d',
    '250u 20d',
    '10u 50d',
    '0u 50d',
    '0d 50d', # Handle both 0u and 0d
    '10d 50d'
}


ALL_PARAMETERS = USER_INPUT | SOURCE_DATA | REQUIRED_PARAMETERS | OPTIONAL_PARAMETERS


class ConfigError(Exception):
    pass


class ArgumentError(Exception):
    pass


class ExternalToolError(Exception):
    pass


class CommandExecutionError(Exception):
    pass


class PwmConverter:
    """
    A class to convert motif matrices between PCM (Position Count Matrix),
    PFM (Position Frequency Matrix), and PWM (Position Weight Matrix) formats.

    Attributes:
        filename (str): The input file name.
        pseudocount (str or float): The pseudocount method ('log', 'sqrt') or a numeric value.
        word_count (float): The word count used for PFM to PCM conversion.
        num_columns (int): The expected number of columns in the input matrix.
        name (str): The name of the motif.
        matrix (list of lists): The matrix data.
        format (str): The format of the input matrix ('pcm', 'pfm', 'pwm').

    Methods:
        to_pwm():
            Converts the input matrix to PWM format and returns it as a dictionary.
    """

    def __init__(self, filename, pseudocount='log', word_count=1000.0, num_columns=4):
        """
        Initializes the PwmConverter with the given filename and conversion parameters.

        Parameters:
            filename (str): The input file name.
            pseudocount (str or float, optional): The pseudocount method ('log', 'sqrt') or a numeric value.
                                                  Defaults to 'log'.
            word_count (float, optional): The word count used for PFM to PCM conversion.
                                          Defaults to 1000.0.
            num_columns (int, optional): The expected number of columns in the input matrix.
                                         Defaults to 4.

        Raises:
            ValueError: If the input file is not found or the matrix is improperly formatted.
        """
        self.filename = filename
        self.pseudocount = pseudocount
        self.word_count = word_count
        self.num_columns = num_columns
        self.name = ''
        self.matrix = []
        self.format = ''

        self._read_matrix()
        self.format = self._determine_motif_format()

    def _is_float(self, s):
        """
        Checks if a string can be converted to a float.

        Parameters:
            s (str): The string to check.

        Returns:
            bool: True if convertible to float, False otherwise.
        """
        try:
            float(s)
            return True
        except ValueError:
            return False

    def _determine_motif_format(self):
        """
        Determines the motif format based on the file extension.

        Returns:
            str: The motif format ('pcm', 'pfm', or 'pwm').

        Raises:
            ValueError: If the file extension is unrecognized.
        """
        ext = os.path.splitext(self.filename)[1].lower()
        if ext == '.pcm':
            return 'pcm'
        elif ext in ['.pfm', '.ppm']:
            return 'pfm'
        elif ext == '.pwm':
            return 'pwm'
        else:
            raise ValueError(f"Unknown file extension '{ext}' for motif format.")

    def _read_matrix(self):
        """
        Reads the matrix from a file and populates the name and matrix attributes.

        Raises:
            ValueError: If the file is not found, the input is empty, or the matrix is improperly formatted.
        """
        try:
            with open(self.filename, 'r') as f:
                lines = f.read().splitlines()
        except FileNotFoundError:
            raise ValueError(f"File '{self.filename}' not found.")

        # Remove empty lines, strip whitespace, and remove comments
        lines = [line.strip() for line in lines if line.strip() and not line.strip().startswith('#')]

        if not lines:
            raise ValueError("Input is empty.")

        # Check if the first line is a header
        if lines[0].startswith('>'):
            header = lines.pop(0)
            self.name = header[1:].strip().split()[0]
        else:
            # Use the filename as the motif name
            self.name = os.path.splitext(os.path.basename(self.filename))[0]

        if not lines:
            raise ValueError("No data found after header.")

        # Now, split the lines into rows
        rows = [line.split() for line in lines if not line.startswith('#') and line]

        # Now check that each row has the expected number of columns
        for idx, row in enumerate(rows):
            if len(row) != self.num_columns:
                raise ValueError(f"All rows must have exactly {self.num_columns} columns. Error at row {idx + 1}.")
            if not all(self._is_float(x) for x in row):
                raise ValueError(f"Non-numeric value found in the matrix at row {idx + 1}.")

        # Convert all elements to floats
        self.matrix = [[float(x) for x in row] for row in rows]

    def _calculate_pseudocount(self, count):
        """
        Calculates the pseudocount based on the specified method.

        Parameters:
            count (float): The total count for the row.

        Returns:
            float: The calculated pseudocount.

        Raises:
            ValueError: If an invalid pseudocount method is provided.
        """
        if self.pseudocount == 'log':
            return math.log(max(count, 2))
        elif self.pseudocount == 'sqrt':
            return math.sqrt(count)
        elif isinstance(self.pseudocount, (int, float)):
            return self.pseudocount
        else:
            raise ValueError(f"Invalid pseudocount value: {self.pseudocount}")

    def _pcm2pfm(self, pcm):
        """
        Converts a PCM (Position Count Matrix) to PFM (Position Frequency Matrix).

        Parameters:
            pcm (dict): A dictionary with 'name' and 'matrix' keys.

        Returns:
            dict: A dictionary with 'name' and 'matrix' keys representing the PFM.

        Raises:
            ValueError: If a row sum is zero.
        """
        pfm_matrix = []
        for row in pcm['matrix']:
            total = sum(row)
            if total == 0:
                raise ValueError("Row sum is zero, cannot normalize.")
            pfm_matrix.append([x / total for x in row])
        return {'name': pcm['name'], 'matrix': pfm_matrix}

    def _pfm2pcm(self, pfm):
        """
        Converts a PFM (Position Frequency Matrix) to PCM (Position Count Matrix).

        Parameters:
            pfm (dict): A dictionary with 'name' and 'matrix' keys.

        Returns:
            dict: A dictionary with 'name' and 'matrix' keys representing the PCM.
        """
        pcm_matrix = []
        for row in pfm['matrix']:
            pcm_matrix.append([el * self.word_count for el in row])
        return {'name': pfm['name'], 'matrix': pcm_matrix}

    def _pcm2pwm(self, pcm):
        """
        Converts a PCM (Position Count Matrix) to PWM (Position Weight Matrix).

        Parameters:
            pcm (dict): A dictionary with 'name' and 'matrix' keys.

        Returns:
            dict: A dictionary with 'name' and 'matrix' keys representing the PWM.
        """
        pwm_matrix = []
        for row in pcm['matrix']:
            count = sum(row)
            pseudo = self._calculate_pseudocount(count)
            denominator = 0.25 * (count + pseudo)
            pwm_row = []
            for el in row:
                numerator = el + 0.25 * pseudo
                pwm_value = math.log(numerator / denominator)
                pwm_row.append(pwm_value)
            pwm_matrix.append(pwm_row)
        return {'name': pcm['name'], 'matrix': pwm_matrix}

    def to_pwm(self):
        """
        Converts the input matrix to PWM format.

        Returns:
            dict: A dictionary with 'name' and 'matrix' keys representing the PWM.

        Raises:
            ValueError: If the input format is unknown or unsupported for conversion.
        """
        if self.format == 'pwm':
            # Input is already PWM
            return {'name': self.name, 'matrix': self.matrix}
        elif self.format == 'pfm':
            # Convert PFM -> PCM -> PWM
            pfm = {'name': self.name, 'matrix': self.matrix}
            pcm = self._pfm2pcm(pfm)
            pwm = self._pcm2pwm(pcm)
            return pwm
        elif self.format == 'pcm':
            # Convert PCM -> PWM
            pcm = {'name': self.name, 'matrix': self.matrix}
            pwm = self._pcm2pwm(pcm)
            return pwm
        else:
            raise ValueError(f"Unknown motif format '{self.format}'.")

    def matrix_as_string(self, model, transpose_output=False):
        """
        Converts a matrix model to a string representation.

        Parameters:
            model (dict): A dictionary with 'name' and 'matrix' keys.
            transpose_output (bool, optional): Whether to transpose the matrix. Defaults to False.

        Returns:
            str: The string representation of the matrix.
        """
        name = model['name']
        matrix = model['matrix']
        if transpose_output:
            matrix = list(map(list, zip(*matrix)))
        lines = [f">{name}"]
        lines += ["\t".join(f"{x:.6f}" for x in row) for row in matrix]
        return "\n".join(lines)


class ThresholdPValueTable:
    def __init__(self):
        self._filename = "thresholdPValueTable"
        self._minScore = None
        self._maxScore = None
        self._scoreFactor = None
        self._granularity = None

    def init(self, min_score, max_score, factor, granularity):
        self._minScore = min_score
        self._maxScore = max_score
        self._scoreFactor = factor
        self._granularity = granularity

    def get_threshold_pvalue_table(self, matrix):
        if not matrix:
            raise ValueError("Matrix cannot be empty")

        e = self._granularity
        round_matrix = self.get_round_matrix(matrix, e)
        distribution = self.get_score_distribution(round_matrix)
        threshold_pvalue_table = [distribution.pop()]

        for i in range(len(distribution)):
            pair = distribution.pop()
            score = pair[0]
            probability = pair[1]
            previous_pvalue = threshold_pvalue_table[i][1]
            threshold_pvalue_table.append([score, previous_pvalue + probability])

        return threshold_pvalue_table[::-1]

    def get_score_distribution(self, matrix):
        if not matrix:
            raise ValueError("Matrix cannot be empty")

        previous_quantities = {0: 1}
        new_quantities = {}
        e = self._granularity
        total_combinations = 4 ** len(matrix)

        for position in range(len(matrix)):
            scores = list(previous_quantities.keys())
            new_quantities = {}

            for score in scores:
                previous_quantity = previous_quantities[score]

                for letter in range(4):
                    new_score = score + matrix[position][letter]
                    quantity = new_quantities.get(new_score, 0)
                    new_quantities[new_score] = quantity + previous_quantity

            previous_quantities = new_quantities

        sorted_scores = sorted(new_quantities.keys())
        score_distribution = [
            [score / 10 ** e, new_quantities[score] / total_combinations]
            for score in sorted_scores
        ]

        return score_distribution

    def get_round_matrix(self, matrix, granularity):
        if not isinstance(granularity, int) or granularity < 0:
            raise ValueError("Granularity must be a non-negative integer")
        if not matrix:
            raise ValueError("Matrix cannot be empty")
        round_matrix = []
        for row in matrix:
            round_row = []
            for value in row:
                if not isinstance(value, (int, float)):
                    raise ValueError(f"Matrix values must be numeric. Invalid value: {value}")
                round_row.append(int(value * 10 ** granularity))
            round_matrix.append(round_row)
        return round_matrix

    def get_best_score(self, matrix):
        if not matrix:
            raise ValueError("Matrix cannot be empty")
        return sum(max(row) for row in matrix)

    def get_worst_score(self, matrix):
        if not matrix:
            raise ValueError("Matrix cannot be empty")
        return sum(min(row) for row in matrix)


class Config:
    def __init__(self, config_file):
        self.config_file = config_file
        self.data = self.load_config()
        self.validate_config()

    def load_config(self):
        try:
            with open(self.config_file, 'r') as file:
                return json.load(file)
        except FileNotFoundError:
            raise ConfigError(f"Config file {self.config_file} not found.")
        except json.JSONDecodeError:
            raise ConfigError(f"Error decoding JSON from {self.config_file}.")

    def validate_config(self):
        for param in USER_INPUT | SOURCE_DATA | REQUIRED_PARAMETERS:
            if param not in self.data:
                raise ConfigError(f"Missing required configuration key: {param}")

        for param in self.data:
            if param not in ALL_PARAMETERS:
                raise ConfigError(f"Unknown parameter {param}")

    def get(self, key, default=None):
        return self.data.get(key, default)


class ArgParser:
    def __init__(self):
        self.parser = argparse.ArgumentParser()
        self.add_arguments()

    def add_arguments(self):
        self.parser.add_argument(
            '--config',
            type=str,
            default='./config.json',
            help='Path to the config JSON file')

        for param in USER_INPUT | SOURCE_DATA | REQUIRED_PARAMETERS | OPTIONAL_PARAMETERS:
            if param == 'custom_motifs':
                self.parser.add_argument(
                    '--custom_motifs',
                    action='store_true',
                    help='Use custom motif files and skip downloading data'
                )
            else:
                self.parser.add_argument(
                    f'--{param}',
                    type=str,
                    help=f'Value for {param.replace("_", " ")}'
                )

    def parse(self):
        args = self.parser.parse_args()
        self.validate_arguments(args)
        return args

    def validate_arguments(self, args):
        if args.config and not os.path.exists(args.config):
            raise ArgumentError(f"Config file {args.config} does not exist.")


class SettingsValidator:
    def __init__(self, settings: dict[str, any]):
        self.settings = settings

    def validate_settings(self):
        all_required_keys = USER_INPUT | SOURCE_DATA | REQUIRED_PARAMETERS

        # Step 1: Check for missing required parameters
        self._check_required_parameters(all_required_keys)

        # Step 2: Validate that required files exist
        self._validate_files_exist(USER_INPUT)

        # Step 3: Validate scoring_mode
        self._validate_scoring_mode()

        # Step 4: Validate and re-assign flanks
        self._validate_flanks()

    def get_flanks(self): # TODO Maybe transform into cast_settings later
        return {'flank_5': self.settings['flank_5'], 'flank_3': self.settings['flank_3']}

    def _check_required_parameters(self, required_keys: set[str]):
        missing_keys = [key for key in required_keys
                        if key not in self.settings.keys() or self.settings[key] is None]
        if missing_keys:
            missing = ', '.join(missing_keys)
            raise ArgumentError(f"Missing required parameter(s): {missing}")

    def _validate_files_exist(self, file_keys: set[str]):
        non_existent_files = [
            key for key in file_keys
            if key in self.settings and
               self.settings[key] != 'Unknown' and
               not os.path.exists(self.settings[key])
        ]
        if non_existent_files:
            errors = [
                f"No such file or directory for {key.replace('_', ' ')}: {self.settings[key]}"
                for key in non_existent_files
            ]
            error_message = "; ".join(errors)
            # raise ArgumentError(error_message) # TODO Create validation stage

    def _validate_scoring_mode(self):
        valid_modes = {'besthit', 'pfm-sum-occupancy', 'custom-thresholds'}
        scoring_mode = self.settings.get('scoring_mode')
        if scoring_mode not in valid_modes and not scoring_mode:
            raise ArgumentError(f"Invalid scoring mode: {scoring_mode}. "
                                f"Expected one mode of the folowing: {', '.join(valid_modes)}.")

    def _validate_flanks(self):
        flank_pairs_text = '\n'.join([pair for pair in AVAILABLE_FLANK_PAIRS])

        try:
            abs_flank_5 = int(self.settings['flank_5'][:-1])
            abs_flank_3 = int(self.settings['flank_3'][:-1])

            letter_flank_5 = self.settings['flank_5'][-1]
            letter_flank_3 = self.settings['flank_3'][-1]

            if not {letter_flank_5, letter_flank_5}.issubset({'u', 'd'}):
                raise ArgumentError(f'Flank pairs must be one of the following:\n{flank_pairs_text}')

            self.settings['flank_5'] = flank_5 = -1 * abs_flank_5 if letter_flank_5 == 'u' else abs_flank_5
            self.settings['flank_3'] = flank_3 = -1 * abs_flank_3 if letter_flank_3 == 'u' else abs_flank_3

        except (ValueError, IndexError):
            raise ArgumentError(f'Flank pairs must be one of the following:\n{flank_pairs_text}')


        # In case of reverting back intervals
        MIN_FLANK_5 = -250
        MAX_FLANK_5 = 10
        MIN_FLANK_3 = 0
        MAX_FLANK_3 = 50


class MaraPreprocessingApp:
    def __init__(self, config=None, args=None):
        self.config = config
        self.args = args
        self.settings = self.merge_settings()
        # self.settings = args # From main.nf
        self.validate_settings()
        self.setup_logging()
        # self.check_executables()
        self.define_directories()
        self.create_directories()
        self.define_output_prefixes()

    def merge_settings(self):
        settings = {}

        # Merge settings from config and command-line arguments
        for param in ALL_PARAMETERS:
            arg_value = getattr(self.args, param, None)
            config_value = self.config.get(param)
            value = arg_value if arg_value is not None else config_value
            if value is not None:
                settings[param] = value

        # Set defaults
        settings['log_level'] = settings.get('log_level', 'INFO')
        settings['num_processes'] = settings.get('num_processes', os.cpu_count())
        settings['scoring_mode'] = settings.get('scoring_mode', 'besthit').lower()
        settings['stage'] = settings.get('stage', 'all')
        settings['custom_motifs'] = settings.get('custom_motifs', False)

        if isinstance(settings['custom_motifs'], str):
            settings['custom_motifs'] = settings['custom_motifs'].lower() == 'true'

        # Prepend BASE_DIR to relevant paths
        for key in USER_INPUT.union(SOURCE_DATA):
            if key in settings:
                settings[key] = os.path.join(BASE_DIR, settings[key])

        return settings

    def validate_settings(self):
        validator = SettingsValidator(self.settings)
        validator.validate_settings()
        self.settings['flank_5'] = validator.settings['flank_5']
        self.settings['flank_3'] = validator.settings['flank_3']

    # def validate_settings(self):
    #     # Check required parameters
    #     for key in USER_INPUT | SOURCE_DATA | REQUIRED_PARAMETERS:
    #         if key not in self.settings or self.settings[key] is None:
    #             raise ArgumentError(f"Missing required parameter: {key}")

    #     # Validate files exist
    #     for key in USER_INPUT | SOURCE_DATA:
    #         if key in self.settings and self.settings[key] != 'Unknown' and not os.path.exists(self.settings[key]):
    #             raise ArgumentError(
    #                 f"No such file or directory for {key.replace('_', ' ')}: {self.settings[key]}"
    #             )

    #     # Validate scoring_mode

    #     # Validate and re-assign flanks
    #     # Attempt to convert flank_5 to integer
    #     try:
    #         self.settings['flank_5'] = int(self.settings.flank_5)
    #     except ValueError as e:
    #         raise ArgumentError("flank_5 must be a valid integer.") from e

    #     # Attempt to convert flank_3 to integer
    #     try:
    #         self.settings['flank_3'] = int(self.settings.flank_3)
    #     except ValueError as e:
    #         raise ArgumentError("flank_3 must be a valid integer.") from e

    #     # Define minimum and maximum values for flank_5 and flank_3
    #     MIN_FLANK_5 = 0       # Example minimum value for flank_5
    #     MAX_FLANK_5 = 100     # Example maximum value for flank_5
    #     MIN_FLANK_3 = 0       # Example minimum value for flank_3
    #     MAX_FLANK_3 = 100     # Example maximum value for flank_3

    #     # Check if both flanks are negative
    #     if flank_5 < 0 and flank_3 < 0:
    #         raise ArgumentError("Expected an even number of integers for flank pairs.")

    #     # Validate flank_5 within its range
    #     if not (MIN_FLANK_5 <= flank_5 <= MAX_FLANK_5):
    #         raise ArgumentError(
    #             f"Invalid flank_5 value: {flank_5}. "
    #             f"Expected a value between {MIN_FLANK_5} and {MAX_FLANK_5}."
    #         )

    #     # Validate flank_3 within its range
    #     if not (MIN_FLANK_3 <= flank_3 <= MAX_FLANK_3):
    #         raise ArgumentError(
    #             f"Invalid flank_3 value: {flank_3}. "
    #             f"Expected a value between {MIN_FLANK_3} and {MAX_FLANK_3}."
    #         )

    #     # Assign the validated values back to settings
    #     self.settings['flank_5'] = flank_5
    #     self.settings['flank_3'] = flank_3

    def setup_logging(self):
        log_level = self.settings.get('log_level', 'INFO').upper()
        numeric_level = getattr(logging, log_level, logging.INFO)
        logging.basicConfig(level=numeric_level, format='%(asctime)s - %(levelname)s - %(message)s')

    def check_executables(self):
        required_executables = ['bedtools', 'wget', 'java']
        for exe in required_executables:
            if shutil.which(exe) is None:
                raise ExternalToolError(f"Required executable '{exe}' not found in PATH.")

    def define_directories(self):
        self.base_dir = os.getcwd()
        self.stage_dirs = {
            'stage_00': os.path.join(BASE_DIR, 'source_data'),
            'stage_01': os.path.join(BASE_DIR, 'stages', 'stage_01'),
            'stage_02': os.path.join(BASE_DIR, 'stages', 'stage_02'),
            'stage_03': os.path.join(BASE_DIR, 'stages', 'stage_03'),
            'stage_04': os.path.join(BASE_DIR, 'stages', 'stage_04')
        }

    def create_directories(self):
        for dir_path in self.stage_dirs.values():
            os.makedirs(dir_path, exist_ok=True)

    def define_output_prefixes(self):
        self.output_prefixes = {
            'tss_clusters': os.path.join(self.stage_dirs['stage_02'], 'TSS_clusters'),
            'promoters': os.path.join(self.stage_dirs['stage_02'], 'hg38_promoters_v1')
        }

    def run(self):
        logging.info("Starting MARA Preprocessing Application")
        logging.info("Settings:")

        for key, value in self.settings.items():
            logging.info(f"  {key}: {value}")

        stage = self.settings['stage']
        custom_motifs = self.settings.get('custom_motifs', False)

        if stage == 'download_data' or stage == 'all':
            if not custom_motifs:
                self.download_data()
                # self.prepare_motif_files() # TODO maybe remove motif preparation after data downloading
            else:
                logging.info("Custom motifs flag is set. Skipping data download.")
                # Ensure motif files are prepared if not downloading data
                self.prepare_motif_files()

        if stage == 'preprocess' or stage == 'all':
            self.preprocess()

        if stage == 'flanking_regions' or stage == 'all':
            self.flanking_regions()

        if stage == 'extract_fasta' or stage == 'all':
            self.extract_fasta()

        if stage == 'motif_analysis' or stage == 'all':
            self.motif_analysis()

        logging.info("MARA Preprocessing Application finished successfully")

    def prepare_motif_files(self):
        """Prepare motif files by generating missing PWM and thresholds if necessary."""
        logging.info("Preparing motif files")

        motif_dir = self.settings['motif_dir']
        pwm_dir = os.path.join(motif_dir, 'pwm')
        pcm_dir = os.path.join(motif_dir, 'pcm')
        pfm_dir = os.path.join(motif_dir, 'pfm')
        thresholds_dir = os.path.join(motif_dir, 'thresholds')

        # Ensure directories exist
        os.makedirs(pwm_dir, exist_ok=True)
        os.makedirs(pcm_dir, exist_ok=True)
        os.makedirs(pfm_dir, exist_ok=True)
        os.makedirs(thresholds_dir, exist_ok=True)

        # Collect motif names from pcm_dir and pfm_dir
        motif_names = set()
        if os.path.exists(pcm_dir):
            pcm_files = [f for f in os.listdir(pcm_dir) if f.endswith('.pcm')]
            motif_names.update([os.path.splitext(f)[0] for f in pcm_files])
        if os.path.exists(pfm_dir):
            pfm_files = [f for f in os.listdir(pfm_dir) if f.endswith('.pfm')]
            motif_names.update([os.path.splitext(f)[0] for f in pfm_files])
        if os.path.exists(pwm_dir):
            pwm_files = [f for f in os.listdir(pwm_dir) if f.endswith('.pwm')]
            motif_names.update([os.path.splitext(f)[0] for f in pwm_files])

        if not motif_names:
            logging.error("No motif files found in motif_dir")
            raise FileNotFoundError("No motif files found in motif_dir")

        # For each motif, check if PWM and threshold files exist
        for motif_name in motif_names:
            pwm_file = os.path.join(pwm_dir, f"{motif_name}.pwm")
            if not os.path.exists(pwm_file):
                # Try to generate PWM file from PCM or PFM
                logging.info(f"Generating PWM for motif {motif_name}")
                pcm_file = os.path.join(pcm_dir, f"{motif_name}.pcm")
                pfm_file = os.path.join(pfm_dir, f"{motif_name}.pfm")
                if os.path.exists(pcm_file):
                    # Use PCM file to generate PWM
                    converter = PwmConverter(pcm_file)
                    pwm_model = converter.to_pwm()
                    pwm_content = converter.matrix_as_string(pwm_model)
                    # Save PWM file
                    with open(pwm_file, 'w') as f:
                        f.write(pwm_content)
                elif os.path.exists(pfm_file):
                    # Use PFM file to generate PWM
                    converter = PwmConverter(pfm_file)
                    pwm_model = converter.to_pwm()
                    pwm_content = converter.matrix_as_string(pwm_model)
                    # Save PWM file
                    with open(pwm_file, 'w') as f:
                        f.write(pwm_content)
                else:
                    logging.error(f"No PCM or PFM file found for motif {motif_name}")
                    continue
            else:
                logging.info(f"PWM file for motif {motif_name} already exists")

            # Now check if threshold file exists
            threshold_file = os.path.join(thresholds_dir, f"{motif_name}.threshold")
            if not os.path.exists(threshold_file):
                # Generate threshold file from PWM
                logging.info(f"Generating threshold for motif {motif_name}")
                if os.path.exists(pwm_file):
                    # Read PWM file
                    converter = PwmConverter(pwm_file)
                    pwm_model = converter.to_pwm()
                    pwm_matrix = pwm_model['matrix']
                    # Generate threshold p-value table
                    threshold_table = ThresholdPValueTable()
                    # Initialize threshold_table with parameters
                    min_score = threshold_table.get_worst_score(pwm_matrix)
                    max_score = threshold_table.get_best_score(pwm_matrix)
                    score_factor = 1.0
                    granularity = 2  # Example value; adjust as needed
                    threshold_table.init(min_score, max_score, score_factor, granularity)
                    # Generate threshold p-value table
                    threshold_pvalue_table = threshold_table.get_threshold_pvalue_table(pwm_matrix)
                    # Save threshold table to file
                    with open(threshold_file, 'w') as f:
                        for score, pvalue in threshold_pvalue_table:
                            f.write(f"{score}\t{pvalue}\n")
                else:
                    logging.error(f"PWM file not found for motif {motif_name}")
                    continue
            else:
                logging.info(f"Threshold file for motif {motif_name} already exists")
        logging.info("Motif files preparation completed")

    def download_data(self):
        """Stage 00: Download required files."""
        logging.info("Starting Stage 00: Download required files")

        genome_file = self.settings['genome_file']
        motif_dir = self.settings['motif_dir']
        genome_address = self.settings.get('genome_url', 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz')
        hocomoco_invivo = self.settings.get('hocomoco_url', 'https://hocomoco12.autosome.org/final_bundle/hocomoco12/H12INVIVO')
        file_names = [
            "H12INVIVO_annotation.jsonl",
            "H12INVIVO_pwm.tar.gz",
            "H12INVIVO_pcm.tar.gz",
            "H12INVIVO_pfm.tar.gz",
            "H12INVIVO_thresholds.tar.gz"
        ]

        # Ensure genome directory exists
        genome_dir = os.path.dirname(genome_file)
        os.makedirs(genome_dir, exist_ok=True)

        # Download and extract genome file
        if not os.path.exists(genome_file):
            try:
                self.run_command(f'wget -P {genome_dir} {genome_address}')
                self.run_command(f'gunzip {genome_file}.gz')
                logging.info(f"Genome file downloaded and extracted to {genome_file}")
                # Remove the genome archive after extraction
                genome_archive = f"{genome_file}.gz"
                if os.path.exists(genome_archive):
                    os.remove(genome_archive)
                    logging.info(f"Removed genome archive: {genome_archive}")
            except CommandExecutionError as e:
                logging.error(f"Error downloading genome file: {e}")
                raise e
        else:
            logging.info(f"Genome file {genome_file} already exists. Skipping download.")

        # Ensure motif directory exists
        os.makedirs(motif_dir, exist_ok=True)

        # Download and extract motif files
        for file_name in file_names:
            file_path = os.path.join(motif_dir, file_name)

            if file_name.endswith('.tar.gz'):
                unpacked_name = file_name[10:-7]  # removes 'H12INVIVO_' and '.tar.gz'
                unpacked_path = os.path.join(motif_dir, unpacked_name)

                if not os.path.exists(unpacked_path):
                    url = f"{hocomoco_invivo}/{file_name}"
                    try:
                        self.run_command(f"wget -P {motif_dir} {url}")
                        self.run_command(f"tar -zxf {file_path} -C {motif_dir}")
                        logging.info(f"Downloaded and extracted {file_name} to {unpacked_path}")
                        # Remove the motif archive after extraction
                        if os.path.exists(file_path):
                            os.remove(file_path)
                            logging.info(f"Removed motif archive: {file_path}")
                    except CommandExecutionError as e:
                        logging.error(f"Error downloading or extracting motif file {file_name}: {e}")
                        raise e
                else:
                    logging.info(f"Directory {unpacked_path} already exists. Skipping download and extraction.")
            else:
                if not os.path.exists(file_path):
                    url = f"{hocomoco_invivo}/{file_name}"
                    try:
                        self.run_command(f"wget -P {motif_dir} {url}")
                        logging.info(f"Downloaded {file_name} to {motif_dir}")
                    except CommandExecutionError as e:
                        logging.error(f"Error downloading motif file {file_name}: {e}")
                        raise e
                else:
                    logging.info(f"{file_path} already exists. Skipping download.")

        logging.info("Stage 00 completed successfully")

    def preprocess(self):
        """Stage 01: Preprocess input BED files."""
        logging.info("Starting Stage 01: Preprocess input BED files")

        input_bed = self.settings['tss_clusters_file']
        output_dir = self.stage_dirs['stage_01']
        output_bed = os.path.join(output_dir, 'tss_clusters_preprocessed.bed')

        if os.path.exists(output_bed):
            logging.info(f"Output file {output_bed} already exists. Skipping preprocessing.")
            return

        # Read input BED file and process
        try:
            with open(input_bed, 'r') as infile, open(output_bed, 'w') as outfile:
                for line in infile:
                    row = line.strip().split('\t')
                    chrom = row[0]
                    chrom_start = row[1]
                    chrom_end = row[2]
                    score = row[4]
                    strand = row[5]
                    name = f"hg38_v1_{chrom}_{strand}_{int(chrom_start) + 1}_{chrom_end}"
                    outfile.write(f"{chrom}\t{chrom_start}\t{chrom_end}\t{name}\t{score}\t{strand}\n")
            logging.info(f"Preprocessed promoters file written to {output_bed}")
        except Exception as e:
            logging.error(f"Error during preprocessing: {e}")
            raise e

        logging.info("Stage 01 completed successfully")

    def flanking_regions(self):
        """Stage 02: Generate flanking regions for the input BED file."""
        logging.info("Starting Stage 02: Generate flanking regions")

        tss_clusters_file = os.path.join(self.stage_dirs['stage_01'], 'tss_clusters_preprocessed.bed')
        flank_5 = self.settings.get('flank_5')
        flank_3 = self.settings.get('flank_3')
        output_dir = self.stage_dirs['stage_02']
        output_prefix = os.path.join(output_dir, 'flanks')

        if not tss_clusters_file or flank_5 is None or flank_3 is None:
            raise ArgumentError("tss_clusters_file, flank_5, and flank_3 must be provided for flanking_regions stage")

        if not os.path.exists(tss_clusters_file):
            raise ArgumentError(f"TSS clusters file {tss_clusters_file} does not exist. Please run preprocess stage first.")

        self.make_flanks(flank_5, flank_3, tss_clusters_file, output_prefix)

        logging.info("Stage 02 completed successfully")

    def make_flanks(self, flank_5, flank_3, input_file, output_prefix):
        """Utility function to create flanking regions."""
        output_file = f"{output_prefix}.bed"

        if os.path.exists(output_file):
            logging.info(f"Flanking regions file {output_file} already exists. Skipping.")
            return

        try:
            with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
                for line in self.process_bed_file(flank_5, flank_3, infile):
                    outfile.write(line + '\n')
            logging.info(f"Generated flanking regions file {output_file}")
        except Exception as e:
            logging.error(f"Error processing BED file {input_file}: {e}")
            raise e

    def process_bed_file(self, flank_5, flank_3, infile):
        """Processes a BED file based on flanking distances."""
        lines = infile.readlines()
        regions = [self.parse_line(line) for line in lines if not line.startswith('#')]

        for region in regions:
            chrom = region['chrom']
            strand = region['strand']
            name = region['name']
            position = region['chrom_start'] if strand == '+' else region['chrom_end']

            if strand == '+':
                chrom_start = position + flank_5
                chrom_end = position + flank_3
            elif strand == '-':
                chrom_start = position - flank_3
                chrom_end = position - flank_5
            else:
                raise ValueError(f"Unknown strand {strand}")

            # Ensure chrom_start is not negative
            chrom_start = max(chrom_start, 0)

            info = [
                chrom,
                int(chrom_start),
                int(chrom_end),
                name,
                '.',
                strand
            ]
            yield "\t".join(map(str, info))

    def extract_fasta(self):
        """Stage 03: Extract FASTA sequences for the flanks BED file."""
        logging.info("Starting Stage 03: Extract FASTA sequences")

        flanks_bed = os.path.join(self.stage_dirs['stage_02'], 'flanks.bed')
        genome_file = self.settings['genome_file']
        output_dir = self.stage_dirs['stage_03']
        output_fasta = os.path.join(output_dir, 'sequences.fa')

        if not os.path.exists(flanks_bed):
            raise ArgumentError(f"Flanks BED file {flanks_bed} does not exist. Please run flanking_regions stage first.")

        if os.path.exists(output_fasta):
            logging.info(f"FASTA file {output_fasta} already exists. Skipping.")
            return

        try:
            self.run_command(f'bedtools getfasta -bed {flanks_bed} -fi {genome_file} -name+ -s > {output_fasta}')
            logging.info(f"Extracted FASTA sequences to {output_fasta}")
            # self.run_command(f'bedtools getfasta -bed {flanks_bed} -fi {genome_file} -name+ -s')
            # raise CommandExecutionError
        except CommandExecutionError as e:
            logging.error(f"Error extracting FASTA sequences: {e}")
            raise e

        logging.info("Stage 03 completed successfully")

    def motif_analysis(self):
        """Stage 04: Perform motif occupancy analysis on the extracted FASTA sequences."""
        logging.info("Starting Stage 04: Perform motif occupancy analysis")

        scoring_mode = self.settings['scoring_mode']
        fasta_file = os.path.join(self.stage_dirs['stage_03'], 'sequences.fa')
        sarus_jar = 'src/sarus-2.1.0.jar'
        motifs_dir = self.settings['motif_dir']
        motif = self.settings.get('motif')
        output_dir = self.stage_dirs['stage_04']
        
        # if motif:
        #     logging.info(f'Processing motif: {motif}')
        #     pwm_file = os.path.join(motifs_dir, 'pwm', f'{motif}.pwm')
        #     threshold_file = os.path.join(motifs_dir, 'thresholds', f'{motif}.thr')
        #     sarus_out = os.path.join(output_dir, f'{motif}_occupancies')
        # else:
        #     logging.info("Processing all motifs")
        #     pwm_file = os.path.join(motifs_dir, 'pwm')
        #     threshold_file = os.path.join(motifs_dir, 'thresholds')
        #     sarus_out = os.path.join(output_dir, 'motif_occupancies')

        logging.info(f'Processing motif: {motif}')
        pwm_file = os.path.join(motifs_dir, 'pwm', f'{motif}.pwm')
        threshold_file = os.path.join(motifs_dir, 'thresholds', f'{motif}.thr')

        output_prefix = os.path.join(output_dir, f'{motif}_occupancies')
        sarus_out = f'{output_prefix}.fasta'
        output_tsv = f"{output_prefix}.tsv"

        if not os.path.exists(fasta_file):
            raise ArgumentError(f"FASTA file {fasta_file} does not exist. Please run extract_fasta stage first.")

        if os.path.exists(output_tsv):
            logging.info(f"Motif occupancy file {output_tsv} already exists. Skipping.")
            return

        try:
            # if motif:
            #     scoring_mode = f'{motifs_dir}/thresholds/{motif}.thr' if scoring_mode == 'custom-thresholds' else scoring_mode
            #     self.run_command(f'java -cp {sarus_jar} ru.autosome.SARUS {fasta_file} {motifs_dir}/pwm/{motif}.pwm {scoring_mode} > {sarus_out}')
            #     self.convert_sarus_scores_to_tsv(sarus_out, output_tsv)
            #     logging.info(f"Motif occupancy analysis completed for motif: {motif}")
            # else:
            #     scoring_mode = f'{motifs_dir}/thresholds' if scoring_mode == 'custom-thresholds' else scoring_mode
            #     self.run_command(f'java -cp {sarus_jar} ru.autosome.SARUS {fasta_file} {motifs_dir}/pwm {scoring_mode} > {sarus_out}')
            #     self.convert_sarus_scores_to_tsv(sarus_out, output_tsv)
            #     logging.info(f"Motif occupancy analysis completed for {fasta_file}")

            scoring_mode = f'{motifs_dir}/thresholds/{motif}.thr' if scoring_mode == 'custom-thresholds' else scoring_mode
            self.run_command(f'java -cp {sarus_jar} ru.autosome.SARUS {fasta_file} {motifs_dir}/pwm/{motif}.pwm {scoring_mode} > {sarus_out}')
            self.convert_sarus_scores_to_tsv(sarus_out, output_tsv)
            logging.info(f"Motif occupancy analysis completed for motif: {motif}")

        except CommandExecutionError as e:
            logging.error(f"Error during motif occupancy analysis for {fasta_file}: {e}")
            raise e

        logging.info("Stage 04 completed successfully")

    def convert_sarus_scores_to_tsv(self, scores_fasta_file, output_tsv_file):
        """Converts SARUS scores from FASTA format to TSV format."""
        # Define the regular expression pattern
        pattern = re.compile(r"^>(?P<name>.+)::(?P<chr>chr\w+):(?P<from>\d+)-(?P<to>\d+)\((?P<strand>[+-])\)$")

        try:
            with open(scores_fasta_file, 'r') as infile, open(output_tsv_file, 'w') as outfile:
                lines = infile.readlines()

                for i in range(0, len(lines), 2):
                    hdr = lines[i].strip()
                    score = lines[i + 1].strip()

                    m = pattern.match(hdr)
                    if m:
                        chr, from_, to, name, strand = m.group('chr'), m.group('from'), m.group('to'), m.group('name'), m.group('strand')
                        score = score.split("\t")[0]  # Get the first score
                        info = [chr, from_, to, name, score, strand]
                        outfile.write("\t".join(info) + "\n")
                    else:
                        logging.warning(f"Header format does not match expected pattern: {hdr}")
            logging.info(f"Converted SARUS scores to TSV format: {output_tsv_file}")
        except Exception as e:
            logging.error(f"Error converting SARUS scores to TSV: {e}")
            raise e

    @staticmethod
    def run_command(command):
        """Utility function to run shell commands."""
        logging.debug(f"Running command: {command}")
        try:
            subprocess.run(command, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            raise CommandExecutionError(f"Command failed: {command}\n{e}")

    @staticmethod
    def parse_line(line):
        """Parses a line from a BED file into a dictionary."""
        row = line.strip().split("\t")
        return {
            "chrom": row[0],
            "chrom_start": int(row[1]),
            "chrom_end": int(row[2]),
            "name": row[3],
            "strand": row[5],
        }

def main():
    arg_parser = ArgParser()
    args = arg_parser.parse()
    config = Config(args.config)

    try:
        app = MaraPreprocessingApp(config=config, args=args)
        app.run()
    except (ConfigError, ArgumentError, ExternalToolError, CommandExecutionError) as e:
        logging.error(e)
        sys.exit(1)

if __name__ == "__main__":
    main()


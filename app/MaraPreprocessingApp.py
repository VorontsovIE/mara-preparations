#!/usr/bin/env python

# All available parameters
USER_INPUT = {
    'tss_clusters_file',
    'bed_chunk',
    'fasta_chunk'
}

REQUIRED_DATA = {
    'genome_file',
    'motif_dir'
}

REQUIRED_PARAMETERS = {
    'scoring_mode',
    'flank_pairs'
}

OPTIONAL_PARAMETERS = {
    'num_processes',
    'log_level',
    'stage',
    'chunk_id',
    'flank_5',
    'flank_3'
}

ALL_PARAMETERS = USER_INPUT | REQUIRED_DATA | REQUIRED_PARAMETERS | OPTIONAL_PARAMETERS

import os
import sys
import json
import re
import subprocess
import argparse
import shutil
import logging

class ConfigError(Exception):
    pass


class ArgumentError(Exception):
    pass


class ExternalToolError(Exception):
    pass


class CommandExecutionError(Exception):
    pass


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
        for param in USER_INPUT | REQUIRED_DATA | REQUIRED_PARAMETERS:
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

        for param in USER_INPUT | REQUIRED_DATA | REQUIRED_PARAMETERS | OPTIONAL_PARAMETERS:
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

class MaraPreprocessingApp:
    def __init__(self, config=None, args=None):
        self.config = config
        self.args = args
        self.settings = self.merge_settings()
        self.validate_settings()
        self.setup_logging()
        self.check_executables()
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
        settings['scoring_mode'] = settings.get('scoring_mode', 'upstream').lower()
        settings['stage'] = settings.get('stage', 'all')

        return settings

    def validate_settings(self):
        # Check required parameters
        for key in USER_INPUT | REQUIRED_DATA | REQUIRED_PARAMETERS:
            if key not in self.settings or self.settings[key] is None:
                raise ArgumentError(f"Missing required parameter: {key}")

        # Validate files exist
        for key in USER_INPUT | REQUIRED_DATA:
            if key in self.settings and self.settings[key] != 'Unknown' and not os.path.exists(self.settings[key]):
                raise ArgumentError(
                    f"No such file or directory for {key.replace('_', ' ')}: {self.settings[key]}"
                )

        # Validate scoring_mode
        if self.settings['scoring_mode'] not in ['upstream', 'downstream']:
            raise ArgumentError(
                f"Invalid scoring mode: {self.settings['scoring_mode']}"
            )

        # Parse flank_pairs
        try:
            pairs = self.settings['flank_pairs'].split()

            if len(pairs) % 2 != 0:
                raise ArgumentError("Expected an even number of integers for flank pairs.")

            self.settings['flank_pairs'] = [(int(pairs[i]), int(pairs[i + 1])) for i in range(0, len(pairs), 2)]

        except ValueError as e:
            raise ArgumentError(
                f"Invalid flank pairs format: {str(e)} Each pair must consist of two \
integers in the format ‘x y’ (e.g., ‘1 2 3 4’). Please ensure \
that the input is correctly formatted with pairs of \
space-separated numbers."
            )

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
            'stage_00': os.path.join(self.base_dir, 'stages', 'stage_00'),
            'stage_01': os.path.join(self.base_dir, 'stages', 'stage_01'),
            'stage_02': os.path.join(self.base_dir, 'stages', 'stage_02'),
            'stage_03': os.path.join(self.base_dir, 'stages', 'stage_03'),
            'stage_04': os.path.join(self.base_dir, 'stages', 'stage_04')
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

        if stage == 'download_data' or stage == 'all':
            self.download_data()

        if stage == 'preprocess' or stage == 'all':
            self.preprocess()

        if stage == 'flanking_regions' or stage == 'all':
            self.flanking_regions()

        if stage == 'extract_fasta' or stage == 'all':
            self.extract_fasta()

        if stage == 'motif_analysis' or stage == 'all':
            self.motif_analysis()

        logging.info("MARA Preprocessing Application finished successfully")

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

                if not os.path.exists(file_path) and not os.path.exists(unpacked_path):
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
                elif os.path.exists(unpacked_path):
                    logging.info(f"Directory {unpacked_path} already exists. Skipping download and extraction.")
                else:
                    logging.info(f"Archive {file_path} already exists. Skipping download.")
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
        """Stage 02: Generate flanking regions for a BED file chunk."""
        logging.info("Starting Stage 02: Generate flanking regions")

        bed_chunk = self.settings['bed_chunk']
        chunk_id = self.settings['chunk_id']
        flank_5 = self.settings.get('flank_5')
        flank_3 = self.settings.get('flank_3')
        output_dir = self.stage_dirs['stage_02']
        output_prefix = os.path.join(output_dir, f'chunk_{chunk_id}')
        scoring_mode = self.settings['scoring_mode']

        if not bed_chunk or flank_5 is None or flank_3 is None or chunk_id is None:
            raise ArgumentError("bed_chunk, flank_5, flank_3, and chunk_id must be provided for flanking_regions stage")


        self.make_flanks(flank_5, flank_3, bed_chunk, output_prefix, scoring_mode)

        logging.info(f"Stage 02 completed successfully for chunk {chunk_id}")

    def make_flanks(self, flank_5, flank_3, input_file, output_prefix, scoring_mode):
        """Utility function to create flanking regions."""
        output_file = f"{output_prefix}_flanks.bed"

        if os.path.exists(output_file):
            logging.info(f"Flanking regions file {output_file} already exists. Skipping.")
            return

        try:
            with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
                for line in self.process_bed_file(flank_5, flank_3, infile, scoring_mode):
                    outfile.write(line + '\n')
            logging.info(f"Generated flanking regions file {output_file}")
        except Exception as e:
            logging.error(f"Error processing BED file {input_file}: {e}")
            raise e

    def process_bed_file(self, flank_5, flank_3, infile, scoring_mode):
        """Processes a BED file based on flanking distances and scoring mode."""
        lines = infile.readlines()
        regions = [self.parse_line(line) for line in lines if not line.startswith('#')]

        for region in regions:
            chrom = region['chrom']
            strand = region['strand']
            name = region['name']
            position = region['chrom_start'] if strand == '+' else region['chrom_end']

            if scoring_mode == 'upstream':
                if strand == '+':
                    chrom_start = position - flank_5
                    chrom_end = position - flank_3
                elif strand == '-':
                    chrom_start = position + flank_3
                    chrom_end = position + flank_5
                else:
                    raise ValueError(f"Unknown strand {strand}")
            elif scoring_mode == 'downstream':
                if strand == '+':
                    chrom_start = position + flank_3
                    chrom_end = position + flank_5
                elif strand == '-':
                    chrom_start = position - flank_5
                    chrom_end = position - flank_3
                else:
                    raise ValueError(f"Unknown strand {strand}")
            else:
                raise ValueError(f"Unknown scoring mode {scoring_mode}")

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
        """Stage 03: Extract FASTA sequences for a BED file chunk."""
        logging.info("Starting Stage 03: Extract FASTA sequences")

        bed_chunk = self.settings['bed_chunk']
        chunk_id = self.settings['chunk_id']
        genome_file = self.settings['genome_file']
        output_dir = self.stage_dirs['stage_03']
        output_fasta = os.path.join(output_dir, f'chunk_{chunk_id}_sequences.fa')

        if not bed_chunk or chunk_id is None:
            raise ArgumentError("bed_chunk and chunk_id must be provided for extract_fasta stage")

        if os.path.exists(output_fasta):
            logging.info(f"FASTA file {output_fasta} already exists. Skipping.")
            return

        try:
            self.run_command(f'bedtools getfasta -bed {bed_chunk} -fi {genome_file} -name+ -s > {output_fasta}')
            logging.info(f"Extracted FASTA sequences to {output_fasta}")
        except CommandExecutionError as e:
            logging.error(f"Error extracting FASTA sequences: {e}")
            raise e

        logging.info(f"Stage 03 completed successfully for chunk {chunk_id}")

    def motif_analysis(self):
        """Stage 04: Perform motif occupancy analysis on a FASTA file chunk."""
        logging.info("Starting Stage 04: Perform motif occupancy analysis")

        fasta_chunk = self.settings['fasta_chunk']
        chunk_id = self.settings['chunk_id']
        sarus_jar = 'app/sarus-2.1.0.jar'
        motifs_dir = self.settings['motif_dir']
        output_dir = self.stage_dirs['stage_04']
        sarus_out = os.path.join(output_dir, f'chunk_{chunk_id}_motif_occupancies')
        sarus_out_scores = f"{sarus_out}.scoresFasta"
        output_tsv = f"{sarus_out}.tsv"

        if not fasta_chunk or chunk_id is None:
            raise ArgumentError("fasta_chunk and chunk_id must be provided for motif_analysis stage")

        if os.path.exists(output_tsv):
            logging.info(f"Motif occupancy file {output_tsv} already exists. Skipping.")
            return

        try:
            # Run SARUS for motif occupancy analysis
            self.run_command(f'java -jar {sarus_jar} --pwm {motifs_dir}/pwm --threshold {motifs_dir}/thresholds --out {sarus_out} --fasta {fasta_chunk}')

            # Convert SARUS output to TSV format
            self.convert_sarus_scores_to_tsv(sarus_out_scores, output_tsv)

            logging.info(f"Motif occupancy analysis completed for {fasta_chunk}")

        except CommandExecutionError as e:
            logging.error(f"Error during motif occupancy analysis for {fasta_chunk}: {e}")
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
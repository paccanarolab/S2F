import configparser
import os

CONFIG = configparser.ConfigParser()
RUN_CONFIG = configparser.ConfigParser()


def load_run(run_config):
    RUN_CONFIG.read(run_config)
    config_file = os.path.expanduser(RUN_CONFIG.get('configuration', 'config_file'))
    CONFIG.read(config_file)


def save_run(run_config):
    RUN_CONFIG.write(run_config)


def load_configuration(config_file):
    CONFIG.read(config_file)


def save(config_file):
    CONFIG.write(config_file)



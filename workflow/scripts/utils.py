import logging
import argparse

def set_logging(log_path, log_level):
    logging.basicConfig(
        level=log_level, 
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler()
        ]
    )

def get_default_parser(description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--log', type=str, default="log.log", help='Log file')
    parser.add_argument('--log-level', type=str, default="INFO", help='Log level')
    return parser
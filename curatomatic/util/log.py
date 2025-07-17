import logging
from pathlib import Path

from curatomatic.model.log import LogLevel

logger = logging.getLogger('curatomatic')


class CentreFormatter(logging.Formatter):
    def format(self, record):
        if record.levelname == 'WARNING':
            record.levelname = ' WARN'
        elif record.levelname == 'CRITICAL':
            record.levelname = ' CRIT'
        elif record.levelname == 'INFO':
            record.levelname = ' INFO'
        return super().format(record)


def log_setup(out_dir: Path, level: LogLevel) -> logging.Logger:
    logger.setLevel(level.to_int())

    # create file handler which logs even debug messages
    fh = logging.FileHandler(out_dir / 'curatomatic.log')
    fh.setLevel(level.to_int())

    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(level.to_int())

    # create formatter and add it to the handlers
    formatter = CentreFormatter("[%(asctime)s] [%(levelname)s] - %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger

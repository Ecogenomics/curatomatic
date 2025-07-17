import logging
from enum import Enum


class LogLevel(str, Enum):
    DEBUG = 'debug'
    INFO = 'info'
    WARNING = 'warning'
    ERROR = 'error'
    CRITICAL = 'critical'

    def to_int(self):
        if self is LogLevel.DEBUG:
            return logging.DEBUG
        elif self is LogLevel.INFO:
            return logging.INFO
        elif self is LogLevel.WARNING:
            return logging.WARNING
        elif self is LogLevel.ERROR:
            return logging.ERROR
        elif self is LogLevel.CRITICAL:
            return logging.CRITICAL
        else:
            raise ValueError(f'Unknown log level: {self}')

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 10:31:19 2021

Simple logger.

@author: omarchehaimi
"""

__author__ = 'Omar Chehaimi'
__copyright__ = 'Prof. M.C. Righi, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 2nd, 2021'

import logging
import os
import json
from pathlib import Path, PurePosixPath

# Global variables 
levels = {'DEBUG': 10, 'INFO': 20, 'WARNING': 30, 
          'ERROR': 40, 'CRITICAL': 50}
log_format = '%(process)d - %(asctime)s - %(levelname)s - %(message)s'

class LoggingBase:
    """
    Base class which implements all the logging operations.
    The logging messages are coloured based as follow:
        DEBUG: green
        INFO: white on black screen and viceversa
        WARNING: yellow (default level of logging library)
        ERROR: red
        CRITICAL: bold red

    Attributes
    ----------
    name : str
        Name of the logger.
    console_level : int
        Console debug level.
    path : str (optional)
        Path to save the logfile.
    file_level : int
        Logfile debug level.
    log_format : str
        Format fot the debug messages.
    
    Methods
    -------
    __check_debug_level(level)
        Checks if the level is a valid debug level.

    __initialize_console_logger_handler(logger, console_level, log_format)
        Initializes a logger handler for the console.

    __initialize_file_logger_handler(logger, file_level, log_format)
        Initializes a logger handler for the file.

    __get_config
        Private method that gets the project configurations from a json file 
        saved in /core/config.json.

    initialize_logfolder()
        Initializes folder in which save logfiles. In case it already exists it 
        does not create anything. The log folder is one level out the main
        folder of the project (../log).

    debug(message)

    info(message)

    warnign(message)

    error(message)

    critical(message)

    get_config()
        Static method that loads the json file containing the configurations
        settings.

    logging_path()
        Static method that gets the path where to save the log files. In case 
        the path was not previously generated the function creates a new one. 
    
    """

    def __init__(self, name, console_level=30, path=None, 
                 file_level=30, log_format=log_format):
        """
        Parameters
        ----------
        name : str
            Name of the logger.
        console_level : int, optional
            Console debug level. The default is 30.
        path : str or None, optional
            Path to save the logfile. The default is None.
        file_level : int, optional
            Logfile debug level. The default is 30.
        log_format : str, optional
            Format fot the debug messages. The default is log_format (gloabl 
            variable).

        """

        # Check debug levels
        self.__check_debug_level(console_level)
        self.__check_debug_level(file_level)

        self.name = name
        self.console_level = console_level
        self.path = path
        self.file_level = file_level

        # Set the root level to DEBUG (the default value is WARNING)
        logging.root.setLevel(logging.DEBUG)

        # Initialize console logger        
        self.__initialize_console_logger_handler(
            name=self.name,
            console_level=self.console_level,
            log_format=log_format)
        # Create console logger
        self.console_logger = logging.getLogger(self.name+'_console_logger')
        
        if path:
            # Initialize file logger
            self.__initialize_file_logger_handler(
                name=self.name,
                file_level=file_level,
                path=path,
                log_format=log_format)
            # Create file logger
            self.file_logger = logging.getLogger(self.name+'_file_logger')

    def __check_debug_level(self, level):
        """
        Checks if the debug level is a valid one.
        
        Parameters
        ----------
        level : int
            Debug level.

        """

        if level not in levels.values():
            raise ValueError("'{}' level is not a valid level."
                             " Valid levels are: DEBUG: 10, INFO: 20, "
                             "WARNING: 30, ERROR: 40, "
                             "CRITICAL: 50.".format(level))

    def __initialize_console_logger_handler(self, name, console_level, 
                                            log_format):
        """
        Initializes a logger handler for the console.

        Parameters
        ----------
        name : str
            Name of the console logger.
        console_level : int
            Console level debug.
        log_format : str
            Log format.

        """

        console_logger = logging.getLogger(name + '_console_logger')

        console_handler = logging.StreamHandler()
        console_handler.setLevel(console_level)

        console_formatter = logging.Formatter(log_format)
        console_handler.setFormatter(console_formatter)

        console_logger.addHandler(console_handler)

    def __initialize_file_logger_handler(self, name, file_level, path, 
                                         log_format):
        """
        Initializes a logger handler for the file.

        Parameters
        ----------
        name : str
            Name of the file logger.
        file_level : int
            File level debug.
        path : str
            Path to save the log file.
        log_format : str
            Log format.

        """

        file_logger = logging.getLogger(name+'_file_logger')
        
        file_handler = logging.FileHandler(path)
        file_handler.setLevel(file_level)

        file_formatter = logging.Formatter(log_format)
        file_handler.setFormatter(file_formatter)

        file_logger.addHandler(file_handler)
    
    def debug(self, message):
        """ 
        Prints the debug level message.
        
        Parameters
        ----------
        message : str
            Message to print at debug level.

        """

        CGREEN = '\33[32m'
        CEND = '\033[0m'
        console_message = CGREEN + message + CEND
        self.console_logger.debug(console_message)

        if self.path:
            self.file_logger.debug(message)

    def info(self, message):
        """ 
        Prints the info level message.
        
        Parameters
        ----------
        message : str
            Message to print at debug level.

        """

        self.console_logger.info(message)

        if self.path:
            self.file_logger.info(message)

    def warning(self, message):
        """ 
        Prints the warning level message.
        
        Parameters
        ----------
        message : str
            Message to print at debug level.

        """

        CYELLOW = '\33[33m'
        CEND = '\033[0m'
        console_message = CYELLOW + message + CEND
        self.console_logger.warning(console_message)

        if self.path:
            self.file_logger.warning(message)

    def error(self, message):
        """ 
        Prints the error level message.
        
        Parameters
        ----------
        message : str
            Message to print at debug level.

        """

        CRED = '\33[31m'
        CEND = '\033[0m'
        console_message = CRED + message + CEND
        self.console_logger.error(console_message)

        if self.path:
            self.file_logger.error(message)

    def critical(self, message):
        """ 
        Prints the critical level message.
        
        Parameters
        ----------
        message : str
            Message to print at debug level.

        """

        CBOLD = '\033[1m'
        CRED = '\33[31m'
        CEND = '\033[0m'
        console_message = CBOLD + CRED + message + CEND
        self.console_logger.critical(console_message)

        if self.path:
            self.file_logger.critical(message)
    
    @staticmethod
    def get_config():
        """
        Static method that gets the project configurations from a json file 
        saved in /core/config.json.

        Return
        ------
        config : dict
            Dictionary containing the configurations.
        """

        project_folder = os.getcwd()
        config_path = project_folder + '/triboflow/core/config.json' 
        with open(config_path, 'r') as config:
            config = json.load(config)
        
        return config

    @staticmethod
    def logging_path():
        """
        Static method that gets the path where to save the logs file.
        The path is one level up the triboflow folder (triboflow../log) and it 
        should be named log. In case this folder does not exist a new one 
        is generated.

        Return
        ------
        log_path : str
            The path to log folder.
        """
        
        project_folder = os.getcwd()
        # PurePosixPath gets the first level parten directory
        log_folder_object = PurePosixPath(project_folder)
        log_folder = str(log_folder_object.parent) + '/log/'
        log_path = Path(log_folder)

        if not log_path.is_dir():
            print("WARNING: There is no folder for log files.")
            print("Creating a new log folder in " + log_folder)
            os.mkdir(log_folder)
            log_path = Path(log_folder)
            if not log_path.is_dir():
                raise RuntimeError('The creation of log path has failed!')

        return str(log_path)

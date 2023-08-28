#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 11:31:19 2021

Simple logger.

@author: omarchehaimi
"""

__author__ = "Omar Chehaimi"
__copyright__ = "Prof. M.C. Righi, University of Bologna"
__contact__ = "clelia.righi@unibo.it"
__date__ = "February 2nd, 2021"

from triboflow.core.logging import LoggingBase

configurations = LoggingBase.get_config()
log_path = LoggingBase.logging_path()

# General log for failures caused by the computer
log = LoggingBase(
    name="log",
    console_level=configurations["logging"]["console_logging_level"],
    path=log_path + "/log.log",
    file_level=configurations["logging"]["file_logging_level"],
)
log.debug("Test debug")
log.info("Test info")
log.warning("Test warning")
log.error("Test error")
log.critical("Test critical")

# Physics log
phys_log = LoggingBase(
    name="phys_log",
    console_level=configurations["logging"]["console_logging_level"],
    path=log_path + "/phys_log.log",
    file_level=configurations["logging"]["file_logging_level"],
)
phys_log.debug("Test debug phys")
phys_log.info("Test info phys")
phys_log.warning("Test warning phys")
phys_log.error("Test error phys")
phys_log.critical("Test critical phys")

# Writing general log messages
log.debug("Test debug after...")
log.info("Test info after...")
log.warning("Test warning after...")
log.error("Test error after...")
log.critical("Test critical after...")

# Another log for physics which write on the previous one.
# The log messages are not doubled because the name of the logger is different!
phys_log_test = LoggingBase(
    name="phys_log_test",
    console_level=configurations["logging"]["console_logging_level"],
    path=log_path + "/phys_log_test.log",
    file_level=configurations["logging"]["file_logging_level"],
)

phys_log_test.debug("Test debug phys phys_log_test")
phys_log_test.info("Test info phys phys_log_test")
phys_log_test.warning("Test warning phys phys_log_test")
phys_log_test.error("Test error phys phys_log_test")
phys_log_test.critical("Test critical phys phys_log_test")

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Local theory interfaces for the ASH workflow pipeline.

This module contains custom interface implementations for various quantum
chemistry and machine learning theories that are not yet part of the
upstream ASH library.
"""

from .interface_DPA3 import DPA3Theory

__all__ = ["DPA3Theory"]

#!/usr/bin/env python3
"""Aggregate collagen IV profiles by gland and plot the BB-94 titration."""
from pathlib import Path
from titration_analysis import cli
from titration_collagen_io import load_collagen_data

if __name__ == "__main__":
    cli("collagen", load_collagen_data, Path(__file__).with_name("collagen-inputs.json"))

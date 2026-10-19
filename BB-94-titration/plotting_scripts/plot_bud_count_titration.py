#!/usr/bin/env python3
"""Decode paired counts and plot the gland-level BB-94 titration response."""
from pathlib import Path
from titration_analysis import cli
from titration_bud_io import load_bud_data

if __name__ == "__main__":
    cli("bud", load_bud_data, Path(__file__).with_name("bud-inputs.json"))

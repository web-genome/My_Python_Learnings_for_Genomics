# fasta_picker.py
"""
A simple module to pick a FASTA file from your device using Tkinter.
Can be imported into any script and reused.
"""

import tkinter as tk
from tkinter import filedialog

def pick_fasta_file():
    """
    Opens a file picker dialog and returns the selected FASTA file path.
    Returns an empty string if the user cancels.
    """
    root = tk.Tk()
    root.withdraw()  # Hide the main Tkinter window

    file_path = filedialog.askopenfilename(
        title="Select a FASTA file",
        filetypes=[("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*")]
    )

    return file_path

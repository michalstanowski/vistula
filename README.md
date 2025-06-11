<p align="center">
    <img src="images/logo.png" alt="" width="50%"/>
</p>


# 🧬 VistulaSwitch

**VistulaSwitch** is a command-line tool for designing RNA toehold switches with thermodynamic evaluation, intended for applications in synthetic biology, diagnostics, and gene regulation.

---

## 🚀 Overview

VistulaSwitch automates the design of synthetic RNA toehold switches based on user-provided trigger and reporter sequences. These switches are capable of regulating translation in response to specific RNA inputs.

The tool integrates **RNAcofold** and **RNAup** from the ViennaRNA package to evaluate candidate constructs by analyzing:

- Free energy of complex formation
- Opening/accessibility energy
- Structural ensemble distribution

---

## 🧠 Features

- 🔬 **Toehold design**: Generates switch constructs tailored to input sequences  
- ⚡ **Thermodynamic analysis**: Uses ViennaRNA tools to evaluate folding and interaction energies  
- 🧪 **Composite scoring**: Ranks candidates based on a combination of energy and structure metrics  
- 🛠️ **CLI interface**: Designed for use in custom pipelines or batch analyses  

---

## 📦 Requirements

- Python ≥ 3.8  
- [ViennaRNA package](https://www.tbi.univie.ac.at/RNA/) (must include `RNAcofold` and `RNAup`)  
- Optional: Biopython, tqdm, matplotlib (if extended output is used)

Install required Python packages:
```bash
pip install -r requirements.txt

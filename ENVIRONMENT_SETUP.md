# TCCN - Environment Setup

This repository contains code for computing condition numbers in transcript space. To run the code successfully, you need to set up the proper Python environment with specific package versions.

## Requirements

- Python 3.10
- pandas 2.3.1
- numpy 2.2.5
- matplotlib 3.10.5
- brokenaxes 0.6.2
- sympy 1.14.0
- scipy 1.15.3
- seaborn 0.13.2
- scikit-learn 1.7.1

## Installation Options

### Option 1: Using conda (Recommended)

If you have conda installed, you can create the environment directly from the provided `environment.yml` file:

```bash
conda env create -f environment.yml
conda activate tccn
```

### Option 2: Using pip

If you prefer using pip, you can install the required packages using the `requirements.txt` file:

```bash
# Create a new virtual environment (optional but recommended)
python -m venv tccn_env
source tccn_env/bin/activate  # On Windows: tccn_env\Scripts\activate

# Install required packages
pip install -r requirements.txt
```

### Option 3: Manual Installation

You can also install the packages manually:

```bash
# Using conda
conda create -n tccn python=3.10
conda activate tccn
conda install pandas=2.3.1 numpy=2.2.5 matplotlib=3.10.5 scipy=1.15.3 seaborn=0.13.2 scikit-learn=1.7.1
pip install brokenaxes==0.6.2 sympy==1.14.0

# Or using pip
pip install pandas==2.3.1 numpy==2.2.5 matplotlib==3.10.5 brokenaxes==0.6.2 sympy==1.14.0 scipy==1.15.3 seaborn==0.13.2 scikit-learn==1.7.1
```

## Running the Code

After setting up the environment, you can run the examples:

1. Make sure your environment is activated
2. Open the Jupyter notebook: `jupyter notebook Example.ipynb`
3. Run the Python scripts: `python Process_gene_data.py`

## File Descriptions

- `Example.ipynb`: Jupyter notebook containing examples of computing condition numbers
- `Process_gene_data.py`: Main script for processing gene data
- `ComputeC.py`: Core computation functions
- `Gene_Data/`: Directory containing example gene data files
- `requirements.txt`: pip requirements file
- `environment.yml`: conda environment file

## Troubleshooting

If you encounter any issues:

1. Make sure you're using the exact package versions specified above
2. Verify that Python 3.10 is being used
3. Ensure all data files in the `Gene_Data/` directory are present
4. Check that your Jupyter notebook kernel is using the correct environment

For any questions or issues, please refer to the paper or contact the authors.

# gene_half_life
A Python code for estimating the half-life of decaying genes using decay time course datasets

Programming Language and Version: Python 3.9.13

Required packages/Libraries: 
1, Pandas
- To install Pandas package, type the following command in command prompt
	py -m pip install "pandas"
- To import the package to start using it, use the following command
	import pandas as pd

2, Numpy
- To install Numpy package, type the following command in command prompt
	py -m pip install "Numpy"
- To import the package to start using it, use the following command
	import Numpy as np

3, scikit_learn
- To install scikit-learn, type the following command in command prompt
	py -m pip install scikit-learn
- To import the Linear Regression model, use the following command
	from sklearn.linear_model import LinearRegression

Required input files: 
1, DecayTimecourse.txt

Description: 
-The code defines a function to calculate the half life of given genes by taking the log of gene decay over time and fitting it to the linear regression model. From the linear regression model, it takes the slope to calculate half life using formula t1/2= ln(2)/k

Execution: 
1,Open Command Prompt.
2, Change the working directory to the location of the file.
	cd path\to\the\file\location
3, Run the script using the following command.
	half_life_calculation.py
4, Once execution is complete check output.
	
Output Files:
None

Author: Amulya Saini
Date: 04/03/2023

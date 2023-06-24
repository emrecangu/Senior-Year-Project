# Senior-Year-Project
## This is the Graduation Project's ,Comparison of Rare Variant Association Tests in Obesity Cohort, main explanation.

This repository contains Python and R scripts for analyzing obesity-related data. The Python script searches for specific texts in Excel files and generates results, while the R script performs SKAT analysis and generates QQ plots.

### R Script
The R script (Senior-Year-Project/RCodes/MainCode.R) performs SKAT analysis and generates QQ plots. 

Follow the steps below to run the script:

Install the required packages by running install.packages("SKAT") and install.packages("fdrtool") in an R environment.

Open the script file and modify the following variables:

File.Bed: Path to the .bed file.

File.Bim: Path to the .bim file.

File.Fam: Path to the .fam file.

File.SetID: Path to the setID file.

File.SSD: Path to the .SSD file.

File.Info: Path to the .SSD.info file.

File.Kin: Path to the emmax.hBN.kinf file.

Other variables as required.

Run the script in an R environment using the source("RCodes/MainCode.R").

The script will perform SKAT analysis and generate results. QQ plots will be generated and saved as PNG images.

### Python Script
The Python script (Senior-Year-Project/PythonCodes/OBGenes_Excel_AvgVal_LogVal_BarPlot_1.py) searches for specific texts in Excel files and generates results. 

Follow the steps below to run the script:

Install the required dependencies by running pip install pandas xlsxwriter matplotlib numpy tqdm.

Open the script file and modify the following variables:

texts_to_search: List of texts to search in the Excel files.

directory_path: Path to the directory containing the Excel files.

Run the script using python Senior-Year-Project/PythonCodes/OBGenes_Excel_AvgVal_LogVal_BarPlot_1.py.

The script will search for the specified texts in the Excel files within the provided directory. It will generate results and save them to a text file (results.txt) and an Excel file (results.xlsx). Bar plots will also be generated and saved as images.

### For further questions, contact the writer at emrecangunaydin@gmail.com or emrecangunaydin@hotmail.com

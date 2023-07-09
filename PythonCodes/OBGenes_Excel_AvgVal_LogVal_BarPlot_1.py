import os
import pandas as pd
import xlsxwriter
import matplotlib.pyplot as plt
import math
import numpy as np
from tqdm import tqdm  # Import tqdm for the progress bar

# Example usage, to search for texts, change the list
texts_to_search = ["LEP", "PCSK1", "LEPR", "POMC", "MC4R", "SIM1", "NTRK2", "SH2B1", "KSR2", "ADCY3"]
directory_path = "Path/to/Files"


def search_excel_files(texts, directory):
    results = {}

    # Iterate over all files in the directory
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path) and (filename.endswith(".xlsx") or filename.endswith(".xls")):
            try:
                df = pd.read_excel(file_path)
            except Exception as e:
                print(f"Error reading {filename}: {str(e)}")
                continue

            if "SetID" in df.columns and "P.value" in df.columns:
                setid_col = df["SetID"]
                pvalue_col = df["P.value"]

                for text in texts:
                    # Find the text in the SetID column
                    mask = setid_col.eq(text)
                    if mask.any():
                        # Get the row indices where the text is found
                        row_indices = mask[mask].index
                        # Get the corresponding P.values from the P.value column
                        pvalues = pvalue_col[row_indices]
                        # Append the results to the dictionary
                        if filename not in results:
                            results[filename] = []
                        results[filename].append((text, pvalues))

    return results


results = search_excel_files(texts_to_search, directory_path)

# Classify the texts and group the results based on file names
classification = {}
for filename, text_results in results.items():
    section = None
    prevalence = None

    if "Burden" in filename:
        section = "Burden"
        if "25" in filename:
            prevalence = 25
        elif "40" in filename:
            prevalence = 40
        elif "50" in filename:
            prevalence = 50
        else:
            continue

    elif "SKAT" in filename and "SKATO" not in filename:
        section = "SKAT"
        if "25" in filename:
            prevalence = 25
        elif "40" in filename:
            prevalence = 40
        elif "50" in filename:
            prevalence = 50
        else:
            continue

    elif "SKATO" in filename:
        section = "SKATO"
        if "25" in filename:
            prevalence = 25
        elif "40" in filename:
            prevalence = 40
        elif "50" in filename:
            prevalence = 50
        else:
            continue

    if section:
        if section not in classification:
            classification[section] = {}

        if prevalence not in classification[section]:
            classification[section][prevalence] = {}

        for text, pvalues in text_results:
            if text not in classification[section][prevalence]:
                classification[section][prevalence][text] = []

            classification[section][prevalence][text].extend(pvalues)

for filename, text_results in results.items():
    for text, pvalues in text_results:
        if len(pvalues) >= 2:
            # Remove duplicate and triplicate values
            pvalues = list(set(pvalues))
            # Calculate average of the values
            average = sum(pvalues) / len(pvalues)
            # Update the results with the average value
            results[filename].append((text, [average]))

# Save results to a text file
text_file_path = os.path.join(os.path.expanduser("~"), "Writing", "Path", "results.txt")
with open(text_file_path, "w") as file:
    for section, prevalence_results in classification.items():
        file.write(f"{section} Section:\n")
        for prevalence, text_results in prevalence_results.items():
            file.write(f"Prevalence: {prevalence}\n")
            for text, pvalues in text_results.items():
                average = sum(pvalues) / len(pvalues)
                file.write(f"- Text: {text}\n")
                file.write(f"  Average Value: {average}\n")
            file.write("\n")

# Save results to an Excel file
excel_file_path = os.path.join(os.path.expanduser("~"), "Writing", "Path", "results.xlsx")
workbook = xlsxwriter.Workbook(excel_file_path)
for section, prevalence_results in tqdm(classification.items(), desc="Saving results to Excel"):
    for prevalence, text_results in prevalence_results.items():
        worksheet_name = f"{section}_{prevalence}"
        worksheet = workbook.add_worksheet(worksheet_name)
        worksheet.write_row(0, 0, ["Text", "Value"])
        row = 1
        for text, pvalues in text_results.items():
            worksheet.write(row, 0, text)
            worksheet.write_column(row, 1, pvalues)
            row += len(pvalues) + 1
workbook.close()

for section, prevalence_results in tqdm(classification.items(), desc="Generating bar plots"):
    if section == "Burden":
        bar_plot_title = "Burden"
    elif section == "SKAT":
        bar_plot_title = "SKAT"
    elif section == "SKATO":
        bar_plot_title = "SKATO"
    else:
        continue

    for prevalence, text_results in prevalence_results.items():
        plt.figure()
        plt.title(f"{bar_plot_title} (Weights.Beta: {prevalence})")  # Add prevalence to the plot title

        texts = []
        averages = []

        for text, pvalues in text_results.items():
            texts.append(text)
            average = sum(pvalues) / len(pvalues)
            averages.append(average)

        plt.bar(texts, averages, label=f"Weights.Beta = c(1, {prevalence})")

        plt.xlabel("Known Obesity Genes")
        plt.ylabel("Value")
        plt.legend()

        plt.ylim(0.0, 1.0)
        plt.yticks([i * 0.1 for i in range(11)])

        # Save the bar plot as an image with prevalence in the file name
        plot_file_path = os.path.join(os.path.expanduser("~"), "Writing", "Path", f"{bar_plot_title}_BarPlot_{prevalence}_1.png")
        plt.savefig(plot_file_path)
        plt.close()

# Convert values in the results.xlsx file to log10 and take the absolute values
converted_file_path = os.path.join(os.path.expanduser("~"), "Writing", "Path", "converted_results.xlsx")
converted_workbook = xlsxwriter.Workbook(converted_file_path)

for section, prevalence_results in classification.items():
    for prevalence, text_results in prevalence_results.items():
        converted_worksheet_name = f"converted_{section}_{prevalence}"
        converted_worksheet = converted_workbook.add_worksheet(converted_worksheet_name)
        converted_worksheet.write_row(0, 0, ["Text", "Converted Value"])
        row = 1
        for text, pvalues in text_results.items():
            converted_pvalues = [abs(math.log10(val)) for val in pvalues]
            converted_worksheet.write(row, 0, text)
            converted_worksheet.write_column(row, 1, converted_pvalues)
            row += len(pvalues) + 1

converted_workbook.close()

# Generate bar plots from the converted values
for section, prevalence_results in classification.items():
    if section == "Burden":
        bar_plot_title = "Burden"
    elif section == "SKAT":
        bar_plot_title = "SKAT"
    elif section == "SKATO":
        bar_plot_title = "SKATO"
    else:
        continue

    for prevalence, text_results in prevalence_results.items():
        plt.figure()
        plt.title(f"{bar_plot_title} (Weights.Beta: {prevalence})")  # Add prevalence to the plot title

        # Define color map for prevalence
        colors = {25: 'blue', 40: 'orange', 50: 'green'}

        texts = []
        converted_averages = []

        for text, pvalues in text_results.items():
            texts.append(text)
            converted_average = sum([abs(math.log10(val)) for val in pvalues]) / len(pvalues)
            converted_averages.append(converted_average)

        # Get the color for the current prevalence
        color = colors.get(prevalence, 'gray')

        # Plot the bar plot with the specified color
        plt.bar(texts, converted_averages, label=f"Weights.Beta = c(1, {prevalence})", color=color)

        plt.xlabel("Known Obesity Genes")
        plt.ylabel("Converted Value")
        plt.legend()

        plt.ylim(0.0, 6.0)
        plt.yticks([i * 1 for i in range(7)])

        # Save the bar plot as an image with prevalence in the file name
        plot_file_path = os.path.join(os.path.expanduser("~"), "Writing", "Path", f"{bar_plot_title}_BarPlot_{prevalence}_Converted.png")
        plt.savefig(plot_file_path)
        plt.close()

print("Results saved successfully!")

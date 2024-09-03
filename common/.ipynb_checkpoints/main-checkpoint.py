import pandas as pd
import matplotlib as plt



def plotcolumn(df,colname):
    plt.hist(df[first_column_name], bins=10, edgecolor='black')

# Add labels and title
plt.xlabel(first_column_name)
plt.ylabel('Frequency')
plt.title(f'Histogram of {first_column_name}')


def main():


    df = pd.read_csv('/Users/pf/sw/smokeDet/ai_detections.csv')

    









if __name__ == "main":
    main()
